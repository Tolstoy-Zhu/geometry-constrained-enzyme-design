"""Auditable NCBI remote BLAST access for small NR snapshots and novelty checks."""
from __future__ import annotations

import hashlib
import json
import re
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable

import requests
from Bio.Blast import NCBIXML
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .io_utils import ensure_dir
from .prediction_sequence import FastaRecord
from .results_parser import parse_ncbi_protein_accession


BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
EUTILS_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


class NcbiConsentRequired(PermissionError):
    pass


@dataclass
class NcbiBlastHit:
    query: str
    target: str
    target_accession: str
    target_title: str
    pident: float
    fident: float
    nident: int
    alnlen: int
    qlen: int
    tlen: int | None
    qstart: int
    qend: int
    tstart: int
    tend: int
    qcov: float
    tcov: float | None
    bits: float
    evalue: float
    qaln: str
    taln: str
    target_sequence: str = ""
    organism: str = ""
    taxid: str = ""
    annotation_error: str = ""


def _valid_contact_email(value: str) -> bool:
    return bool(re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", value.strip()))


def _parse_submit_response(text: str) -> tuple[str, int]:
    rid_match = re.search(r"^\s*RID\s*=\s*(\S+)", text, flags=re.MULTILINE)
    rtoe_match = re.search(r"^\s*RTOE\s*=\s*(\d+)", text, flags=re.MULTILINE)
    if not rid_match:
        message = re.sub(r"\s+", " ", text).strip()[:500]
        raise RuntimeError(f"NCBI BLAST submission did not return an RID: {message}")
    return rid_match.group(1), int(rtoe_match.group(1)) if rtoe_match else 10


def _parse_fasta_records(text: str) -> dict[str, tuple[str, str, str]]:
    records: dict[str, tuple[str, str, str]] = {}
    header = ""
    sequence: list[str] = []

    def commit() -> None:
        if not header:
            return
        clean = "".join(sequence).replace(" ", "").upper()
        first_header = header.split("\x01", 1)[0]
        accession = parse_ncbi_protein_accession(first_header)
        if not accession:
            token = first_header.split(maxsplit=1)[0]
            accession = token.split("|")[-1] if "|" in token else token
        organism_match = re.search(r"\[([^\[\]]+)\]\s*$", first_header)
        organism = organism_match.group(1) if organism_match else ""
        title = re.sub(r"\s+\[[^\[\]]+\]\s*$", "", first_header).strip()
        records[accession.upper()] = (clean, title, organism)
        records[accession.split(".", 1)[0].upper()] = (clean, title, organism)

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            commit()
            header = line[1:]
            sequence = []
        else:
            sequence.append(line)
    commit()
    return records


class NcbiRemoteBlastClient:
    """Serial NCBI BLAST URL-API client with disk cache and guideline pacing."""

    def __init__(
        self,
        email: str,
        cache_dir: str | Path,
        *,
        poll_interval_seconds: float = 60.0,
        timeout_seconds: int = 7200,
        request_timeout_seconds: int = 120,
        session: requests.Session | None = None,
        sleep: Callable[[float], None] = time.sleep,
    ) -> None:
        if not _valid_contact_email(email):
            raise ValueError("A valid contact email is required by the NCBI BLAST API")
        self.email = email.strip()
        self.cache_dir = ensure_dir(cache_dir)
        self.poll_interval_seconds = max(0.0, poll_interval_seconds)
        self.timeout_seconds = timeout_seconds
        self.request_timeout_seconds = request_timeout_seconds
        self.sleep = sleep
        self.session = session or requests.Session()
        if session is None:
            retry = Retry(
                total=4,
                connect=4,
                read=4,
                status=4,
                backoff_factor=1.0,
                status_forcelist=(429, 500, 502, 503, 504),
                allowed_methods=frozenset({"GET", "POST"}),
                respect_retry_after_header=True,
            )
            adapter = HTTPAdapter(max_retries=retry)
            self.session.mount("https://", adapter)
            self.session.mount("http://", adapter)
        self.session.headers.update(
            {
                "User-Agent": (
                    "AquaEsterNet-StructMine/0.4.0 "
                    f"(NCBI BLAST research workflow; contact={self.email})"
                )
            }
        )

    def _cache_key(self, record: FastaRecord, max_hits: int, evalue: float) -> str:
        material = json.dumps(
            {
                "program": "blastp",
                "database": "nr",
                "sequence_sha256": record.sha256,
                "max_hits": max_hits,
                "evalue": evalue,
            },
            sort_keys=True,
        )
        return hashlib.sha256(material.encode("ascii")).hexdigest()

    def _write_json(self, path: Path, payload: object) -> None:
        tmp = path.with_suffix(path.suffix + ".tmp")
        tmp.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
        tmp.replace(path)

    def search(
        self,
        record: FastaRecord,
        *,
        max_hits: int = 2500,
        evalue: float = 1e-5,
        confirm_remote: bool = False,
    ) -> tuple[list[NcbiBlastHit], dict[str, object]]:
        if not confirm_remote:
            raise NcbiConsentRequired(
                "Explicit confirmation is required before sending a sequence to NCBI BLAST"
            )
        if not 1 <= max_hits <= 5000:
            raise ValueError("max_hits must be between 1 and 5000")

        key = self._cache_key(record, max_hits, evalue)
        hits_path = self.cache_dir / f"{key}.hits.json"
        provenance_path = self.cache_dir / f"{key}.provenance.json"
        if hits_path.is_file() and provenance_path.is_file():
            hits = [NcbiBlastHit(**item) for item in json.loads(hits_path.read_text(encoding="utf-8"))]
            provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
            provenance["cache_hit"] = True
            return hits, provenance

        submit = self.session.post(
            BLAST_URL,
            data={
                "CMD": "Put",
                "PROGRAM": "blastp",
                "DATABASE": "nr",
                "QUERY": record.to_fasta(),
                "EXPECT": str(evalue),
                "HITLIST_SIZE": str(max_hits),
                "FORMAT_TYPE": "XML",
                "TOOL": "AquaEsterNet-StructMine",
                "EMAIL": self.email,
            },
            timeout=self.request_timeout_seconds,
        )
        submit.raise_for_status()
        rid, estimated_seconds = _parse_submit_response(submit.text)
        if self.poll_interval_seconds:
            self.sleep(max(self.poll_interval_seconds, float(estimated_seconds)))

        deadline = time.monotonic() + self.timeout_seconds
        has_hits = False
        while True:
            status_response = self.session.get(
                BLAST_URL,
                params={
                    "CMD": "Get",
                    "FORMAT_OBJECT": "SearchInfo",
                    "RID": rid,
                    "TOOL": "AquaEsterNet-StructMine",
                    "EMAIL": self.email,
                },
                timeout=self.request_timeout_seconds,
            )
            status_response.raise_for_status()
            status_text = status_response.text
            status_match = re.search(r"Status=(\w+)", status_text)
            status = status_match.group(1).upper() if status_match else ""
            if status == "READY":
                has_hits = "ThereAreHits=yes" in status_text
                break
            if status in {"FAILED", "UNKNOWN"}:
                raise RuntimeError(f"NCBI BLAST RID {rid} ended with status {status}")
            if time.monotonic() >= deadline:
                raise TimeoutError(f"NCBI BLAST RID {rid} exceeded {self.timeout_seconds} seconds")
            self.sleep(self.poll_interval_seconds)

        hits: list[NcbiBlastHit] = []
        xml_path = self.cache_dir / f"{key}.xml"
        if has_hits:
            result_response = self.session.get(
                BLAST_URL,
                params={
                    "CMD": "Get",
                    "FORMAT_TYPE": "XML",
                    "RID": rid,
                    "TOOL": "AquaEsterNet-StructMine",
                    "EMAIL": self.email,
                },
                timeout=self.request_timeout_seconds,
            )
            result_response.raise_for_status()
            xml_path.write_text(result_response.text, encoding="utf-8")
            hits = self._parse_xml(result_response.text, record)

        provenance: dict[str, object] = {
            "backend": "NCBI BLAST URL API",
            "program": "blastp",
            "database": "nr",
            "rid": rid,
            "query_name": record.name,
            "query_sha256": record.sha256,
            "max_hits": max_hits,
            "evalue": evalue,
            "contact_email_recorded": True,
            "snapshot_date": datetime.now(timezone.utc).date().isoformat(),
            "retrieved_at": datetime.now(timezone.utc).isoformat(),
            "cache_hit": False,
            "xml_path": str(xml_path.resolve()) if xml_path.is_file() else "",
        }
        self._write_json(hits_path, [asdict(hit) for hit in hits])
        self._write_json(provenance_path, provenance)
        return hits, provenance

    @staticmethod
    def _parse_xml(xml_text: str, record: FastaRecord) -> list[NcbiBlastHit]:
        from io import StringIO

        blast_record = NCBIXML.read(StringIO(xml_text))
        hits: list[NcbiBlastHit] = []
        for alignment in blast_record.alignments:
            if not alignment.hsps:
                continue
            hsp = alignment.hsps[0]
            accession = str(alignment.accession or "").strip()
            target = str(alignment.hit_id or accession).strip()
            if not accession:
                accession = parse_ncbi_protein_accession(target) or target
            pident = 100.0 * hsp.identities / hsp.align_length if hsp.align_length else 0.0
            qcov = (hsp.query_end - hsp.query_start + 1) / len(record.sequence)
            target_title = str(alignment.hit_def or "")
            organism_match = re.search(r"\[([^\[\]]+)\]\s*$", target_title)
            hits.append(
                NcbiBlastHit(
                    query=record.name,
                    target=target,
                    target_accession=accession,
                    target_title=target_title,
                    pident=pident,
                    fident=pident / 100.0,
                    nident=int(hsp.identities),
                    alnlen=int(hsp.align_length),
                    qlen=len(record.sequence),
                    tlen=None,
                    qstart=int(hsp.query_start),
                    qend=int(hsp.query_end),
                    tstart=int(hsp.sbjct_start),
                    tend=int(hsp.sbjct_end),
                    qcov=qcov,
                    tcov=None,
                    bits=float(hsp.bits),
                    evalue=float(hsp.expect),
                    qaln=str(hsp.query),
                    taln=str(hsp.sbjct),
                    organism=organism_match.group(1) if organism_match else "",
                )
            )
        return hits

    def fetch_target_sequences(
        self,
        hits: list[NcbiBlastHit],
        *,
        api_key: str = "",
        batch_size: int = 200,
    ) -> list[NcbiBlastHit]:
        accessions = list(
            dict.fromkeys(hit.target_accession for hit in hits if hit.target_accession)
        )
        fetched: dict[str, tuple[str, str, str]] = {}
        delay = 0.11 if api_key else 0.34
        for start in range(0, len(accessions), batch_size):
            batch = accessions[start : start + batch_size]
            response = self.session.get(
                EUTILS_EFETCH_URL,
                params={
                    "db": "protein",
                    "id": ",".join(batch),
                    "rettype": "fasta",
                    "retmode": "text",
                    "tool": "AquaEsterNet-StructMine",
                    "email": self.email,
                    **({"api_key": api_key} if api_key else {}),
                },
                timeout=self.request_timeout_seconds,
            )
            response.raise_for_status()
            fetched.update(_parse_fasta_records(response.text))
            if start + batch_size < len(accessions):
                self.sleep(delay)

        for hit in hits:
            key = hit.target_accession.upper()
            record = fetched.get(key) or fetched.get(key.split(".", 1)[0])
            if record is None:
                hit.annotation_error = "NCBI efetch did not return the target sequence"
                continue
            sequence, title, organism = record
            hit.target_sequence = sequence
            hit.tlen = len(sequence)
            hit.tcov = hit.alnlen / len(sequence) if sequence else None
            hit.target_title = title or hit.target_title
            hit.organism = organism
        return hits


def hits_to_rows(hits: Iterable[NcbiBlastHit]) -> list[dict[str, object]]:
    return [asdict(hit) for hit in hits]
