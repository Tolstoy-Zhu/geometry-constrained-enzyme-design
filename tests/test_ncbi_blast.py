from __future__ import annotations

from pathlib import Path

import pytest
import requests

from bur01_structmine.ncbi_blast import (
    NcbiBlastHit,
    NcbiConsentRequired,
    NcbiRemoteBlastClient,
    _parse_fasta_records,
)
from bur01_structmine.prediction_sequence import FastaRecord


class FakeResponse:
    def __init__(self, text: str) -> None:
        self.text = text
        self.status_code = 200

    def raise_for_status(self) -> None:
        return None


class FakeBlastSession(requests.Session):
    def __init__(self) -> None:
        super().__init__()
        self.posts: list[dict[str, object]] = []
        self.gets: list[dict[str, object]] = []

    def post(self, url, data=None, timeout=None, **kwargs):  # type: ignore[override]
        self.posts.append({"url": url, "data": data, "timeout": timeout})
        return FakeResponse("RID = TEST_RID_1\nRTOE = 0\n")

    def get(self, url, params=None, timeout=None, **kwargs):  # type: ignore[override]
        self.gets.append({"url": url, "params": params, "timeout": timeout})
        if params and params.get("FORMAT_OBJECT") == "SearchInfo":
            return FakeResponse("Status=READY\nThereAreHits=yes\n")
        if params and params.get("rettype") == "fasta":
            return FakeResponse(
                ">ref|WP_000001.1| hydrolase [Test bacterium]\nMSTAG\n"
            )
        return FakeResponse("<BlastOutput />")


def _hit() -> NcbiBlastHit:
    return NcbiBlastHit(
        query="Bur01",
        target="ref|WP_000001.1|",
        target_accession="WP_000001.1",
        target_title="hydrolase",
        pident=80.0,
        fident=0.8,
        nident=4,
        alnlen=5,
        qlen=5,
        tlen=None,
        qstart=1,
        qend=5,
        tstart=1,
        tend=5,
        qcov=1.0,
        tcov=None,
        bits=50.0,
        evalue=1e-20,
        qaln="MSTAG",
        taln="MSTAG",
    )


def test_ncbi_requires_contact_and_per_request_consent(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="contact email"):
        NcbiRemoteBlastClient("", tmp_path)

    client = NcbiRemoteBlastClient(
        "researcher@example.org",
        tmp_path,
        session=FakeBlastSession(),
        poll_interval_seconds=0,
    )
    with pytest.raises(NcbiConsentRequired):
        client.search(FastaRecord("query", "MSTAG"), confirm_remote=False)


def test_ncbi_search_is_serial_cached_and_auditable(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = FakeBlastSession()
    client = NcbiRemoteBlastClient(
        "researcher@example.org",
        tmp_path,
        session=session,
        poll_interval_seconds=0,
        sleep=lambda _: None,
    )
    monkeypatch.setattr(
        NcbiRemoteBlastClient,
        "_parse_xml",
        staticmethod(lambda xml, record: [_hit()]),
    )
    record = FastaRecord("Bur01", "MSTAG")

    hits, provenance = client.search(
        record,
        max_hits=2500,
        evalue=1e-5,
        confirm_remote=True,
    )
    cached_hits, cached_provenance = client.search(
        record,
        max_hits=2500,
        evalue=1e-5,
        confirm_remote=True,
    )

    assert len(session.posts) == 1
    assert session.posts[0]["data"]["EMAIL"] == "researcher@example.org"  # type: ignore[index]
    assert session.posts[0]["data"]["HITLIST_SIZE"] == "2500"  # type: ignore[index]
    assert hits[0].target_accession == cached_hits[0].target_accession
    assert provenance["database"] == "nr"
    assert cached_provenance["cache_hit"] is True


def test_ncbi_efetch_sequence_mapping_and_rate_limit(tmp_path: Path) -> None:
    session = FakeBlastSession()
    sleeps: list[float] = []
    client = NcbiRemoteBlastClient(
        "researcher@example.org",
        tmp_path,
        session=session,
        poll_interval_seconds=0,
        sleep=sleeps.append,
    )
    second = _hit()
    second.target_accession = "WP_000002.1"

    fetched = client.fetch_target_sequences([_hit(), second], batch_size=1)

    assert fetched[0].target_sequence == "MSTAG"
    assert fetched[0].organism == "Test bacterium"
    assert fetched[1].annotation_error
    assert sleeps == [0.34]


def test_fasta_parser_indexes_versioned_and_unversioned_accessions() -> None:
    parsed = _parse_fasta_records(
        ">ref|WP_012345678.1| alpha/beta hydrolase [Example species]\nMSTAG\n"
    )

    assert parsed["WP_012345678.1"][0] == "MSTAG"
    assert parsed["WP_012345678"][2] == "Example species"
