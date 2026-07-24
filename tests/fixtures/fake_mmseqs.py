#!/usr/bin/env python3
"""Small deterministic MMseqs2 stand-in for CI workflow tests."""
from __future__ import annotations

import sys
from pathlib import Path


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name = ""
    sequence: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name:
                records.append((name, "".join(sequence).upper()))
            name = line[1:]
            sequence = []
        else:
            sequence.append(line)
    if name:
        records.append((name, "".join(sequence).upper()))
    return records


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n{sequence}\n")


def option(name: str, default: str = "") -> str:
    try:
        return sys.argv[sys.argv.index(name) + 1]
    except (ValueError, IndexError):
        return default


def easy_search() -> None:
    query_path = Path(sys.argv[2])
    target_path = Path(sys.argv[3])
    output_path = Path(sys.argv[4])
    fields = option("--format-output").split(",")
    queries = read_fasta(query_path)
    targets = read_fasta(target_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows: list[str] = []
    for query_name, query_sequence in queries:
        for target_header, target_sequence in targets:
            target_name = target_header.split(maxsplit=1)[0]
            length = min(len(query_sequence), len(target_sequence))
            qaln = query_sequence[:length]
            taln = target_sequence[:length]
            identities = sum(left == right for left, right in zip(qaln, taln))
            fident = identities / length if length else 0.0
            pident = fident * 100.0
            values = {
                "query": query_name.split(maxsplit=1)[0],
                "target": target_name,
                "pident": f"{pident:.4f}",
                "fident": f"{fident:.6f}",
                "nident": str(identities),
                "alnlen": str(length),
                "qlen": str(len(query_sequence)),
                "tlen": str(len(target_sequence)),
                "qstart": "1",
                "qend": str(length),
                "tstart": "1",
                "tend": str(length),
                "qcov": f"{length / len(query_sequence):.6f}",
                "tcov": f"{length / len(target_sequence):.6f}",
                "bits": f"{20.0 + 2.0 * identities:.2f}",
                "evalue": "1e-80" if fident >= 0.50 else "1e-2",
                "qaln": qaln,
                "taln": taln,
                "tseq": target_sequence,
                "theader": target_header,
            }
            rows.append("\t".join(values.get(field, "") for field in fields))
    output_path.write_text("\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")


def easy_cluster() -> None:
    input_path = Path(sys.argv[2])
    prefix = Path(sys.argv[3])
    records = read_fasta(input_path)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    representatives: list[tuple[str, str]] = []
    representative_by_sequence: dict[str, str] = {}
    cluster_rows: list[str] = []
    for name, sequence in records:
        representative = representative_by_sequence.get(sequence)
        if representative is None:
            representative = name.split(maxsplit=1)[0]
            representative_by_sequence[sequence] = representative
            representatives.append((representative, sequence))
        cluster_rows.append(f"{representative}\t{name.split(maxsplit=1)[0]}")
    write_fasta(representatives, Path(f"{prefix}_rep_seq.fasta"))
    Path(f"{prefix}_cluster.tsv").write_text(
        "\n".join(cluster_rows) + ("\n" if cluster_rows else ""),
        encoding="utf-8",
    )


def databases() -> None:
    target = Path(sys.argv[3])
    target.parent.mkdir(parents=True, exist_ok=True)
    Path(f"{target}.dbtype").write_bytes(b"\x00\x00\x00\x00")


def main() -> None:
    if len(sys.argv) == 1 or sys.argv[1] == "version":
        print("18.8cc5c")
    elif sys.argv[1] == "easy-search":
        easy_search()
    elif sys.argv[1] == "easy-cluster":
        easy_cluster()
    elif sys.argv[1] == "databases":
        databases()
    else:
        raise SystemExit(f"Unsupported fake MMseqs2 command: {sys.argv[1]}")


if __name__ == "__main__":
    main()
