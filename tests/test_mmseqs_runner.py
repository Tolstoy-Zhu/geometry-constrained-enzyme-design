from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from bur01_structmine.mmseqs_runner import (
    MMSEQS_FIELDS,
    build_easy_search_args,
    get_mmseqs_version,
    run_easy_cluster,
    run_easy_search,
)


def test_hong_preset_command_is_explicit() -> None:
    args = build_easy_search_args(
        "query.fasta",
        "nr",
        "hits.tsv",
        "tmp",
        threads=32,
        iterations=2,
        max_seqs=2500,
        search_type=1,
        evalue=1e-5,
    )

    assert args[:5] == ["easy-search", "query.fasta", "nr", "hits.tsv", "tmp"]
    assert args[args.index("--num-iterations") + 1] == "2"
    assert args[args.index("--max-seqs") + 1] == "2500"
    assert args[args.index("--search-type") + 1] == "1"
    assert args[args.index("--format-output") + 1] == ",".join(MMSEQS_FIELDS)


def test_invalid_mmseqs_parameters_fail_before_execution() -> None:
    with pytest.raises(ValueError, match="positive"):
        build_easy_search_args(
            "query",
            "target",
            "output",
            "tmp",
            threads=0,
            iterations=2,
            max_seqs=2500,
            search_type=1,
            evalue=1e-5,
        )


def test_fake_mmseqs_search_cluster_and_version(
    tmp_path: Path,
    fake_mmseqs_executable: Path,
) -> None:
    query = tmp_path / "query.fasta"
    target = tmp_path / "target.fasta"
    query.write_text(">query\nMSTAG\n", encoding="ascii")
    target.write_text(">WP_000001.1 [Test organism]\nMSTAG\n", encoding="ascii")
    output = tmp_path / "hits.tsv"

    run_easy_search(
        query,
        target,
        output,
        tmp_path / "search_tmp",
        executable=str(fake_mmseqs_executable),
        backend="native",
    )
    representatives, clusters = run_easy_cluster(
        target,
        tmp_path / "clusters" / "nr_s50",
        tmp_path / "cluster_tmp",
        executable=str(fake_mmseqs_executable),
        backend="native",
    )

    assert output.is_file() and output.read_text(encoding="utf-8").strip()
    assert representatives.is_file()
    assert clusters.read_text(encoding="utf-8").startswith("WP_000001.1")
    assert (
        get_mmseqs_version(
            str(fake_mmseqs_executable),
            backend="native",
        )
        == "18.8cc5c"
    )
