from __future__ import annotations

from pathlib import Path

import pytest
import yaml
from pydantic import ValidationError

from bur01_structmine.config import load_config, project_root
from bur01_structmine.results_parser import filter_hits, parse_pdb_chain, parse_uniprot_accession


def test_all_shipped_configs_validate() -> None:
    root = project_root()
    paths = [
        root / "config.yaml",
        *sorted(path for path in (root / "configs").glob("*.yaml") if not path.name.endswith(".lock.yaml")),
    ]
    for path in paths:
        cfg = load_config(path)
        assert cfg["query"]["chain"] == "A"
        assert cfg["output"]["run_root"] == "results/runs"


def test_prediction_runtime_lock_is_pinned() -> None:
    lock = yaml.safe_load((project_root() / "configs" / "prediction_runtime.lock.yaml").read_text(encoding="utf-8"))
    assert lock["project_version"] == "0.3.0"
    assert lock["packages"]["colabfold"] == "1.6.1"
    assert lock["packages"]["jax"] == "0.5.3"
    assert lock["packages"]["openmm"] == "8.5.2"


def test_invalid_config_lists_multiple_errors(tmp_path: Path) -> None:
    config = {
        "query": {"pdb_file": "missing.pdb"},
        "search_databases": [{"name": "x", "enabled": True}],
        "filter": {"min_alntmscore": 2, "max_rmsd": 0, "min_alnlen": 0},
    }
    path = tmp_path / "invalid.yaml"
    path.write_text(yaml.safe_dump(config), encoding="utf-8")
    with pytest.raises(ValidationError) as exc_info:
        load_config(path)
    errors = exc_info.value.errors()
    assert len(errors) == 3


def test_identifier_parsing_and_default_filter() -> None:
    assert parse_pdb_chain("1xkl_chain_C") == ("1XKL", "C")
    assert parse_pdb_chain("pdb|6tj2|A") == ("6TJ2", "A")
    assert parse_uniprot_accession("AF-Q6RYA0-F1-model_v6") == "Q6RYA0"

    import pandas as pd

    df = pd.DataFrame(
        [
            {"target": "positive", "alntmscore": 0.7, "rmsd": 2, "qcov": 0.8, "tcov": 0.8, "alnlen": 220, "fident": 0.3},
            {"target": "P00520", "alntmscore": 0.1, "rmsd": 6, "qcov": 0.3, "tcov": 0.1, "alnlen": 100, "fident": 0.1},
        ]
    )
    filtered = filter_hits(df)
    assert filtered["target"].tolist() == ["positive"]
