from __future__ import annotations

from unittest.mock import Mock

import pandas as pd
import requests

from bur01_structmine.annotator import annotate_hits
from bur01_structmine.http_utils import CachedSession


class JsonResponse:
    status_code = 200
    text = '{"value": 7}'

    def raise_for_status(self) -> None:
        return None

    def json(self) -> dict[str, int]:
        return {"value": 7}


def test_cached_session_reuses_disk_cache_and_has_retry_policy(tmp_path) -> None:
    session = CachedSession(tmp_path, max_age_seconds=3600)
    session.session.get = Mock(return_value=JsonResponse())
    assert session.get_json("https://example.test/item") == {"value": 7}
    assert session.get_json("https://example.test/item") == {"value": 7}
    assert session.session.get.call_count == 1
    retry = session.session.get_adapter("https://").max_retries
    assert retry.total == 4
    assert 429 in retry.status_forcelist


def test_stale_cache_is_used_after_network_failure(tmp_path) -> None:
    session = CachedSession(tmp_path, max_age_seconds=-1)
    session.session.get = Mock(return_value=JsonResponse())
    assert session.get_json("https://example.test/stale") == {"value": 7}
    session.session.get = Mock(side_effect=requests.ConnectionError("offline"))
    assert session.get_json("https://example.test/stale") == {"value": 7}
    assert "ConnectionError" in session.last_error


def test_annotation_disabled_is_explicit(tmp_path) -> None:
    df = pd.DataFrame([{"target": "AF-Q6RYA0-F1-model_v6", "annotation_source": "alphafold"}])
    out = annotate_hits(
        df,
        enable_rcsb=False,
        enable_uniprot=False,
        enable_ncbi=False,
        enable_alphafold=False,
        cache_dir=tmp_path,
    )
    assert out.loc[0, "annotation_status"] == "disabled"
    assert out.loc[0, "annotation_errors"] == ""
    assert out.loc[0, "final_protein_sequence"] == ""
