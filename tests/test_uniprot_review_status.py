import pytest
from bur01_structmine.uniprot_client import UniProtClient


@pytest.mark.parametrize("entry_type,expected", [
    ("UniProtKB reviewed (Swiss-Prot)", "True"),
    ("UniProtKB unreviewed (TrEMBL)", "False"),
    ("UNIPROTKB UNREVIEWED (TREMBL)", "False"),
    ("", ""),
    ("unknown", ""),
])
def test_review_status_does_not_match_substring(tmp_path, entry_type, expected):
    client = UniProtClient(cache_dir=tmp_path, sleep_seconds=0)
    result = client._to_record({"primaryAccession": "A0A286T7C7", "entryType": entry_type})
    assert result.is_reviewed == expected


def test_inactive_record_has_explicit_error_and_unknown_review_status(tmp_path, monkeypatch):
    from bur01_structmine.annotator import _set_uniprot
    client = UniProtClient(cache_dir=tmp_path, sleep_seconds=0)
    monkeypatch.setattr(client.http, "get_json", lambda *args, **kwargs: {
        "results": [{"primaryAccession": "A0A351EQR4", "entryType": "Inactive"}]})
    record = client.fetch_many(["A0A351EQR4"])["A0A351EQR4"]
    assert "inactive" in client.errors["A0A351EQR4"]
    output = {}
    _set_uniprot(output, [record])
    assert output["uniprot_is_reviewed"] == ""
    assert output["uniprot_protein_sequence"] == ""
