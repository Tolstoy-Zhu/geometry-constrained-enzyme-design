import json
import pytest
from bur01_structmine.unified.msa import extend_msa,read_monomer_a3m

def test_extension_preserves_insertions_and_records_missing_evidence(tmp_path):
    source=tmp_path/'core.a3m';source.write_text('#4\t1\n>query\nACDE\n>homolog\nACq-E\n')
    out=tmp_path/'full.a3m'
    receipt=extend_msa(source,'ACDE','MMACDFG',2,out)
    assert read_monomer_a3m(out)==[('query','MMACDFG'),('homolog','--ACq-E-')]
    assert receipt['uncovered_terminal_columns']=={'n_terminal':2,'c_terminal':1}
    assert receipt['query_changes']==[{'position':6,'original':'E','query':'F'}]
    assert not receipt['new_public_sequence_submission']
    with pytest.raises(FileExistsError):extend_msa(source,'ACDE','MMACDFG',2,out)
    with pytest.raises(ValueError,match='original core'):extend_msa(source,'ACDF','MMACDFG',2,out)

@pytest.mark.parametrize('text',['#4,4\t1,1\n>q\nACDEACDE\n','>q\nACDE\n>h\nACD\n','>q\nACdE\n'])
def test_ambiguous_alignments_fail_closed(tmp_path,text):
    path=tmp_path/'bad.a3m';path.write_text(text)
    with pytest.raises(ValueError):read_monomer_a3m(path)
