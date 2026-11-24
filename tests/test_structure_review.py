from pathlib import Path
from bur01_structmine.config import project_root
from bur01_structmine.structure_search import load_ca_trace
from bur01_structmine.unified.structure_review import review_structure

def test_review_checks_actual_construct_and_specific_atoms(tmp_path):
    pdb=project_root()/'data/query/Bur01m.pdb'
    trace=load_ca_trace(pdb,chain_id='A')
    result=review_structure(pdb,'M'+trace.sequence,tmp_path,atom_rules=[{
        'positions':[100,263],'first_atoms':['OG'],'second_atoms':['NE2'],'minimum_A':.1,'maximum_A':.2}])
    assert result['status']=='blocked'
    assert 'full_construct_sequence_mismatch' in result['blocking_issues']
    assert 'protocol_specific_catalytic_geometry' in result['blocking_issues']
    assert result['production_ready'] is False


def test_assembly_rebuilds_hydrogens_for_selected_protonation(tmp_path):
    from bur01_structmine.unified.assembly_relaxation import preparation_pdb_lines
    pdb=tmp_path/'source.pdb'
    def atom(serial,name,residue,position,element):
        return f'ATOM  {serial:5d} {name:^4s} {residue:3s} A{position:4d}    {1.:8.3f}{2.:8.3f}{3.:8.3f}  1.00  0.00          {element:>2s}'
    pdb.write_text('\n'.join([atom(1,'HD1','HIS',13,'H'),atom(2,'ND1','HIS',13,'N'),
                              atom(3,'HE2','HIS',274,'H'),atom(4,'NE2','HIS',274,'N'),
                              atom(5,'H','MET',1,'H'),atom(6,'N','MET',1,'N'),'TER','END']))
    result=preparation_pdb_lines(pdb,[274])
    assert len(result)==5
    assert result[0][17:20]=='HIE' and result[1][17:20]=='HID'
    assert [line[12:16].strip() for line in result[:3]]==['ND1','NE2','N']
