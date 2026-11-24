import pytest
from bur01_structmine.unified.ligand_parameters import acid_fractions,molecule_from_spec
from bur01_structmine.unified.ligand_parameters import normalize_rounded_charge

def test_charge_population_is_normalized_and_monotone():
    f=acid_fractions(4,[3.13,4.76,6.40])
    assert sum(f)==pytest.approx(1)
    assert f[1]>f[2]>f[0]>f[3]
    assert acid_fractions(3.13,[3.13])==pytest.approx([.5,.5])
    with pytest.raises(ValueError):acid_fractions(float('nan'),[3])

@pytest.mark.parametrize('smiles,charge',[('CCO',1),('[Zn+2]',2),('CCO.CCO',0),('[CH3]',0)])
def test_unsupported_chemistry_is_stopped(smiles,charge):
    pytest.importorskip('rdkit')
    with pytest.raises(ValueError):molecule_from_spec({'name':'test','smiles':smiles,'formal_charge':charge})

def test_citrate_isomers_remain_distinct():
    pytest.importorskip('rdkit')
    from rdkit import Chem
    central=molecule_from_spec({'name':'central','smiles':'OC(=O)CC(O)(C(=O)[O-])CC(=O)O','formal_charge':-1})
    terminal=molecule_from_spec({'name':'terminal','smiles':'[O-]C(=O)CC(O)(C(=O)O)CC(=O)O','formal_charge':-1})
    assert Chem.MolToSmiles(central)!=Chem.MolToSmiles(terminal)
    assert central.GetNumAtoms()==terminal.GetNumAtoms()

def test_only_small_charge_rounding_can_be_corrected(tmp_path):
    path=tmp_path/'ligand.mol2'
    original='@<TRIPOS>ATOM\n1 C1 0 0 0 c3 1 LIG -0.10005\n2 H1 1 0 0 hc 1 LIG 0.10000\n@<TRIPOS>BOND\n1 1 2 1\n'
    path.write_text(original)
    report=normalize_rounded_charge(path,0)
    assert path.with_name('ligand.raw.mol2').read_text()==original
    assert report['raw_charge_sum']==pytest.approx(-.00005)
    assert sum(float(x.split()[8]) for x in path.read_text().splitlines()[1:3])==pytest.approx(0,abs=1e-8)
    with pytest.raises(ValueError,match='exceeds'):normalize_rounded_charge(path,1)


def test_parameter_bundle_contains_only_checked_files_and_rejects_changes(tmp_path):
    import json,zipfile
    from bur01_structmine.unified.datasets import file_hash
    from bur01_structmine.unified.ligand_parameters import package_parameters
    folder=tmp_path/'ethanol';folder.mkdir();parameter=folder/'ligand.mol2';parameter.write_text('fixture')
    molecule={'name':'ethanol','status':'numerical_checks_passed','sha256':{'ligand.mol2':file_hash(parameter)}}
    (folder/'parameter_report.json').write_text(json.dumps(molecule))
    (tmp_path/'parameterization_report.json').write_text(json.dumps({'status':'numerical_checks_passed','molecules':[molecule]}))
    with zipfile.ZipFile(package_parameters(tmp_path)) as archive:
        assert archive.read('ethanol/ligand.mol2')==b'fixture'
        assert 'ethanol/parameter_report.json' in archive.namelist()
    parameter.write_text('changed')
    with pytest.raises(ValueError,match='checksum'):package_parameters(tmp_path)
