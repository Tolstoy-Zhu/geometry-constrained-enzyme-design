import json
import numpy as np
import pytest
from bur01_structmine.unified.protonation import assign_microstate,normalize_protein_pdb,apply_microstate
from bur01_structmine.unified.protein_system import component_counts,verify_construct
from bur01_structmine.unified.equilibration import equilibration_schedule,assess_equilibration,aligned_rmsd,write_amber_restart


def group(residue,position,pka,chain='A',kind=None,label=None):
    return {'residue':residue,'position':position,'pka':pka,'chain':chain,
            'type':kind or ('HIS' if residue=='HIS' else 'COO'),'label':label or residue,'coupled':[]}


def protocol():
    return {'temperature_K':313.15,'timestep_fs':2.,'checkpoint_ps':50.,'nvt_ps':200.,
        'restrained_npt_ps':1000.,'unrestrained_npt_ps':1000.,'restraint_kj_mol_nm2':[1000.,500.,200.,50.],
        'qc':{'temperature_mean_tolerance_K':10.,'volume_cv_max':.1,'core_rmsd_max_A':3.,'geometry_min_A':1.8,'geometry_max_A':4.5}}


def test_pka_terminal_coo_is_not_a_serine_sidechain():
    report={'groups':[group('SER',99,3.1,kind='COO',label='C-  99 A')]}
    policy=assign_microstate(report,pH=4.)
    assert policy['assignments'][0]['variant']=='charged_C_terminus'
    with pytest.raises(ValueError,match='Neutral C terminus'):assign_microstate(report,pH=2.)


def test_symmetric_sites_and_parent_conservation_are_explicit():
    report={'groups':[group('ASP',2,3.8),group('ASP',2,4.4,'B')]}
    parent=assign_microstate(report,pH=4.)
    assert parent['assignments'][0]['variant']=='ASH'
    changed={'groups':[group('ASP',2,3.7),group('ASP',2,3.8,'B')]}
    assert assign_microstate(changed,pH=4.,reference=parent)['assignments'][0]['variant']=='ASH'
    assert assign_microstate(changed,pH=4.)['assignments'][0]['variant']=='ASP'


def test_catalytic_microstate_is_a_recorded_override():
    report={'groups':[group('HIS',9,6.9)]}
    default=assign_microstate(report,pH=4.)
    assert default['assignments'][0]['variant']=='HIP'
    override={'position':9,'residue':'HIS','variant':'HID','reason':'Conditional catalytic hypothesis'}
    policy=assign_microstate(report,pH=4.,overrides=[override])
    assert policy['assignments'][0]['variant']=='HID'
    assert policy['assignments'][0]['mean_pka']==6.9
    assert policy['production_ready'] is False
    with pytest.raises(ValueError):assign_microstate(report,pH=4.,overrides=[{**override,'position':8}])
    with pytest.raises(ValueError):assign_microstate(report,pH=4.,overrides=[{**override,'residue':'ASP'}])
    with pytest.raises(ValueError):assign_microstate(report,pH=4.,overrides=[override,override])
    with pytest.raises(ValueError):assign_microstate(report,pH=4.,reference={**default,'pH':7.})


def test_atom_name_alignment_keeps_cd1_carbon_and_coordinates(tmp_path):
    source=tmp_path/'source.pdb';target=tmp_path/'normal.pdb'
    source.write_text('ATOM      1 CD1  ILE A   1       1.000   2.000   3.000  1.00  0.00           C\nEND\n')
    normalize_protein_pdb(source,target)
    line=target.read_text().splitlines()[0]
    assert line[12:16]==' CD1'
    assert line[30:54]==source.read_text().splitlines()[0][30:54]
    assert line[76:78].strip()=='C'


def test_chosen_histidine_template_survives_canonical_names(tmp_path):
    path=tmp_path/'protein.pdb'
    path.write_text('ATOM      1  HD1 HIS A   9       1.000   2.000   3.000  1.00  0.00           H\nEND\n')
    policy=assign_microstate({'groups':[group('HIS',9,6.9)]},pH=4.)
    apply_microstate(path,path,policy)
    line=path.read_text().splitlines()[0]
    assert line[17:20]=='HIP' and line[12:16].strip()=='HD1'


def test_concentrations_use_actual_box_and_reject_vanishing_components():
    assert component_counts([{'name':'ethanol','concentration_M':1.}],343)=={'ethanol':207}
    with pytest.raises(ValueError):component_counts([{'name':'trace','concentration_M':.000001}],343)
    with pytest.raises(ValueError):component_counts([{'name':'x','concentration_M':float('nan')}],343)


def test_equilibration_is_bounded_and_ends_unrestrained():
    schedule=equilibration_schedule(protocol())
    assert len(schedule)==44 and schedule[-1]['end_step']*.002==2200
    assert all(s['restraint_k']==0 for s in schedule[-20:])
    assert schedule[3]['temperature_K']==313.15
    with pytest.raises(ValueError):equilibration_schedule({**protocol(),'unrestrained_npt_ps':4000})
    with pytest.raises(ValueError):equilibration_schedule({**protocol(),'nvt_ps':201})


def test_geometry_failure_is_not_hidden_by_temperature_or_rmsd():
    rows=[{'stage':'unrestrained_npt','temperature_K':313.2,'volume_nm3':600.,'core_rmsd_A':1.,'geometry_A':{'A:site':3.1}} for _ in range(10)]
    assert assess_equilibration(rows,protocol())['status']=='passed'
    rows[-1]={**rows[-1],'geometry_A':{'A:site':4.6}}
    result=assess_equilibration(rows,protocol())
    assert result['status']=='blocked' and 'final_geometry:A:site' in result['blocking_issues']
    assert result['geometry']['A:site']['in_range_fraction']==.8
    rows[-1]['temperature_K']=float('nan')
    assert assess_equilibration(rows,protocol())['status']=='blocked'


def test_alignment_removes_rotation_not_structural_deformation():
    xyz=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,2]])
    rotation=np.array([[0,1,0],[-1,0,0],[0,0,1]])
    assert aligned_rmsd(xyz,xyz@rotation+20)<1e-12
    deformed=xyz.copy();deformed[3,2]=5
    assert aligned_rmsd(xyz,deformed)>.5


def test_amber_restart_roundtrip_with_box(tmp_path):
    app=pytest.importorskip('openmm.app');unit=pytest.importorskip('openmm.unit')
    xyz=np.array([[.5,10.,20.],[2.,3.,4.],[5.,6.,7.]])
    path=tmp_path/'state.inpcrd';write_amber_restart(path,xyz,np.diag([60.,70.,80.]))
    result=app.AmberInpcrdFile(str(path))
    assert np.allclose(result.positions.value_in_unit(unit.angstrom),xyz)
    assert np.allclose(np.array(result.boxVectors.value_in_unit(unit.angstrom)),np.diag([60,70,80]))
    with pytest.raises(ValueError):write_amber_restart(path,xyz,[[60,1,0],[0,70,0],[0,0,80]])


def test_construct_check_rejects_an_interchain_bond():
    app=pytest.importorskip('openmm.app')
    top=app.Topology();chain=top.addChain()
    first=top.addResidue('ALA',chain);second=top.addResidue('GLY',chain)
    a=top.addAtom('C',app.element.carbon,first);b=top.addAtom('N',app.element.nitrogen,second)
    verify_construct(top,{'A':'A','B':'G'})
    top.addBond(a,b)
    with pytest.raises(ValueError,match='between protein chains'):verify_construct(top,{'A':'A','B':'G'})
