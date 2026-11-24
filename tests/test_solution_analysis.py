import numpy as np
import pytest
from bur01_structmine.unified.solution_analysis import periodic_distances,largest_cluster,solution_status


def test_contacts_cross_periodic_faces_and_whole_box_translations():
    first=np.array([[.2,0,0],[3,3,3]])
    second=np.array([[9.8,0,0]])
    result=periodic_distances(first,second,[10,10,10])
    assert result[0,0]==pytest.approx(.4)
    assert np.allclose(result,periodic_distances(first+[10,-20,30],second,[10,10,10]))
    with pytest.raises(ValueError):periodic_distances(first,second,[0,10,10])


def test_cluster_connections_are_transitive_without_double_counting():
    assert largest_cluster(5,[(0,1),(1,2),(2,0),(3,4)])==3
    assert largest_cluster(0,[])==0
    assert largest_cluster(5,[])==1


def test_stable_density_cannot_override_failed_reference_comparison():
    report={'relative_density_drift':.0001,'density_comparison':{'within_tolerance':False}}
    result=solution_status(report)
    assert result['status']=='blocked' and result['blocking_issues']==['density_reference_mismatch']
    report['density_comparison']['within_tolerance']=True
    assert solution_status(report)['status']=='passed'
