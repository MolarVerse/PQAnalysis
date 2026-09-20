from . import pytestmark

from PQAnalysis.traj import MDEngineFormat



def test_is_qmcfc_type():
    assert MDEngineFormat.is_qmcfc_type(MDEngineFormat.QMCFC) == True
    assert MDEngineFormat.is_qmcfc_type(MDEngineFormat.PQ) == False

    assert MDEngineFormat.is_qmcfc_type("qmcfc") == True
    assert MDEngineFormat.is_qmcfc_type("pq") == False
