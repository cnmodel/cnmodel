import pytest
import numpy as np
from cnmodel.util.matlab_proc import MatlabProcess


def test_matlab():
    global proc
    try:
        proc = MatlabProcess()
    except RuntimeError:
        # no matlab available; skip this test
        pytest.skip("MATLAB unavailable")

    base_vcount = proc.who().shape[0]
    
    e4 = proc.eye(4)
    assert isinstance(e4, np.ndarray)
    assert np.all(e4 == np.eye(4))
    assert proc.who().shape[0] == base_vcount

    o6_ref = proc.ones(6, _transfer=False)
    o6 = o6_ref.get()
    assert np.all(o6 == np.ones(6))
    assert proc.who().shape[0] == base_vcount + 1
    
    del o6_ref
    assert proc.who().shape[0] == base_vcount

    proc.close()

    
if __name__ == '__main__':
    test_matlab()
