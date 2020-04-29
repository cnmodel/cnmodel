import pytest
import numpy as np
# from cnmodel.util.matlab_proc import MatlabProcess
try:
    import matlab.engine
    MATLAB_FOUND = True
except:
    MATLAB_FOUND = False

def test_matlab():
    global proc
    try:
        # proc = MatlabProcess()
        proc = matlab.engine.start_matlab()
    except RuntimeError:
        # no matlab available; skip this test
        pytest.skip("MATLAB unavailable")

    base_vcount = proc.who(nargout=1) # .shape[0]
    assert len(base_vcount) == 0
    
    e4 = np.array(proc.eye(4, nargout=1))

    assert isinstance(e4, np.ndarray)
    assert np.all(e4 == np.eye(4))
    # assert proc.who().shape[0] == base_vcount

    o6_ref = np.array(proc.ones(6, nargout=1)) #  _transfer=False))
    # o6 = o6_ref.get()
    o6 = np.array(proc.ones(6, nargout=1))
    assert np.all(o6 == np.ones(6))
    # print(proc.who(nargout=1))
    # assert proc.who(nargout=1) == base_vcount + 1
    
    del o6_ref
    # assert proc.who().shape[0] == base_vcount

    proc.close()

    
if __name__ == '__main__':
    if MATLAB_FOUND:
        test_matlab()
        print('passed')
