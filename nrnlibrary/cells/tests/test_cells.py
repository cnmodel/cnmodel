import os, pickle
import numpy as np

import nrnlibrary
import nrnlibrary.cells as cells
from nrnlibrary.util.testing import IVCurve, run_vc, run_democlamp

path = os.path.dirname(__file__)
audit = nrnlibrary.AUDIT_TESTS

def cell_info(cell):
    """
    Run I/V protocol on cell and generate a dict of results.
    """
    iv = IVCurve()
    iv.run([-1.0, 1.0, 0.1], cell)
    info = dict(
        spikes=iv.spike_times(),
        rmp=iv.rest_vm(),
        rm=iv.input_resistance(),
        vpeak=iv.peak_vm(),
        vss=iv.steady_vm(),
        )
    return info

def result_file(key):
    """
    Return a file name to be used for storing / retrieving test results
    given *key*.
    """
    return os.path.join(path, 'cell_data', key + '.pk')

def load_cell_info(key):
    """
    Load prior test results for *key*.
    If there are no prior results, return None.
    """
    fn = result_file(key)
    if os.path.isfile(fn):
        return pickle.load(open(fn, 'rb'))
    return None

def save_cell_info(info, key):
    """
    Store test results for *key*.
    """
    fn = result_file(key)
    dirname = os.path.dirname(fn)
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    pickle.dump(info, open(fn, 'wb'))
    
    
def assert_cell_info(cell, key):
    """
    Test *cell* and raise exception if the results do not match prior
    data.
    """
    info = cell_info(cell)
    expect = load_cell_info(key)
    
    if expect is not None:
        
        # Check test structures are the same
        assert len(info) == len(expect)
        for k in info:
            assert k in expect
            
        # Check data matches
        for k in info:
            assert np.all(info[k] == expect[k])
    else:
        if not audit:
            raise Exception("No prior test results for cell type '%s'. "
                            "Run test.py --audit store new test data." % key)
        
        print "New test results for %s:" % key
        print info
        print "Store? [y/n]",
        yn = raw_input()
        if yn.lower().startswith('y'):
            save_cell_info(info, key)

def test_bushy():
    cell = cells.Bushy(species='mouse')
    assert_cell_info(cell, 'bushy_mouse')
    

def test_stellate():
    cell = cells.TStellate(nav11=True)
    assert_cell_info(cell, 'tstellate_nav11')
    
