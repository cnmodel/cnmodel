import numpy as np
from neuron import h
from cnmodel.util import reset


def test_max_open_probability():
    reset()
    sec = h.Section()
    
    # Create AMPA and NMDA mechanisms
    # AMPA uses mode=0; no rectification
    apsd = h.AMPATRUSSELL(0.5, sec=sec)
    # For NMDA we will hold the cell at +40 mV
    npsd = h.NMDA_Kampa(0.5, sec=sec)
    
    # And a presynaptic terminal to provide XMTR input
    term = h.MultiSiteSynapse(0.5, sec=sec)
    term.nZones = 1
    term.setpointer(term._ref_XMTR[0], 'XMTR', apsd)
    term.setpointer(term._ref_XMTR[0], 'XMTR', npsd)
    
    h.celsius = 34.0
    h.finitialize()
    op = [[], []]
    for i in range(100):
        # force very high transmitter concentration for every timestep
        term.XMTR[0] = 10000
        sec.v = 40.0
        h.fadvance()
        op[0].append(apsd.Open)
        op[1].append(npsd.Open)
        
    assert np.allclose(max(op[0]), apsd.MaxOpen)
    assert np.allclose(max(op[1]), npsd.MaxOpen)
    
    
if __name__ == '__main__':
    test_max_open_probability()
