import numpy as np

from nrnlibrary.util import testing
from neuron import h

h.dt = 0.025

def test_make_pulse():
    params = dict(
        delay=10,
        Sfreq=50,
        dur=1,
        amp=15,
        PT=0,
        NP=5,
        )
    w, maxt, times = testing.make_pulse(params)
    
    assert w.min() == 0.0
    assert w.max() == 15
    assert w.dtype == np.float64
    
    triggers = np.argwhere(np.diff(w) > 0)[:,0] + 1
    assert np.all(triggers == times)
