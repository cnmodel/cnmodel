import numpy as np
from numpy.testing import assert_raises

from cnmodel.util import stim
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
    
    assert_raises(Exception, lambda: stim.make_pulse(params))
    params['dt'] = 0.025
    
    w, maxt, times = stim.make_pulse(params)
    
    assert w.min() == 0.0
    assert w.max() == 15
    assert w.dtype == np.float64
    
    triggers = np.argwhere(np.diff(w) > 0)[:,0] + 1
    assert np.all(triggers == times)
    assert w.sum() == 15 * len(times) * int(1/h.dt)
    
    params['PT'] = 100
    w, maxt, times = stim.make_pulse(params)
    triggers = np.argwhere(np.diff(w) > 0)[:,0] + 1
    assert triggers[-1] - triggers[-2] == 100/h.dt
    