"""
Test using sound stimulation to generate SGC spike trains.

This script uses an_model.get_spiketrain(), which internally calls MATLAB to 
generate spike trains and caches the output. A higher-level approach is to use
DummySGC, which will automatically feed the spike train into a synapse for 
input to downstream neurons (see test_sgc_input.py). Lower-level access to the
auditory nerve model is demonstrated in test_an_model.py.
"""
import numpy as np
import pyqtgraph as pg
from cnmodel import an_model
from cnmodel.util import sound

cf = 1.5e3
levels = range(-10, 101, 10)
seed = 34978

win = pg.GraphicsWindow()
p1 = win.addPlot(title='Rate-level function')

for sr in 1,2,3:
    spikes = []
    for level in levels:
        print "Level:", level
        stim = sound.TonePip(rate=100e5, duration=0.5, f0=cf, dbspl=level, 
                             pip_duration=0.5, pip_start=[0], ramp_duration=2.5e-3)
        spikes.append(an_model.get_spiketrain(cf=cf, sr=sr, seed=seed, stim=stim))
        seed += 1

    p1.plot(levels, [s.size for s in spikes], pen=(sr, 6))


import sys
if sys.flags.interactive == 0:
    pg.QtGui.QApplication.exec_()
