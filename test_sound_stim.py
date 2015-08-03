"""
Test using sound stimulation to generate input to bushy cell
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
