"""
Test the embedded auditory nerve model with a set of tone pips.
(Zilany et al. 2014; requires MATLAB)

Adapted from Manis (makeANF_CF_RI.m)
"""
import time
import numpy as np
import pyqtgraph as pg
from cnmodel import an_model
from cnmodel.util import sound

# model fiber parameters
CF    = 1.5e3   # CF in Hz   
cohc  = 1.0    # normal ohc function
cihc  = 1.0    # normal ihc function
species = 1    # 1 for cat (2 for human with Shera et al. tuning 3 for human with Glasberg & Moore tuning)
noiseType = 1  # 1 for variable fGn (0 for fixed fGn)
fiberType = 3  # spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects "1" = Low "2" = Medium "3" = High
implnt = 0     # "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse

# stimulus parameters
F0 = CF     # stimulus frequency in Hz
Fs = 100e3  # sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = 150e-3  # stimulus duration in seconds
pdur = 0.02  # pip duration
pstart = [0.01, 0.035]  # pip start times
rt = 2.5e-3 # rise/fall time in seconds
stimdb = 65 # stimulus intensity in dB SPL

# PSTH parameters
nrep = 1000            # number of stimulus repetitions (e.g., 50)
psthbinwidth = 0.5e-3 # binwidth in seconds

stim = sound.TonePip(rate=Fs, duration=T, f0=F0, dbspl=stimdb, 
                     pip_duration=pdur, pip_start=pstart, ramp_duration=rt)
t = stim.time
pin = stim.sound
db = stim.measure_dbspl(rt, T-rt)


an_model.seed_rng(34978)

start = time.time()
vihc = an_model.model_ihc(pin, CF, nrep, 1/Fs, T+1e-3, cohc, cihc, species, _transfer=False) 
print "IHC time:", time.time() - start

start = time.time()
m, v, psth = an_model.model_synapse(vihc, CF, nrep, 1/Fs, fiberType, noiseType, implnt)
print "Syn time:", time.time() - start

win = pg.GraphicsWindow()
p1 = win.addPlot(title='Input Stimulus (%0.1f dBSPL)' % db)
p1.plot(t, pin)

p2 = win.addPlot(col=0, row=1, title='IHC voltage')
p2.setXLink(p1)
vihc = vihc.get()[0]
vihc = vihc[:len(vihc) // nrep]
t = np.arange(len(vihc)) * 1e-5
p2.plot(t, vihc)

p3 = win.addPlot(col=0, row=2, title='ANF PSTH')
p3.setXLink(p2)
ds = 100
size = psth.size // ds
psth = psth[:size*ds].reshape(size, ds).sum(axis=1)
t = np.arange(len(psth)) * 1e-5 * ds
p3.plot(t, psth[:-1], stepMode=True, fillLevel=0, fillBrush='w')


import sys
if sys.flags.interactive == 0:
    pg.QtGui.QApplication.exec_()
