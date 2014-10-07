"""
Test the embedded auditory nerve model with a set of tone pips.
(Zilany et al. 2014; requires MATLAB)

Adapted from Manis (makeANF_CF_RI.m)
"""
import numpy as np
import pyqtgraph as pg
from nrnlibrary import an_model
from nrnlibrary.util import audio

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
T  = 50e-3  # stimulus duration in seconds
rt = 2.5e-3 # rise/fall time in seconds
stimdb = 65 # stimulus intensity in dB SPL

# PSTH parameters
nrep = 100            # number of stimulus repetitions (e.g., 50)
psthbinwidth = 0.5e-3 # binwidth in seconds

t = np.arange(0, T, 1/Fs) # time vector
mxpts = len(t)
irpts = rt * Fs  # number of samples to ramp

pin = np.sqrt(2) * 20e-6 * 10**(stimdb/20) * np.sin(2*np.pi*F0*t) # unramped stimulus
pin[:irpts] *= np.linspace(0, 1, irpts) 
pin[-irpts:] *= np.linspace(1, 0, irpts)

an_model.seed_rng(34978)
import time
start = time.time()
vihc = an_model.model_ihc(pin,CF,nrep,1/Fs,T*2,cohc,cihc,species, _transfer=False) 
print "IHC time:", time.time() - start

start = time.time()
psth = []
m = an_model.get_matlab()
for i in range(1):
    m, v, p = an_model.model_synapse(vihc,CF,nrep,1/Fs,fiberType,noiseType,implnt, _transfer=False) 
    psth.append(p.get())
psth = np.concatenate(psth, axis=0).sum(axis=0)
print "Syn time:", time.time() - start

win = pg.GraphicsWindow()
p1 = win.addPlot(title='Input Stimulus')
p1.plot(t, pin)

p2 = win.addPlot(col=0, row=1, title='IHC voltage')
vihc = vihc.get()[0]
vihc = vihc[:len(vihc) // nrep]
t = np.arange(len(vihc)) * 1e-5
p2.plot(t, vihc)

p3 = win.addPlot(col=0, row=2, title='PSTH')
ds = 100
size = psth.size // ds
psth = psth[:size*ds].reshape(size, ds).sum(axis=1)
t = np.arange(len(psth)) * 1e-5 * ds
p3.plot(t, psth)
