"""
Test using sound stimulation to generate SGC spike trains.

This script uses an_model.get_spiketrain(), which internally calls MATLAB to 
generate spike trains and caches the output. A higher-level approach is to use
DummySGC, which will automatically feed the spike train into a synapse for 
input to downstream neurons (see test_sgc_input.py). Lower-level access to the
auditory nerve model is demonstrated in test_an_model.py.

This script also can measure the time required to generate the spike trains, or
to retrieve them from the cache. In addition, it can access a second interface
to the Zilany et al. mode that is in pure python ("cochlea") for comparison.
In general the speed to compute the spike trains for the same sets of stimuli
is faster with the pure Python interface than with the Matlab interface, and
fastest for retrieval of pre-computed trains. Note that changing the value
of "seed" will force a recomputation of the spike trains.
"""
import numpy as np
import pyqtgraph as pg
from cnmodel import an_model
from cnmodel.util import sound
import cochlea

import time
seed = 34986
    

def time_usage(func):
    def wrapper(*args, **kwargs):
        beg_ts = time.time()
        res = func(*args, **kwargs)
        end_ts = time.time()
        print("** Elapsed time: %f" % (end_ts - beg_ts))
        return res
    return wrapper


def set_dbspl(signal, dbspl):
    """Scale the level of `signal` to the given dB_SPL."""
    p0 = 20e-6
    rms = np.sqrt(np.sum(signal**2) / signal.size)
    scaled = signal * 10**(dbspl / 20.0) * p0 / rms
    return scaled


@time_usage
def sound_stim(seed, useMatlab=True):
    cf = 1.5e3
    levels = range(-10, 101, 10)

    result = {}
    if useMatlab:
        simulator = 'matlab'
    else:
        simulator = 'cochlea'
    for sr in 1,2,3:
        spikes = []
        for level in levels:
            stim = sound.TonePip(rate=100e3, duration=0.5, f0=cf, dbspl=level, 
                             pip_duration=0.5, pip_start=[0], ramp_duration=2.5e-3)
            if simulator == 'cochlea':
                stim._sound = set_dbspl(stim.sound, level) # adjust scaling here
            spikes.append(an_model.get_spiketrain(cf=cf, sr=sr, seed=seed, stim=stim, simulator=simulator))
            seed += 1
        result[sr] = {'levels': levels, 'spikes': spikes}
    return result

import sys

def runtest():
    usematlab = True
    if len(sys.argv) > 0:
        if len(sys.argv) == 1:
            print 'Call requires argument, must be either "matlab" or "cochlea"; default is "matlab"'
            exit()
        flag = sys.argv[1]
        if flag not in ['matlab', 'cochlea']:
            print 'Flag must be either "matlab" or "cochlea"; default is "matlab"'
            exit()
        if flag == 'cochlea':
            usematlab=False
    if usematlab:
        print 'Running with matlab simulator'
    else:
        print 'Running with MR cochlea simulator'
            
    result = sound_stim(seed, useMatlab=usematlab)

    win = pg.GraphicsWindow()
    p1 = win.addPlot(title='Rate-level function')
    for i, x in enumerate(result.keys()):
        p1.plot(result[x]['levels'], [s.size for s in result[x]['spikes']], pen=(x, 6))
    return win

if __name__ == '__main__':
    win = runtest()
    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()
