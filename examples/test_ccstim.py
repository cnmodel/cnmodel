"""
Test the ccstim generator.

Usage:   python test_ccstim.py

This script runs ccstim for each of it's potential inputs, and plots the resulting waveforms.
"""

import sys
import numpy as np
from cnmodel.util import ccstim
import pyqtgraph as pg

pulsetypes = ['square', 'hyp', 'timedSpikes', 'exp']

def test_cc_stim():
    """
        stim: a dictionary with keys [required]
            delay (delay to start of pulse train, msec [all]
            duration: duration of pulses in train, msec [all]
            Sfreq: stimulus train frequency (Hz) [timedSpikes]
            PT: post-train test delay [all]
            NP: number of pulses in the train [timedSpikes, exp]
            amp: amplitude of the pulses in the train [all]
            hypamp: amplitude of prehyperpolarizing pulse [hyp]
            hypdur: duration of prehyperpolarizing pulse [hyp]
            spikeTimes" times of spikes [timedSpikes]
    """
    stim = {'square': {'delay': 10, 'duration': 100, 'Sfreq': 10, 'PT': 0, 'NP': 1, 'amp': 100.},
            'hyp': {'delay': 10, 'duration': 100, 'Sfreq': 10, 'PT': 0, 'NP': 1, 'amp': 100.,
                'hypamp': -50, 'hypdur': 50},
            'timedSpikes': {'delay': 10, 'duration': 1, 'PT': 0, 'amp': 100.,
                'spikeTimes': [10., 20, 30, 40., 50.]},
            'exp': {'delay': 10, 'duration': 3, 'Sfreq': 20, 'PT': 0, 'NP': 4, 'amp': 100.}
            }

    dt = 0.1

    for p in pulsetypes:
        (w, tmax, ts) = ccstim.ccstim(stim[p], dt, pulsetype=p)
        tb = np.arange(0, w.shape[0]*dt, dt)
        pg.plot(tb, w, title=p)

if __name__ == "__main__":
    test_cc_stim()
    if sys.flags.interactive == 0 :
        pg.QtGui.QApplication.exec_()
