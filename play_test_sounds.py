"""
Test sounds and plot waveforms.

This script tests the sound waveform generator for a variety of sounds

"""
import numpy as np
import pyqtgraph as pg
from cnmodel.util import sound
from collections import OrderedDict
from scipy import signal
import PySounds
import sys

if len(sys.argv) >= 2:
    stim = sys.argv[1]
else:
    exit()


PS = PySounds.PySounds()
    
cf = 2e3
Fs = 44100  # sample frequency
level = 80.
seed = 34978
fmod = 5.
dmod = 10.
# win = pg.GraphicsWindow()
# pipwin = win.addPlot(title='sound pip', row=0, col=0)
# pipmodwin = win.addPlot(title='100 \% SAM modulated pip', row=1, col=0)
# noisewin = win.addPlot(title='WB noise', row=2, col=0)
# noisemodwin = win.addPlot(title='100 \% SAM Modulated WB Noise', row=3, col=0)
# clickwin = win.addPlot(title='clicks', row=4, col=0)
#
# pipwins = win.addPlot(title='sound pip Spec', row=0, col=1)
# pipmodwins = win.addPlot(title='100 \% SAM modulated pip', row=1, col=1)
# noisewins = win.addPlot(title='WB noise', row=2, col=1)
# noisemodwins = win.addPlot(title='100 \% SAM Modulated WB Noise', row=3, col=1)
# clickwins = win.addPlot(title='click spec', row=4, col=1)
pipwin = None
pipmodwin = None
noisewin = None
noisemodwin = None
clickwin = None
pipwins = None
pipmodwins = None
noisewins = None
noisemodwins = None
clickwins = None
fmwins = None

stims = OrderedDict([('pip', (pipwin, sound.TonePip)),
                     ('pipmod', (pipmodwin, sound.SAMTone)),
                     ('noise', (noisewin, sound.NoisePip)),
                     ('noisemod', (noisemodwin, sound.SAMNoise)),
                     ('clicks', (clickwin, sound.ClickTrain)),
                     ('fmsweep', (fmwins, sound.FMSweep))])
                     
specs = OrderedDict([('pip', (pipwins, sound.TonePip)),
                     ('pipmod', (pipmodwins, sound.SAMTone)),
                     ('noise', (noisewins, sound.NoisePip)),
                     ('noisemod', (noisemodwins, sound.SAMNoise)),
                     ('clicks', (clickwins, sound.ClickTrain)),
                     ('fmsweep', (fmwins, sound.FMSweep))])

if stim in stims:
    print stim
    if stim in ['clicks']:
        wave = stims[stim][1](rate=Fs, duration=1.0, dbspl=level,
                         click_duration=1e-4, click_starts=1e-3*np.linspace(10, 500, 10))
    elif stim in ['fmsweep']:
        wave = stims[stim][1](rate=Fs, duration=0.25, dbspl=level,
                            start=0., ramp='logarithmic', freqs=[3000, 16000])
    elif stim in ['pip', 'pipmod', 'noise', 'noisemod']:
        wave = stims[stim][1](rate=Fs, duration=2.0, f0=cf, dbspl=level, 
                         pip_duration=1.8, pip_start=[10e-3], ramp_duration=2.5e-3,
                         fmod=fmod, dmod=dmod, seed=seed)
    # stims[stim][0].plot(wave.time, wave.sound)
    # f, Pxx_spec = signal.periodogram(wave.sound, Fs) # window='flattop', nperseg=8192, noverlap=512, scaling='spectrum')
    # specs[stim][0].plot(f, np.sqrt(Pxx_spec))
    print ('Playing %s' % stim)
    PS.playSound(wave.sound, wave.sound, Fs)

# if sys.flags.interactive == 0:
#     pg.QtGui.QApplication.exec_()
