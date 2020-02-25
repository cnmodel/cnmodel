"""
Test sounds and plot waveforms.

This script tests the sound waveform generator for a variety of sounds.
It is self-contained - computes sound waveforms, and does nothing
else. 

"""
import sys
import numpy as np
import pyqtgraph as pg
from cnmodel.util import sound
from collections import OrderedDict
from scipy import signal


class test_sounds:
    def __init__(self):
        cf = 2e3
        Fs = 100e3  # sample frequency
        level = 80.0
        seed = 34978
        fmod = 10.0
        dmod = 100.0
        win = pg.GraphicsWindow()
        pipwin = win.addPlot(title="Sound Pip", row=0, col=0)
        pipmodwin = win.addPlot(title="100 % SAM modulated pip", row=1, col=0)
        noisewin = win.addPlot(title="Wideband noise", row=2, col=0)
        noisemodwin = win.addPlot(title="100 % SAM Modulated WB Noise", row=3, col=0)
        clickwin = win.addPlot(title="Clicks", row=4, col=0)
        wavewin = win.addPlot(title="Wavefile (geese)", row=5, col=0)

        pipwins = win.addPlot(title="Sound Pip Spectrum", row=0, col=1)
        pipmodwins = win.addPlot(title="100 % SAM modulated pip spectrum", row=1, col=1)
        noisewins = win.addPlot(title="Wideband noise spectrum", row=2, col=1)
        noisemodwins = win.addPlot(title="100 % SAM Modulated WB Noise", row=3, col=1)
        clickwins = win.addPlot(title="Click spectrum", row=4, col=1)
        wavewins = win.addPlot(title="Wavefile (geese) spectrum", row=5, col=1)

        stims = OrderedDict(
            [
                ("pip", (pipwin, sound.TonePip)),
                ("pipmod", (pipmodwin, sound.SAMTone)),
                ("noise", (noisewin, sound.NoisePip)),
                ("noisemod", (noisemodwin, sound.SAMNoise)),
                ("clicks", (clickwin, sound.ClickTrain)),
                ("wavefile", (wavewin, sound.ReadWavefile)),
            ]
        )

        specs = OrderedDict(
            [
                ("pip", (pipwins, sound.TonePip)),
                ("pipmod", (pipmodwins, sound.SAMTone)),
                ("noise", (noisewins, sound.NoisePip)),
                ("noisemod", (noisemodwins, sound.SAMNoise)),
                ("clicks", (clickwins, sound.ClickTrain)),
                ("wavefile", (wavewins, sound.ReadWavefile)),
            ]
        )

        for stim in stims:
            print(f"Computing stimulus: {stim:>12s}")
            if stim in ["clicks"]:
                wave = stims[stim][1](
                    rate=Fs,
                    duration=1.0,
                    dbspl=level,
                    click_duration=1e-4,
                    click_starts=1e-3 * np.linspace(10, 500, 50),
                )
                # wave = stims[stim][1](rate=Fs, dbspl=level, click_interval=10., nclicks=10,
                #                  click_duration=1e-4, click_start=10.)
            elif stim in ["wavefile"]:
                wave = stims[stim][1](
                    wavefile="examples/stim172_geese.wav", rate=Fs, dbspl=level
                )
            else:
                wave = stims[stim][1](
                    rate=Fs,
                    duration=1.0,
                    f0=cf,
                    dbspl=level,
                    pip_duration=0.8,
                    pip_start=[10e-3],
                    ramp_duration=2.5e-3,
                    fmod=fmod,
                    dmod=dmod,
                    seed=seed,
                )

            wave.sound  # force generation of the waveform
            stims[stim][0].plot(wave.time, wave.sound)
            f, Pxx_spec = signal.periodogram(
                wave.sound, Fs
            )  # window='flattop', nperseg=8192, noverlap=512, scaling='spectrum')
            specs[stim][0].plot(f, np.sqrt(Pxx_spec))
            stims[stim][0].show()
                # print(dir(stims[stim][0]))
        if sys.flags.interactive == 0:
            pg.QtGui.QApplication.exec_()

if __name__ == "__main__":
    test_sounds()
