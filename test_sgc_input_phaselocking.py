import sys
import numpy as np
import pyqtgraph as pg
from neuron import h
from cnmodel.protocols import Protocol
from cnmodel import cells
from cnmodel.util import sound
import cnmodel.util.pynrnutilities as PU


class SGCInputTestPL(Protocol):
    def run(self, temp=34.0, dt=0.025, seed=575982035):
        preCell = cells.DummySGC(cf=4000, sr=2)
        postCell = cells.Bushy.create()
        synapse = preCell.connect(postCell)
        self.pre_cell = preCell
        self.post_cell = postCell
        self.synapse = synapse
        self.run_duration = 1.
        self.pip_duration = 0.8
        self.pip_start = [0.02]
        self.Fs = 100e3
        self.f0 = 4000.
        self.stim = sound.TonePip(rate=self.Fs, duration=self.run_duration, f0=self.f0, dbspl=60,
                                  ramp_duration=2.5e-3, pip_duration=self.pip_duration,
                                  pip_start=self.pip_start)
        
        preCell.set_sound_stim(self.stim, seed=seed)
        
        self['vm'] = postCell.soma(0.5)._ref_v
        #self['prevm'] = preCell.soma(0.5)._ref_v
        for i in range(30):
            self['xmtr%d'%i] = synapse.terminal.relsite._ref_XMTR[i]
            synapse.terminal.relsite.Dep_Flag = False
        self['t'] = h._ref_t
        
        h.tstop = 1e3*self.run_duration # duration of a run
        h.celsius = temp
        h.dt = dt
        
        self.custom_init()
        h.run()

    def show(self):
        self.win = pg.GraphicsWindow()
        Fs = self.Fs
        p1 = self.win.addPlot(title='stim', row=0, col=0)
        p1.plot(self.stim.time * 1000, self.stim.sound)
        p1.setXLink(p1)

        p2 = self.win.addPlot(title='AN spikes', row=1, col=0)
        vt = pg.VTickGroup(self.pre_cell._spiketrain)
        p2.addItem(vt)
        p2.setXLink(p1)

        p3 = self.win.addPlot(title='Bushy Spikes', row=2, col=0)
        bspk = PU.findspikes(self['t'], self['vm'], -30.)
        bspktick = pg.VTickGroup(bspk)
        p3.addItem(bspktick)
        p3.setXLink(p1)

        p4 = self.win.addPlot(title='Bushy Vm', row=3, col=0)
        p4.plot(self['t'], self['vm'])
        p4.setXLink(p1)

        p5 = self.win.addPlot(title='xmtr', row=0, col=1)
        for i in range(30):
            p5.plot(self['t'], self['xmtr%d'%i], pen=(i, 15))
        p5.setXLink(p1)
        
        p6 = self.win.addPlot(title='AN phase', row=1, col=1)
        phasewin = [self.pip_start[0] + 0.25*self.pip_duration, self.pip_start[0] + self.pip_duration]
        prespk = self.pre_cell._spiketrain
        spkin = prespk[np.where(prespk > phasewin[0]*1e3)]
        spikesinwin = spkin[np.where(spkin <= phasewin[1]*1e3)]
        vs = PU.vector_strength(spikesinwin, self.f0)
        print 'AN Vector Strength: %7.3f, d=%.2f (us) Rayleigh: %7.3f  p = %.3e  n = %d' % (vs['r'], vs['d']*1e6, vs['R'], vs['p'], vs['n'])
        (hist, binedges) = np.histogram(vs['ph'])
        curve = p6.plot(binedges, hist, stepMode=True, fillBrush=(100, 100, 255, 150), fillLevel=0)
        p6.setXRange(0., 2*np.pi)

        p7 = self.win.addPlot(title='Bushy phase', row=2, col=1)
        spkin = prespk[np.where(bspk > phasewin[0]*1e3)]
        spikesinwin = spkin[np.where(spkin <= phasewin[1]*1e3)]
        vs = PU.vector_strength(spikesinwin, self.f0)
        print 'BU Vector Strength: %7.3f, d=%.2f (us) Rayleigh: %7.3f  p = %.3e  n = %d' % (vs['r'], vs['d']*1e6, vs['R'], vs['p'], vs['n'])
        (hist, binedges) = np.histogram(vs['ph'])
        curve = p7.plot(binedges, hist, stepMode=True, fillBrush=(100, 100, 255, 150), fillLevel=0)
        p7.setXRange(0., 2*np.pi)
        p7.setXLink(p6)

        self.win.show()



prot = SGCInputTestPL()
prot.run()
prot.show()

import sys
if sys.flags.interactive == 0:
    pg.QtGui.QApplication.exec_()
