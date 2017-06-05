import sys
import numpy as np
import pyqtgraph as pg
from neuron import h
from cnmodel.protocols import Protocol
from cnmodel import cells
from cnmodel.util import sound
import cnmodel.util.pynrnutilities as PU
from cnmodel import data


class SGCInputTestPL(Protocol):
    def set_cell(self, cell='bushy'):
        self.cell = cell
    
    def run(self, temp=34.0, dt=0.025, seed=575982035):
        if self.cell == 'bushy':
            postCell = cells.Bushy.create()
        elif self.cell == 'tstellate':
            postCell = cells.TStellate.create()
        elif self.cell == 'octopus':
            postCell = cells.Octopus.create()
        elif self.cell == 'dstellate':
            postCell = cells.DStellate.create()
        else:
            raise ValueError('cell %s is not yet implemented for phaselocking')
        self.post_cell = postCell
        self.run_duration = 1.  # in seconds
        self.pip_duration = 0.8  # in seconds
        self.pip_start = [0.02]  # in seconds
        self.Fs = 100e3  # in Hz
        self.f0 = 1500.  # stimulus in Hz
        self.cf = 1500.  # SGCs in Hz
        self.stim = sound.TonePip(rate=self.Fs, duration=self.run_duration, f0=self.f0, dbspl=60,
                                  ramp_duration=2.5e-3, pip_duration=self.pip_duration,
                                  pip_start=self.pip_start)
        
        n_sgc = data.get('convergence', species='mouse', post_type=postCell.type, pre_type='sgc')[0]
        self.n_sgc = int(np.round(n_sgc))

        self.pre_cells = []
        self.synapses = []
        j = 0
        for k in range(self.n_sgc):
            seed = seed + k
            preCell = cells.DummySGC(cf=self.cf, sr=2)
            synapse = preCell.connect(postCell)
            for i in range(synapse.terminal.n_rzones):
                self['xmtr%03d'%j] = synapse.terminal.relsite._ref_XMTR[i]
                j = j + 1
            synapse.terminal.relsite.Dep_Flag = False
            preCell.set_sound_stim(self.stim, seed=seed)
            self.pre_cells.append(preCell)
            self.synapses.append(synapse)
        
        self['vm'] = postCell.soma(0.5)._ref_v
        #self['prevm'] = preCell.soma(0.5)._ref_v
        self['t'] = h._ref_t
        postCell.cell_initialize()
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
        vt = pg.VTickGroup(self.pre_cells[0]._spiketrain)
        p2.addItem(vt)
        p2.setXLink(p1)

        p3 = self.win.addPlot(title='%s Spikes' % self.cell, row=2, col=0)
        bspk = PU.findspikes(self['t'], self['vm'], -30.)
        bspktick = pg.VTickGroup(bspk)
        p3.addItem(bspktick)
        p3.setXLink(p1)

        p4 = self.win.addPlot(title='%s Vm' % self.cell, row=3, col=0)
        p4.plot(self['t'], self['vm'])
        p4.setXLink(p1)

        p5 = self.win.addPlot(title='xmtr', row=0, col=1)
        j = 0
        for k in range(self.n_sgc):
            synapse = self.synapses[k]
            for i in range(synapse.terminal.n_rzones):
                p5.plot(self['t'], self['xmtr%03d'%j], pen=(i, 15))
                j = j + 1
        p5.setXLink(p1)
        
        p6 = self.win.addPlot(title='AN phase', row=1, col=1)
        phasewin = [self.pip_start[0] + 0.25*self.pip_duration, self.pip_start[0] + self.pip_duration]
        prespk = self.pre_cells[0]._spiketrain  # just sample one...
        spkin = prespk[np.where(prespk > phasewin[0]*1e3)]
        spikesinwin = spkin[np.where(spkin <= phasewin[1]*1e3)]
        vs = PU.vector_strength(spikesinwin, self.f0)
        print 'AN Vector Strength: %7.3f, d=%.2f (us) Rayleigh: %7.3f  p = %.3e  n = %d' % (vs['r'], vs['d']*1e6, vs['R'], vs['p'], vs['n'])
        (hist, binedges) = np.histogram(vs['ph'])
        curve = p6.plot(binedges, hist, stepMode=True, fillBrush=(100, 100, 255, 150), fillLevel=0)
        p6.setXRange(0., 2*np.pi)

        p7 = self.win.addPlot(title='%s phase' % self.cell, row=2, col=1)
        spkin = bspk[np.where(bspk > phasewin[0]*1e3)]
        spikesinwin = spkin[np.where(spkin <= phasewin[1]*1e3)]
        vs = PU.vector_strength(spikesinwin, self.f0)
        print '%s Vector Strength: %7.3f, d=%.2f (us) Rayleigh: %7.3f  p = %.3e  n = %d' % (self.cell, vs['r'], vs['d']*1e6, vs['R'], vs['p'], vs['n'])
        (hist, binedges) = np.histogram(vs['ph'])
        curve = p7.plot(binedges, hist, stepMode=True, fillBrush=(100, 100, 255, 150), fillLevel=0)
        p7.setXRange(0., 2*np.pi)
        p7.setXLink(p6)

        self.win.show()


if __name__ == '__main__':
    if len(sys.argv) > 1:
        cell = sys.argv[1]
    else:
        cell = 'bushy'
    print 'cell type: ', cell
    prot = SGCInputTestPL()
    prot.set_cell(cell)
    prot.run()
    prot.show()

    import sys
    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()
