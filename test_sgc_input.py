import sys
import numpy as np
import pyqtgraph as pg
from neuron import h
from cnmodel.protocols import Protocol
from cnmodel import cells
from cnmodel.util import sound


class SGCInputTest(Protocol):
    def run(self, temp=34.0, dt=0.025, seed=575982035):
        preCell = cells.DummySGC(cf=4000, sr=2)
        postCell = cells.Bushy.create()
        synapse = preCell.connect(postCell)
        self.pre_cell = preCell
        self.post_cell = postCell
        self.synapse = synapse
 
        self.stim = sound.TonePip(rate=100e3, duration=0.1, f0=4000, dbspl=80,
                                  ramp_duration=2.5e-3, pip_duration=0.04, 
                                  pip_start=[0.02])
        
        preCell.set_sound_stim(self.stim, seed=seed)
        
        self['vm'] = postCell.soma(0.5)._ref_v
        #self['prevm'] = preCell.soma(0.5)._ref_v
        for i in range(30):
            self['xmtr%d'%i] = synapse.terminal.relsite._ref_XMTR[i]
            synapse.terminal.relsite.Dep_Flag = False
        self['t'] = h._ref_t
        
        h.tstop = 100.0 # duration of a run
        h.celsius = temp
        h.dt = dt
        
        self.custom_init()
        h.run()

    def show(self):
        self.win = pg.GraphicsWindow()
        
        p1 = self.win.addPlot(title='Bushy Vm')
        p1.plot(self['t'], self['vm'])
        p2 = self.win.addPlot(title='xmtr', row=1, col=0)
        for i in range(30):
            p2.plot(self['t'], self['xmtr%d'%i], pen=(i, 15))
        p2.setXLink(p1)
        
        p3 = self.win.addPlot(title='AN spikes', row=2, col=0)
        vt = pg.VTickGroup(self.pre_cell._spiketrain)
        p3.addItem(vt)
        p3.setXLink(p1)
        
        p4 = self.win.addPlot(title='stim', row=3, col=0)
        p4.plot(self.stim.time * 1000, self.stim.sound)
        p4.setXLink(p1)
        self.win.show()


prot = SGCInputTest()
prot.run()
prot.show()

import sys
if sys.flags.interactive == 0:
    pg.QtGui.QApplication.exec_()
