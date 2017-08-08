"""
Test generating a series of EPSPs in a bushy cell in response to tone pip.

This script:

1. Creates one SGC and one bushy cell.
2. Connects the two cells with a synapse.
3. Specifies a tone pip input to the SGC cell.
4. Records the bushy cell Vm.

The auditory nerve spike train is generated automatically by the DummySGC class
using the tone pip. For lower-level access to the auditory nerve model, see the
test_an_model.py and test_sound_stim.py examples.

The simulator that is run depends on what is available and how the script is called. 
python examples/test_sgc_input.py [cochlea | matlab] will try to use the specified
simulator. If no simulator is specified, it will try to use cochlea or matlab, in
that order.

"""
import sys
import numpy as np
import pyqtgraph as pg
from neuron import h
from cnmodel.protocols import Protocol
from cnmodel import cells
from cnmodel.util import sound
from cnmodel.util import custom_init


class SGCInputTest(Protocol):
    def run(self, temp=34.0, dt=0.025, seed=575982035, simulator=None):
        preCell = cells.DummySGC(cf=4000, sr=2)
        postCell = cells.Bushy.create()
        synapse = preCell.connect(postCell)
        self.pre_cell = preCell
        self.post_cell = postCell
        self.synapse = synapse
 
        self.stim = sound.TonePip(rate=100e3, duration=0.1, f0=4000, dbspl=80,
                                  ramp_duration=2.5e-3, pip_duration=0.04, 
                                  pip_start=[0.02])
        
        preCell.set_sound_stim(self.stim, seed=seed, simulator=simulator)
        
        self['vm'] = postCell.soma(0.5)._ref_v
        #self['prevm'] = preCell.soma(0.5)._ref_v
        for i in range(30):
            self['xmtr%d'%i] = synapse.terminal.relsite._ref_XMTR[i]
            synapse.terminal.relsite.Dep_Flag = False
        self['t'] = h._ref_t
        
        h.tstop = 100.0 # duration of a run
        h.celsius = temp
        h.dt = dt
        
        custom_init()
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


if __name__ == '__main__':
    simulator = None
    if len(sys.argv) > 1:
        simulator=sys.argv[1]
    prot = SGCInputTest()
    prot.run(simulator=simulator)
    prot.show()

    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()
