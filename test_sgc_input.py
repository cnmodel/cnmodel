import sys
import numpy as np
import pyqtgraph as pg
from nrnlibrary.protocols import Protocol
from nrnlibrary import cells
from neuron import h

class SGCInputTest(Protocol):
    def run(self, temp=34.0, dt=0.025):
        preCell = cells.DummySGC()
        postCell = cells.Bushy.create()
        synapse = preCell.connect(postCell)
        
        preCell.set_spiketrain(np.linspace(10, 50, 10))
        
        self['vm'] = postCell.soma(0.5)._ref_v
        #self['prevm'] = preCell.soma(0.5)._ref_v
        self['xmtr'] = synapse.terminal.relsite._ref_XMTR[0]
        self['t'] = h._ref_t
        
        h.tstop = 100.0 # duration of a run
        h.celsius = temp
        h.dt = dt
        
        self.custom_init()
        h.run()

    def show(self):
        self.win = pg.GraphicsWindow()
        p1 = self.win.addPlot(title='Vm')
        p1.plot(self['t'], self['vm'])
        p2 = self.win.addPlot(title='xmtr', row=1, col=0)
        p2.plot(self['t'], self['xmtr'])
        # this should be completely flat:
        #p3 = self.win.addPlot(title='Vm pre', row=2, col=0)
        #p3.plot(self['t'], self['prevm'])
        self.win.show()


prot = SGCInputTest()
prot.run()
prot.show()

import sys
if sys.flags.interactive == 0:
    pg.QtGui.QApplication.exec_()
