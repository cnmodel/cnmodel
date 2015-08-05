"""
Test physiological response properties.
"""

import numpy as np
import pyqtgraph as pg
from cnmodel import populations
from cnmodel.util import sound, random
from cnmodel.protocols import Protocol
from neuron import h


class CNSoundStim(Protocol):
    def __init__(self, stim, seed, temp=34.0, dt=0.025):
        Protocol.__init__(self)
        
        random.set_seed(seed)
        self.stim = stim
        self.temp = temp
        self.dt = dt

        # Create cell populations.
        # This creates a complete set of _virtual_ cells for each population. No 
        # cells are instantiated at this point.
        self.sgc = populations.SGC(model='dummy')
        self.bushy = populations.Bushy()
        self.dstellate = populations.DStellate()
        #self.tstellate = populations.TStellate()

        # Connect populations. 
        # This only defines the connections between populations; no synapses are 
        # created at this stage.
        self.sgc.connect(self.bushy, self.dstellate)#, self.tstellate)
        self.dstellate.connect(self.bushy)#, self.tstellate)
        #self.tstellate.connect(self.bushy)

        # Select cells to record from.
        # At this time, we actually instantiate the selected cells.
        # select 10 bushy cells closest to 16kHz
        bushy_cell_ids = self.bushy.select(10, cf=16e3, create=True)  
        # select 10 stellate cells closest to 16kHz
        #tstel_cell_ids = self.tstellate.select(10, cf=16e3, create=True)  

        # Now create the supporting circuitry needed to drive the cells we selected.
        # At this time, cells are created in all populations and automatically 
        # connected with synapses.
        self.bushy.resolve_inputs(depth=2)
        #self.tstellate.resolve_inputs(depth=2)
        # Note that using depth=2 indicates the level of recursion to use when 
        # resolving inputs. For example, resolving inputs for the bushy cell population
        # (level 1) creates presynaptic cells in the dstellate population, and resolving
        # inputs for the dstellate population (level 2) creates presynaptic cells in the
        # sgc population. 

        self.sgc.set_sound_stim(stim, seed=seed)

    def run(self):
        self.reset()
        
        # set up recording vectors
        for pop in self.bushy, self.dstellate:
            for ind in pop.real_cells():
                cell = pop.get_cell(ind)
                self[cell] = cell.soma(0.5)._ref_v
        self['t'] = h._ref_t
            
        h.tstop = self.stim.duration * 1000
        h.celsius = self.temp
        h.dt = self.dt
        
        print "init.."
        self.custom_init()
        print "start.."
        while h.t < h.tstop:
            h.fadvance()
            print "%0.2f / %0.2f" % (h.t, h.tstop)
            

class NetworkSimDisplay(pg.QtGui.QWidget):
    def __init__(self, prot):
        pg.QtGui.QWidget.__init__(self)
        
        self.selected = None
        
        self.prot = prot
        self.layout = pg.QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.pw = pg.GraphicsLayoutWidget()
        self.layout.addWidget(self.pw, 0, 0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        
        self.cell_plot = self.pw.addPlot()
        real = self.prot.bushy.real_cells()
        for i, ind in enumerate(real):
            cell = self.prot.bushy.get_cell(ind)
            p = self.cell_plot.plot(self.prot['t'], self.prot[cell], pen=(i, len(real)*1.5), name='bushy-%d' % ind)
            p.curve.setClickable(True)
            p.sigClicked.connect(self.cell_curve_clicked)
            p.cell_ind = ind
        self.cell_plot.addLegend()

        self.input_plot = self.pw.addPlot(row=1, col=0)
        self.input_plot.setXLink(self.cell_plot)
        
    def cell_curve_clicked(self, c):
        if self.selected is not None:
            pen = self.selected.curve.opts['pen']
            pen.setWidth(1)
            self.selected.setPen(pen)
            
        pen = c.curve.opts['pen']
        pen.setWidth(3)
        c.setPen(pen)
        self.selected = c
        print c, c.cell_ind

        self.show_cell(c.cell_ind)
        
    def show_cell(self, ind):
        """Show spike trains of inputs to selected cell.
        """
        self.input_plot.clear()
        rec = self.prot.bushy._cells[ind]
        i = 0
        plots = []
        for j, pop in enumerate((self.prot.dstellate, self.prot.sgc)):
            if pop not in rec['connections']:
                continue
            pre_inds = rec['connections'][pop]
            for ind in pre_inds:
                cell = pop.get_cell(ind)
                if hasattr(cell, '_spiketrain'):
                    spikes = cell._spiketrain
                else:
                    vm = self.prot[cell]
                    spike_inds = np.argwhere((vm[1:]>-20) & (vm[:-1]<=-20))[:,0]
                    spikes = self.prot['t'][spike_inds]
                y = np.ones(len(spikes)) * i
                self.input_plot.plot(spikes, y, pen=None, symbol=['o', 't'][j], symbolBrush=(i, 30))
                i += 1
                    
        
        
if __name__ == '__main__':
    app = pg.mkQApp()
    
    # Create a sound stimulus and use it to generate spike trains for the SGC
    # population
    stim = sound.TonePip(rate=100e3, duration=0.1, f0=16e3, dbspl=80,
                         ramp_duration=2.5e-3, pip_duration=0.04, 
                         pip_start=[0.02])

    prot = CNSoundStim(stim, seed=34657845)
    prot.run()

    #prot.plot_vsoma()
    nd = NetworkSimDisplay(prot)
    nd.show()
