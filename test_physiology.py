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
        
    def plot_vsoma(self):
        self.plot = pg.plot()
        real = self.bushy.real_cells()
        for i, ind in enumerate(real):
            cell = self.bushy.get_cell(ind)
            self.plot.plot(self['t'], self[cell], pen=(i, len(real)*1.5))
        
        # plot ticks for first cell's SGC input
        sgc_ind = self.bushy._cells[real[0]]['connections'][self.sgc]
        self.plot.ticks = []
        for i, ind in enumerate(sgc_ind):
            cell = self.sgc.get_cell(ind)
            spikes = cell._spiketrain
            vticks = pg.VTickGroup(spikes, pen=(i, len(sgc_ind)*1.5), yrange=(0.9, 1))
            self.plot.ticks.append(vticks)
            self.plot.addItem(vticks)
            
        ds_ind = self.bushy._cells[real[0]]['connections'][self.dstellate]
        for i, ind in enumerate(ds_ind):
            cell = self.dstellate.get_cell(ind)
            vm = self[cell]
            spike_inds = np.argwhere((vm[1:]>-20) & (vm[:-1]<=-20))[:,0]
            spikes = self['t'][spike_inds]
            vticks = pg.VTickGroup(spikes, pen=(i, len(ds_ind)*1.5), yrange=(0.8, 0.9))
            self.plot.ticks.append(vticks)
            self.plot.addItem(vticks)
            
        
        
if __name__ == '__main__':
    # Create a sound stimulus and use it to generate spike trains for the SGC
    # population
    stim = sound.TonePip(rate=100e3, duration=0.1, f0=16e3, dbspl=80,
                         ramp_duration=2.5e-3, pip_duration=0.04, 
                         pip_start=[0.02])

    prot = CNSoundStim(stim, seed=34657845)
    prot.run()

    prot.plot_vsoma()
