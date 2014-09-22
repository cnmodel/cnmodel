import numpy as np
import pyqtgraph as pg

from neuron import h

import nrnlibrary.util as util
from .protocol import Protocol


class PopulationTest(Protocol):
    def reset(self):
        super(PopulationTest, self).reset()

    def run(self, pops, cf=4000, temp=34.0, dt=0.025):
        """ 
        1. Connect pop1 => pop2
        2. Instantiate a single cell in pop2
        3. Automatically generate presynaptic cells and synapses from pop1
        4. Stimulate presynaptic cells and record postsynaptically
        """
        
        pre_pop, post_pop = pops
        pre_pop.connect(post_pop)
        
        # start with one cell, selected from the user-selected population, that has
        # a cf close to 4kHz
        post_cell_ind = post_pop.select(1, cf=cf, create=True)
        post_cell = post_pop.get_cell(post_cell_ind)
        post_pop.resolve_inputs(depth=1)
        
        
        post_sec = post_cell.soma
        self.post_cell = post_cell
        
        pre_cell_inds = post_pop.cell_connections(post_cell_ind)
        pre_cells = [pre_pop.get_cell(i) for i in pre_cell_inds]
        pre_secs = [cell.soma for cell in pre_cells]
        self.pre_cells = pre_cells
        
        
        #
        # voltage clamp the target cell
        #
        clampV = 40.0
        vccontrol = h.VClamp(0.5, sec=post_cell.soma)
        vccontrol.dur[0] = 10.0
        vccontrol.amp[0] = clampV
        vccontrol.dur[1] = 100.0
        vccontrol.amp[1] = clampV
        vccontrol.dur[2] = 20.0
        vccontrol.amp[2] = clampV

        #
        # set up stimulation of the presynaptic cells
        #
        self.stim_params = []
        self.istim = []
        for pre_cell in pre_cells:
            istim = h.iStim(0.5, sec=pre_cell.soma)
            stim = {}
            stim['NP'] = 10
            stim['Sfreq'] = 100.0 # stimulus frequency
            stim['delay'] = 10.0
            stim['dur'] = 0.5
            stim['amp'] = 10.0
            stim['PT'] = 0.0
            stim['dt'] = dt
            (secmd, maxt, tstims) = util.make_pulse(stim)
            self.stim_params.append(stim)
        
            istim.delay = 0
            istim.dur = 1e9 # these actually do not matter...
            istim.iMax = 0.0

            # istim current pulse train
            i_stim_vec = h.Vector(secmd)
            i_stim_vec.play(istim._ref_i, dt, 0)
            self.istim.append(istim)
            self['istim'] = istim._ref_i

        # create hoc vectors for each parameter we wish to monitor and display
        self['v_pre'] = pre_cells[0].soma(0.5)._ref_v
        self['t'] = h._ref_t
        self['v_soma'] = pre_cell.soma(0.5)._ref_v

        #
        # Run simulation
        #
        h.tstop = 200.0
        h.celsius = temp
        h.dt = dt
        self.temp = temp
        self.dt = dt
        
        self.custom_init()
        h.run()
            


    def show(self):
        self.win = pg.GraphicsWindow()
        self.win.resize(1000, 1000)
        
        cmd_plot = self.win.addPlot()
        cmd_plot.plot(self['t'], self['istim'])
        
        self.win.nextRow()
        pre_plot = self.win.addPlot()
        cmd_plot.plot(self['t'], self['v_pre'])

        self.win.nextRow()
        post_plot = self.win.addPlot()
        cmd_plot.plot(self['t'], self['v_post'])
        
        