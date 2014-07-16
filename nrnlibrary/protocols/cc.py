from neuron import h
import numpy as np
import scipy
import scipy.integrate
import scipy.stats

try:
    import pyqtgraph as pg
    HAVE_PG = True
except ImportError:
    HAVE_PG = False

from ..util.stim import make_pulse
from .protocol import Protocol

class CurrentClamp(Protocol):
    def __init__(self):
        super(CurrentClamp, self).__init__()
    
    def run(self, cell, cmd, temp=22, dt=0.025):
        """
        Run a single current-clamp recording on *section*.
        
        Parameters:
        cell : Cell
            The Cell instance to test. IClamp will be attached to 
            cell.soma(0.5).
        cmd : array
            Array of current command values
        temp : 
            temperature of simulation (22)
        dt : 
            timestep of simulation (0.025)
        """
        self.reset()
        self.cell = cell
        self.current_cmd = cmd
        self.dt = dt
        self.temp = temp
        
        tend = len(cmd) * dt
        
        # Configure IClamp
        istim = h.iStim(0.5, sec=cell.soma)
        istim.delay = 0
        istim.dur = 1e9 # these actually do not matter...
        istim.iMax = 0.0
        i_stim = h.Vector(cmd)
        i_stim.play(istim._ref_i, h.dt, 0)
        
        # Connect recording vectors
        self['vm'] = cell.soma(0.5)._ref_v
        self['i_inj'] = istim._ref_i
        self['time'] = h._ref_t


        # GO
        h.dt = dt
        h.celsius = temp
        h.tstop = tend
        cell.initialize()
        h.frecord_init()
        while h.t < h.tstop:
            h.fadvance()
            

    def show(self):
        """
        Plot results from run()
        """
        if not HAVE_PG:
            raise Exception("Requires pyqtgraph")
        
        #
        # Generate figure with subplots
        #
        app = pg.mkQApp()
        win = pg.GraphicsWindow()
        win.resize(1000, 800)
        Vplot = win.addPlot(labels={'left': 'Vm (mV)', 'bottom': 'Time (ms)'})
        win.nextRow()
        Iplot = win.addPlot(labels={'left': 'Iinj (nA)', 'bottom': 'Time (ms)'})
        
        win.ci.layout.setRowStretchFactor(0, 10)
        win.ci.layout.setRowStretchFactor(1, 5)

        Vplot.plot(self['time'], self['vm'])
        Iplot.plot(self['time'], self['i_inj'])
        
        self.win = win

    
