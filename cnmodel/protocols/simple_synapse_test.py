from neuron import h
import pyqtgraph as pg

import cnmodel.util as util
from .protocol import Protocol
from .. import cells


class SimpleSynapseTest(Protocol):
    def reset(self):
        super(SimpleSynapseTest, self).reset()

    def run(self, pre_sec, post_sec, temp=34.0, dt=0.025, 
            vclamp=-65.0, iterations=1, tstop=200.0, stim_params=None, **kwds):
        """ 
        """
        Protocol.run(self, **kwds)
        
        pre_cell = cells.cell_from_section(pre_sec)
        post_cell = cells.cell_from_section(post_sec)
        synapse = pre_cell.connect(post_cell, type='multisite')
#        synapse = pre_cell.connect(post_cell, type='simple')
        self.synapse = synapse
        self.pre_sec = synapse.terminal.section
        self.post_sec = synapse.psd.section
        self.pre_cell = pre_cell
        self.post_cell = post_cell
        
        #
        # voltage clamp the target cell
        #
        vccontrol = h.VClamp(0.5, sec=post_cell.soma)
        vccontrol.dur[0] = tstop
        vccontrol.amp[0] = vclamp
        #vccontrol.dur[1] = 100.0
        #vccontrol.amp[1] = clampV
        #vccontrol.dur[2] = 20.0
        #vccontrol.amp[2] = clampV

        #
        # set up stimulation of the presynaptic axon/terminal
        #
        
        istim = h.iStim(0.5, sec=pre_cell.soma)
        stim = {
            'NP': 10,
            'Sfreq': 100.0,
            'delay': 10.0,
            'dur': 0.5,
            'amp': 10.0,
            'PT': 0.0,
            'dt': dt,
        }
        if stim_params is not None:
            stim.update(stim_params)
        (secmd, maxt, tstims) = util.make_pulse(stim)
        self.stim = stim

        if tstop is None:
            tstop = len(secmd) * dt
        
        istim.delay = 0
        istim.dur = 1e9 # these actually do not matter...
        istim.iMax = 0.0

        # istim current pulse train
        i_stim_vec = h.Vector(secmd)
        i_stim_vec.play(istim._ref_i, dt, 0)

        #
        # Run simulation
        #
        h.tstop = tstop # duration of a run
        h.celsius = temp
        h.dt = dt
        self.temp = temp
        self.dt = dt
        for nrep in xrange(iterations): # could do multiple runs.... 
            self.reset()
            self['v_pre'] = pre_cell.soma(0.5)._ref_v
            self['t'] = h._ref_t
            self['v_soma'] = post_cell.soma(0.5)._ref_v
            self['i_soma'] = vccontrol._ref_i
            util.custom_init()
            h.run()

    def show(self):
        self.win = pg.GraphicsWindow()
        self.win.resize(800, 800)
        t = self['t']

        p1 = self.win.addPlot(title=self.pre_cell.status['name'])
        p1.setLabels(left='V pre (mV)', bottom='Time (ms)')
        p1.plot(t, self['v_pre'])
        
        self.win.nextRow()
        p2 = self.win.addPlot(title=self.post_cell.status['name'])
        p2.plot(t[1:], self['i_soma'][1:], pen=pg.mkPen('w', width=2))
        p2.setLabels(left='I post (nA)', bottom='Time (ms)')
        p2.setXLink(p1)
