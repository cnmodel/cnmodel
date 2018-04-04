"""
Test generating a series of EPSPs in a MSO cell in response to tone pip. This
test demonstrates how a binaural circuit can be constructed.

This script:

1. Creates multiple SGCs (base instance has 3) converging onto two bushy cells.
     The 2 bushy cells then converge onto one MSO cell
2. Connects the group of sgc cells from one ear to one bushy cell, and the 
    other sgcs from the other ear to the
    other bushy cell.
3. Specifies the CFs of the SGC cells individually. Also specifies the frequency of 
    and stimuli by ear allowing for "binaural beats"
5. Records the bushy and MSO cell membrane voltages, sgc spike time, and calculates
    vector strengths.

The auditory nerve spike train is generated automatically by the DummySGC class
using the tone pip. For lower-level access to the auditory nerve model, see the
test_an_model.py and test_sound_stim.py examples.

Usage:
    python examples/test_mso_inputs.py [cochlea | matlab]

    The AN simulator that is run depends on what is available and how the script is called. 
    python examples/test_mso_inputs.py [cochlea | matlab] will try to use the specified
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
import cnmodel.util.pynrnutilities as PU

class MSOBinauralTest(Protocol):
    def run(self, temp=38.0, dt=0.025, seed=575982035, simulator=None):
        ears = {'left': [500., 502, 498], 'right': [500., 502., 498]}
        
        self.beatfreq = 0.
        self.f0 = 500.
        f0 = {'left': self.f0, 'right': self.f0+self.beatfreq}
        nsgc = len(ears.keys())
        sgcCell = {}
        bushyCell = {}
        msoCell = {}
        synapse = {}
        self.stim = {}
        self.ears = ears
        self.stimdur = 0.2
        self.stimdelay = 0.02
        self.rundur = self.stimdelay + self.stimdur + 0.02

        for i, ear in enumerate(ears.keys()):
            nsgc = len(ears[ear])  # how many sgcs are specified for this ear
            sgcCell[ear] = [cells.DummySGC(cf=ears[ear][k], sr=2) for k in range(nsgc)]
            bushyCell[ear] = [cells.Bushy.create(temperature=temp)]
            synapse[ear] = [sgcCell[ear][k].connect(bushyCell[ear][0]) for k in range(nsgc)]
            self.stim[ear] = [sound.TonePip(rate=100e3, duration=self.stimdur+0.1, f0=f0[ear], dbspl=80,
                                  ramp_duration=2.5e-3, pip_duration=self.stimdur, 
                                  pip_start=[self.stimdelay]) for k in range(nsgc)]
            for k in range(len(self.stim[ear])):
                sgcCell[ear][k].set_sound_stim(self.stim[ear][k], seed=seed + i*seed + k, simulator=simulator)
            self['vm_bu_%s' % ear] = bushyCell[ear][0].soma(0.5)._ref_v
            for k in range(30):
                self['xmtr%d_%s'%(k, ear)] = synapse[ear][0].terminal.relsite._ref_XMTR[k]
            for k in range(len(synapse[ear])):
                synapse[ear][k].terminal.relsite.Dep_Flag = False  # turn off depression

        msoCell = cells.MSO.create(temperature=temp)  # one target MSO cell
        msosyn = {}
        for ear in ears:
            msosyn[ear] = bushyCell[ear][0].connect(msoCell)
        self.sgc_cells = sgcCell
        self.bushy_cells = bushyCell
        self.synapses = synapse
        self.msyns = msosyn
        self.msoCell = msoCell
        self.all_cells = []  # hold all "real" cells (DummySGC does not have mechanisms)
        for ear in ears.keys():
            self.all_cells.append([c for c in self.bushy_cells[ear]])
        self.all_cells.append([self.msoCell])
        
        self['vm_mso'] = self.msoCell.soma(0.5)._ref_v
        for k, ear in enumerate(ears.keys()):
            for i in range(30):
                self['mso_xmtr%d_%s'%(i, ear)] = msosyn[ear].terminal.relsite._ref_XMTR[i]
            msosyn[ear].terminal.relsite.Dep_Flag = False  # turn off depression

        self['t'] = h._ref_t
        
        h.tstop = self.rundur*1e3 # duration of a run
        h.celsius = temp
        h.dt = dt
        
        custom_init()
        # confirm that all cells are ok
        for cg in self.all_cells:
            for c in cg:
                c.check_all_mechs()
        while h.t < h.tstop:
            h.fadvance()
            

    def show(self):
        self.win = pg.GraphicsWindow()
        
        p5 = self.win.addPlot(title='stim')
        p5.plot(self.stim['left'][0].time * 1000., self.stim['left'][0].sound)

        p1 = self.win.addPlot(title='Bushy Vm', row=1, col=0)
        for k, ear in enumerate(self.ears.keys()):
            p1.plot(self['t'], self['vm_bu_%s' % ear], pen=(k, 15))
        
        p2 = self.win.addPlot(title='SGC-BU xmtr left', row=0, col=1)
        for i in range(30):
            p2.plot(self['t'], self['xmtr%d_left'%i], pen=(i, 15))
        p2.setXLink(p1)
        p2r = self.win.addPlot(title='SGC-BU xmtr right', row=1, col=1)
        for i in range(30):
            p2r.plot(self['t'], self['xmtr%d_right'%i], pen=(i, 15))
        p2r.setXLink(p1)        
        
        p3 = self.win.addPlot(title='MSO Vm', row=2, col=0)
        p3.plot(self['t'], self['vm_mso'])
        p3.setXLink(p1)
        
        p4 = self.win.addPlot(title='BU-MSO xmtr', row=2, col=1)
        for k, ear in enumerate(self.ears.keys()):
            for i in range(30):
                p2.plot(self['t'], self['mso_xmtr%d_%s'%(i, ear)], pen=(i, 15))
        p4.setXLink(p1)        

        p4 = self.win.addPlot(title='AN spikes', row=3, col=0)
        ntrain = len(self.sgc_cells['left'])
        for k in range(ntrain):
            yr = [k/float(ntrain), (k+0.8)/float(ntrain)]
            vt = pg.VTickGroup(self.sgc_cells['left'][k]._spiketrain, yrange = yr, pen=(k, 15))
            p4.addItem(vt)
        p4.setXLink(p1)
        p5.setXLink(p1)
        
        # phaselocking calculations
        phasewin = [self.stimdelay + 0.2*self.stimdur, self.stimdelay + self.stimdur]
        msospk = PU.findspikes(self['t'], self['vm_mso'], -30.)

        spkin = msospk[np.where(msospk > phasewin[0]*1e3)]
        spikesinwin = spkin[np.where(spkin <= phasewin[1]*1e3)[0]]

        # set freq for VS calculation
        f0 = self.f0
        fb = self.beatfreq
        vs = PU.vector_strength(spikesinwin, f0)

        print 'MSO Vector Strength at %.1f: %7.3f, d=%.2f (us) Rayleigh: %7.3f  p = %.3e  n = %d' % (f0, vs['r'], vs['d']*1e6, vs['R'], vs['p'], vs['n'])
        if fb > 0:
            vsb = PU.vector_strength(spikesinwin, fb)
            print 'MSO Vector Strength to beat at %.1f: %7.3f, d=%.2f (us) Rayleigh: %7.3f  p = %.3e  n = %d' % (fb, vsb['r'], vsb['d']*1e6, vsb['R'], vsb['p'], vsb['n'])
        (hist, binedges) = np.histogram(vs['ph'])
        p6 = self.win.addPlot(title='VS', row=3, col=1)
        p6.plot(binedges, hist, stepMode=True, fillBrush=(100, 100, 255, 150), fillLevel=0)
        p6.setXRange(0., 2*np.pi)
        
        self.win.show()


if __name__ == '__main__':
    simulator = None
    if len(sys.argv) > 1:
        simulator=sys.argv[1]
    prot = MSOBinauralTest()
    prot.run(simulator=simulator)
    prot.show()

    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()
