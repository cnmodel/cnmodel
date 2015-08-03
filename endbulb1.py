import sys
import numpy as np
import pyqtgraph as pg
from neuron import h
from cnmodel.protocols import Protocol
from cnmodel import cells
from cnmodel.util import get_anspikes
from cnmodel.util import sound
import pylibrary.pyqtgraphPlotHelpers as pgh
import pyqtgraph.multiprocess as mproc


class Endbulb(Protocol):
    def __init__(self):
        # create a cell
#        print dir(Endbulb)
        self.post_cell = cells.Bushy.create(species='mouse')
        self.post_cell.cell_initialize()  # initialize the cell to it's rmp
        #self.custom_init()
        self.h = h
        self._CF = None  # cell CF in Hz
        self._stimFreq = None  # stim frequency, in Hz
        self._dBSPL = None # dB
        self._stim = None
        self.spikeTrain = None # when using precomputed spike trains
        self.nconverge = 4
        self.srs = [1]*1+[2]*2+[3]*1
#        print dir(self)

    @property
    def CF(self):
        if self._CF is None:
            self._CF = 500.
        return self._CF

    @CF.setter
    def CF(self, cf):
        if not (100 <= cf <= 5e4):
            raise ValueError("CF must be in range 100 to 50,000 Hz!")
        else:
            self._CF = cf

    @property
    def stimFreq(self):
        if self._stimFreq is None:
            self._stimFreq = self.CF
        return self._stimFreq

    @stimFreq.setter
    def stimFreq(self, sfreq):
        if not (100 <= sfreq <= 5e4):
            raise ValueError("Stimulus frequency must be in range 100 to 50,000 Hz!")
        else:
            self._stimFreq = sfreq

    @property
    def dBSPL(self):
        if self._dBSPL is None:
            self._dBSPL = 80
        return self._dBSPL

    @dBSPL.setter
    def dBSPL(self, db):
        if not(0 <= db <= 120.):
            raise ValueError('dB SPL must be in range 0 - 120. dB')
        else:
            self._dBSPL = db

    @property
    def stim(self):
        if self._stim is None:
            self._stim = sound.TonePip(rate=100e3, duration=0.04, f0=self.stimFreq, dbspl=self.dBSPL,
                                  ramp_duration=2.5e-3, pip_duration=0.025,
                                  pip_start=[0.005])
        return self._stim

    def info(self):
        """
        get the run parameters, as a dictionary
        :return:
        """
        runinfo = {'CF': self._CF, 'StimFreq': self._stimFreq, 'SPL': self._dBSPL,
                   'nConverge': self.nconverge, 'SR': self.srs}
        return runinfo

    def run(self, temp=34.0, dt=0.025, seed=575982035, seedoffset=0):
        cfs = np.random.normal(loc=self.CF, scale=20, size=self.nconverge)  # jitter frequencies. 20Hz is adjacent HC at 4k
        self.pre_cell = [None]*self.nconverge
        self.synapse = [None]*self.nconverge
        for i in range(self.nconverge):
            preCell = cells.DummySGC(cf=cfs[i], sr=self.srs[i])
            synapse = preCell.connect(self.post_cell)
            synapse.terminal.relsite.Dep_Flag = False
            #synapse.terminal.relsite.nZones = 5
            synapse.terminal.relsite.dF = 0.2
            self.pre_cell[i] = preCell
            self.synapse[i] = synapse
        # print len(self.synapse[0].psd.ampa_psd)
        # print dir(self.synapse[0].psd.ampa_psd[0])

        sitemax = len(self.synapse[0].psd.ampa_psd)  # for instance...
        if sitemax > 50:
            sitemax = 50
        self.sitemax = sitemax
        self.vecs={}
        self.vecs['vm'] = self.h.Vector()
        self.vecs['vm'].record(self.post_cell.soma(0.5)._ref_v)
        self.vecs['t'] = self.h.Vector()
        self.vecs['t'].record(self.h._ref_t)
        k=2  # monitor the third input
        self.monitor = k
        for i, var in enumerate(['xmtr_%d_%d'%(i,k) for i in range(sitemax)]):
            self.vecs[var] = self.h.Vector()
            self.vecs[var].record(self.synapse[k].terminal.relsite._ref_XMTR[i])
        for i, var in enumerate(['gReceptor_%d_%d'%(i,k) for i in range(sitemax)]):
            self.vecs[var] = self.h.Vector()
            self.vecs[var].record(self.synapse[k].psd.ampa_psd[i]._ref_g)

        for i in range(self.nconverge):
            self.pre_cell[i].set_sound_stim(self.stim, seed=seed+i+seedoffset)

        self.h.tstop = 40.0 # duration of a run, msec
        self.h.celsius = temp
        self.h.dt = dt
        self.post_cell.vm0 = None  # force initialization
        self.post_cell.cell_initialize()  # initialize the cell to it's rmp
        self.custom_init()
#        Rin, tau, v = self.post_cell.measure_rintau(auto_initialize=False)
#        print '    *** Rin: %9.0f  tau: %9.1f   v: %6.1f' % (Rin, tau, v)
#        self.post_cell.cell_initialize()  # re-initialize the cell to it's rmp
#        self.custom_init()
        #self.h.run()
        while self.h.t < self.h.tstop:
            self.h.fadvance()
        return self.vecs_to_dict()


    def run_precomputed(self, temp=34.0, dt=0.025, seed=575982035, seedoffset=0):
        cfs = np.random.normal(loc=self.CF, scale=20, size=self.nconverge)  # jitter frequencies. 20Hz is adjacent HC at 4k
        self.pre_cell = [None]*self.nconverge
        self.synapse = [None]*self.nconverge
        for i in range(self.nconverge):
            preCell = cells.DummySGC() # cf=cfs[i], sr=self.srs[i])
            synapse = preCell.connect(self.post_cell)
            synapse.terminal.relsite.Dep_Flag = False
            #synapse.terminal.relsite.nZones = 5
            synapse.terminal.relsite.dF = 0.2
            self.pre_cell[i] = preCell
            self.synapse[i] = synapse
        # print len(self.synapse[0].psd.ampa_psd)
        # print dir(self.synapse[0].psd.ampa_psd[0])
        manager = get_anspikes.ManageANSpikes()
        self.spikeTrain = manager.getANatFandSPL(spontclass = 'MS', freq=10000., CF=self.CF, SPL=self.dBSPL)
        sitemax = len(self.synapse[0].psd.ampa_psd)  # for instance...
        if sitemax > 50:
            sitemax = 50
        self.sitemax = sitemax
        self.vecs={}
        self.vecs['vm'] = self.h.Vector()
        self.vecs['vm'].record(self.post_cell.soma(0.5)._ref_v)
        self.vecs['t'] = self.h.Vector()
        self.vecs['t'].record(self.h._ref_t)
        k=2  # monitor the third input
        self.monitor = k
        for i, var in enumerate(['xmtr_%d_%d'%(i,k) for i in range(sitemax)]):
            self.vecs[var] = self.h.Vector()
            self.vecs[var].record(self.synapse[k].terminal.relsite._ref_XMTR[i])
        for i, var in enumerate(['gReceptor_%d_%d'%(i,k) for i in range(sitemax)]):
            self.vecs[var] = self.h.Vector()
            self.vecs[var].record(self.synapse[k].psd.ampa_psd[i]._ref_g)

        self.maxt = 0
        for i in range(self.nconverge):
            self.pre_cell[i].set_spiketrain(1000*self.spikeTrain[i][:])
           # print i, 1000*self.spikeTrain[i][:]
            if len(self.spikeTrain[i][:]) > 0:
                self.maxt = np.max([self.maxt, 1000.*np.max(self.spikeTrain[i][:])])

        self.h.tstop = self.maxt # duration of a run, msec
        self.h.celsius = temp
        self.h.dt = dt
        self.post_cell.vm0 = None  # force initialization
        self.post_cell.cell_initialize()  # initialize the cell to it's rmp
        self.custom_init()
#        Rin, tau, v = self.post_cell.measure_rintau(auto_initialize=False)
#        print '    *** Rin: %9.0f  tau: %9.1f   v: %6.1f' % (Rin, tau, v)
#        self.post_cell.cell_initialize()  # re-initialize the cell to it's rmp
#        self.custom_init()
        #self.h.run()
        while self.h.t < self.h.tstop:
            self.h.fadvance()
        return self.vecs_to_dict()

    def getResults(self):
        return self.vecs_to_dict()  # make sure current
        
    def vecs_to_dict(self):
        self.d = {}
        for k in self.vecs.keys():
            #print 'k: ', k
            self.d[k] = np.array(self.vecs[k])  # turn in to numpy arrays
        # also add spike trains of pre cells
        for k in range(self.nconverge):
            self.d['preCell%03d'%k] = self.pre_cell[k]._spiketrain
        return self.d

    def show0(self):
        self.win = pg.GraphicsWindow()
        
        p1 = self.win.addPlot(title='Bushy Vm')
        p1.plot(np.array(self.vecs['t']), np.array(self.vecs['vm']))
        p2 = self.win.addPlot(title='XMTR', row=1, col=0)
        print 'receptor max: ', np.min(np.array(self.vecs['gReceptor_0_%d'%self.monitor]))
        for i in range(self.sitemax):
            p2.plot(self.vecs['t'], self.vecs['xmtr_%d_%d'%(i,self.monitor)], pen=(i, 15))
            #p2.plot(np.array(self.vecs['t']), np.array(self.vecs['gReceptor_%d_%d'%(i,self.monitor)]), pen=(i, 15))
        p2.setXLink(p1)
        p3 = self.win.addPlot(title='gPSD', row=2, col=0)
        #print 'receptor max: ', np.min(np.array(self.vecs['gReceptor_0_%d'%self.monitor]))
        for i in range(self.sitemax):
            #p2.plot(self.vecs['t'], self.vecs['xmtr_%d_%d'%(i,self.monitor)], pen=(i, 15))
            p3.plot(np.array(self.vecs['t']), np.array(self.vecs['gReceptor_%d_%d'%(i,self.monitor)]), pen=(i, 15))
        p3.setXLink(p1)
        p4 = self.win.addPlot(title='AN spikes', row=3, col=0)
        dist = 1./self.nconverge
        ytick = np.linspace(0, dist, self.nconverge+1)
        for i in range(self.nconverge):
            try:
                nsp = len(self.pre_cell[i]._spiketrain)
            except:
                nsp = len([self.pre_cell[i]._spiketrain])
            vt = pg.ScatterPlotItem(self.pre_cell[i]._spiketrain, [ytick[i+1]]*nsp, symbol='+', pen=(i,self.nconverge))
            #vt = pg.VTickGroup(self.pre_cell[i]._spiketrain, yrange=[ytick[i],ytick[i+1]], pen=(i,self.nconverge))
            p4.addItem(vt)
        p4.setXLink(p1)
        
        p5 = self.win.addPlot(title='stim', row=4, col=0)
        p5.plot(self.stim.time * 1000, self.stim.sound)
        p5.setXLink(p1)
        self.win.show()

    def show(self):
        self.win = pgh.figure(title='Endbulb1')
        layout = pgh.LayoutMaker(cols=1, rows=5, win=self.win, labelEdges=True, ticks='talbot')

        layout.title(0, title='Bushy Vm')
        layout.plot((0,0), np.array(self.vecs['t']), np.array(self.vecs['vm']), pen=pg.mkPen('k'))
        layout.title(1, title='XMTR')
        # print 'receptor max: ', np.min(np.array(self.vecs['gReceptor_0_%d'%self.monitor]))
        for i in range(self.sitemax):
            layout.plot(1, self.vecs['t'], self.vecs['xmtr_%d_%d'%(i,self.monitor)], pen=(i, 15))
            #p2.plot(np.array(self.vecs['t']), np.array(self.vecs['gReceptor_%d_%d'%(i,self.monitor)]), pen=(i, 15))
        layout.getPlot(1).setXLink(layout.getPlot(0))
        layout.title(2, title='gPSD')
        #print 'receptor max: ', np.min(np.array(self.vecs['gReceptor_0_%d'%self.monitor]))
        for i in range(self.sitemax):
            #p2.plot(self.vecs['t'], self.vecs['xmtr_%d_%d'%(i,self.monitor)], pen=(i, 15))
            layout.plot(2, np.array(self.vecs['t']), np.array(self.vecs['gReceptor_%d_%d'%(i,self.monitor)]), pen=(i, 15))
        layout.getPlot(2).setXLink(layout.getPlot(0))
        layout.title(3, title='AN spikes')
        dist = 1./self.nconverge
        ytick = np.linspace(0, dist, self.nconverge+1)
        for i in range(self.nconverge):
            try:
                nsp = len(self.pre_cell[i]._spiketrain)
            except:
                nsp = len([self.pre_cell[i]._spiketrain])
            vt = pg.ScatterPlotItem(self.pre_cell[i]._spiketrain, [ytick[i+1]]*nsp, symbol='+', pen=(i,self.nconverge))
            #vt = pg.VTickGroup(self.pre_cell[i]._spiketrain, yrange=[ytick[i],ytick[i+1]], pen=(i,self.nconverge))
            layout.getPlot(3).addItem(vt)
        layout.getPlot(3).setXLink(layout.getPlot(0))
        layout.title(4, title='stim')
        #layout.plot(4, self.stim.time, self.stim.sound, pg.mkPen('k'))
        layout.getPlot(4).setXLink(layout.getPlot(0))
        pgh.show()

if __name__ == '__main__':
    prot = Endbulb()
    prot.CF=4000
    prot.stimFreq = 4000
    print ('dbSPL: {:.1f}  CF: {:.1f}  stim Freq: {:.1f}'.format(prot.dBSPL, prot.CF, prot.stimFreq))
    prot.run_precomputed()
    prot.show()

    import sys
    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()
