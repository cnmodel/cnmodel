import sys
import numpy as np
import pyqtgraph as pg
from neuron import h
from cnmodel.protocols import Protocol
from cnmodel import cells
from cnmodel.util import sound
import cnmodel.util.pynrnutilities as PU
from cnmodel import data


def runTrial(cell, info):
    """
    info is a dict
    """
    if cell == 'bushy':
        post_cell = cells.Bushy.create(species='mouse')
    elif cell == 'tstellate':
        post_cell = cells.TStellate.create(species='mouse')
    elif cell == 'octopus':
        post_cell = cells.Octopus.create(species='mouse')
    elif cell == 'dstellate':
        post_dell = cells.DStellate.create(species='mouse')
    else:
        raise ValueError('cell %s is not yet implemented for PSTH testing' % self.cell)
    pre_cells = []
    synapses = []
    j = 0
    xmtr = {}
    for k in range(info['n_sgc']):
        preCell = cells.DummySGC(cf=info['cf'], sr=2)
        pre_cells.append(preCell)
        synapse = preCell.connect(post_cell)
        synapses.append(synapse)
        for i in range(synapse.terminal.n_rzones):
            xmtr['xmtr%03d'%j] = h.Vector()
            xmtr['xmtr%03d'%j].record(synapse.terminal.relsite._ref_XMTR[i])
            j = j + 1
        synapse.terminal.relsite.Dep_Flag = False  # no depression in these simulations
    seed = info['seed']
    Vm = h.Vector()
    for k in range(info['n_sgc']):
        seed = seed + 1  # independent inputs each trial and sgc
        pre_cells[k].set_sound_stim(info['stim'], seed=seed)
    Vm.record(post_cell.soma(0.5)._ref_v)
    #self['prevm'] = preCell.soma(0.5)._ref_v
    time = h._ref_t
    
    h.tstop = 1e3*info['run_duration'] # duration of a run
    h.celsius = info['temp']
    h.dt = info['dt']
    post_cell.cell_initialize()
    info['init']()
    h.run()
    return {'time': time, 'vm': Vm, 'xmtr': xmtr, 'pre_cells': pre_cells, 'post_cell': post_cell}

class SGCInputTestPSTH(Protocol):
    def set_cell(self, cell='bushy'):
        self.cell = cell
    
    def run(self, temp=34.0, dt=0.025, seed=575982035, reps=4, stimulus='tone', parallelize=False):
        assert stimulus in ['tone', 'SAM', 'clicks']  # cases available
        assert self.cell in ['bushy', 'tstellate', 'octopus', 'dstellate']
        self.nrep = reps
        self.stimulus = stimulus
        self.run_duration = 0.20  # in seconds
        self.pip_duration = 0.05  # in seconds
        self.pip_start = [0.01]  # in seconds
        self.Fs = 100e3  # in Hz
        self.f0 = 4000.  # stimulus in Hz
        self.cf = 4000.  # SGCs in Hz
        self.fMod = 100.  # mod freq, Hz
        self.dMod = 50.  # % mod depth, Hz
        self.dbspl = 60.
        if self.stimulus == 'SAM':
            self.stim = sound.SAMTone(rate=self.Fs, duration=self.run_duration, f0=self.f0, 
                          fmod=self.fMod, dmod=self.dMod, dbspl=self.dbspl,
                          ramp_duration=2.5e-3, pip_duration=self.pip_duration,
                          pip_start=self.pip_start)
        if self.stimulus == 'tone':
            self.f0 = 1000.
            self.cf = 1000.
            self.stim = sound.TonePip(rate=self.Fs, duration=self.run_duration, f0=self.f0, dbspl=self.dbspl,
                          ramp_duration=2.5e-3, pip_duration=self.pip_duration,
                          pip_start=self.pip_start)

        if self.stimulus == 'clicks':
            self.click_rate = 0.020  # msec
            self.stim = sound.ClickTrain(rate=self.Fs, duration=self.run_duration,
                        f0=self.f0, dbspl=self.dbspl, click_start=0.010, click_duration=100.e-6,
                        click_interval=self.click_rate, nclicks=int((self.run_duration-0.01)/self.click_rate),
                        ramp_duration=2.5e-3)
        
        n_sgc = data.get('convergence', species='mouse', post_type=self.cell, pre_type='sgc')[0]
        self.n_sgc = int(np.round(n_sgc))

        self.vms = []
        info = {'n_sgc': self.n_sgc, 'stim': self.stim, 'cf': self.cf,
            'seed': seed, 'run_duration': self.run_duration,
            'temp': temp, 'dt': dt, 'init': self.custom_init}
        if not parallelize:
            self.pre_cells = []
            self.synapses = []
            for n in range(self.nrep):
                res = runTrial(self.cell, info)
                # res contains: {'time': time, 'vm': Vm, 'xmtr': xmtr, 'pre_cells': pre_cells, 'post_cell': post_cell}
                if n == 0:
                    self['t'] = res['time']
                    # print res['xmtr']
                    # print {k: v.to_python() for k, v in res['xmtr'].items()}
                    self.xmtr = {k: v.to_python() for k, v in res['xmtr'].items()}
                    self.pre_cells = res['pre_cells']
                self.vms.append('vm')
                
        #     j = 0
        #     for k in range(self.n_sgc):
        #         preCell = cells.DummySGC(cf=self.cf, sr=2)
        #         self.pre_cells.append(preCell)
        #         synapse = preCell.connect(self.post_cell)
        #         self.synapses.append(synapse)
        #         for i in range(synapse.terminal.n_rzones):
        #             self['xmtr%03d'%j] = synapse.terminal.relsite._ref_XMTR[i]
        #             j = j + 1
        #         synapse.terminal.relsite.Dep_Flag = False  # no depression in these simulations
        #     for n in range(self.nrep):  # do repetitions
        #         for k in range(self.n_sgc):
        #             seed = seed + k  # independent inputs each trial and sgc
        #             self.pre_cells[k].set_sound_stim(self.stim, seed=seed)
        #         self['vm'] = self.post_cell.soma(0.5)._ref_v
        #         #self['prevm'] = preCell.soma(0.5)._ref_v
        #         self['t'] = h._ref_t
        #         self.post_cell.cell_initialize()
        #         h.tstop = 1e3*self.run_duration # duration of a run
        #         h.celsius = temp
        #         h.dt = dt
        #
        #         self.custom_init()
        #         h.run()
        #         self.vms.append(self['vm'])  # save the voltage (t will be the same for all runs)
        if parallelize:
            pass


    def show(self):
        self.win = pg.GraphicsWindow()
        self.win.setBackground('w')
        Fs = self.Fs
        p1 = self.win.addPlot(title='Stimulus', row=0, col=0,
            labels={'bottom': 'T (ms)', 'left': 'V'})
        p1.plot(self.stim.time * 1000, self.stim.sound, pen=pg.mkPen('k', width=0.75))
        p1.setXLink(p1)

        p2 = self.win.addPlot(title='AN spikes', row=1, col=0,
            labels={'bottom': 'T (ms)', 'left': 'AN spikes (first trial)'})
        xan = []
        yan = []
        for k in range(len(self.pre_cells)):
            r = self.pre_cells[k]._spiketrain
            xan.extend(r)
            yr = k + np.zeros_like(r) + 0.2
            yan.extend(yr)
        c = pg.PlotCurveItem()
        xp = np.repeat(np.array(xan), 2)
        yp = np.repeat(np.array(yan), 2)
        yp[1::2] = yp[::2] + 0.6
        c.setData(xp.flatten(), yp.flatten(), connect='pairs', pen=(pg.mkPen('k')))
        p2.addItem(c)
        p2.setXLink(p1)

        p3 = self.win.addPlot(title='%s Spikes' % self.cell, row=2, col=0,
            labels={'bottom': 'T (ms)', 'left': 'Trial #'})
        xcn = []
        ycn = []
        xspks = []
        for k in range(self.nrep):
            bspk = PU.findspikes(self['t'], self.vms[k], -30.)
            xcn.extend(bspk)
            yr = k + np.zeros_like(bspk) + 0.2
            ycn.extend(yr)
        d = pg.PlotCurveItem()
        xp = np.repeat(np.array(xcn), 2)
        yp = np.repeat(np.array(ycn), 2)
        yp[1::2] = yp[::2] + 0.6
        d.setData(xp.flatten(), yp.flatten(), connect='pairs', pen=pg.mkPen('k', width=1.5))
        p3.addItem(d)
        p3.setXLink(p1)

        p4 = self.win.addPlot(title='%s Vm' % self.cell, row=3, col=0,
            labels={'bottom': 'T (ms)', 'left': 'Vm (mV)'})
        p4.plot(self['t'], self['vm'], pen=pg.mkPen('k', width=1.0))
        p4.setXLink(p1)

        p5 = self.win.addPlot(title='xmtr', row=0, col=1,
            labels={'bottom': 'T (ms)', 'left': 'gSyn'})
        j = 0
        for k in range(self.n_sgc):
            synapse = self.synapses[k]
            for i in range(synapse.terminal.n_rzones):
                p5.plot(self['t'], self['xmtr%03d'%j], pen=(i, 15))
                j = j + 1
        p5.setXLink(p1)
        
        p6 = self.win.addPlot(title='AN PSTH', row=1, col=1,
            labels={'bottom': 'T (ms)', 'left': 'Sp/ms/trial'})
        bins = np.arange(0, 150, 1)
        (hist, binedges) = np.histogram(xan, bins)
        curve6 = p6.plot(binedges, hist, stepMode=True,
            fillBrush=(0, 0, 0, 255), brush=pg.mkBrush('k'), fillLevel=0)
        
        p7 = self.win.addPlot(title='%s PSTH' % self.cell, row=2, col=1,
            labels={'bottom': 'T (ms)', 'left': 'Sp/ms/trial'})
        bins = np.arange(0, 150, 1)
        (hist, binedges) = np.histogram(xcn, bins)
        curve7 = p7.plot(binedges, hist, stepMode=True,
            fillBrush=(0, 0, 0, 255), brush=pg.mkBrush('k'), fillLevel=0)
        
        
        # p6 = self.win.addPlot(title='AN phase', row=1, col=1)
        # phasewin = [self.pip_start[0] + 0.25*self.pip_duration, self.pip_start[0] + self.pip_duration]
        # prespk = self.pre_cells[0]._spiketrain  # just sample one...
        # spkin = prespk[np.where(prespk > phasewin[0]*1e3)]
        # spikesinwin = spkin[np.where(spkin <= phasewin[1]*1e3)]
        # print "\nCell type: %s" % self.cell
        # print "Stimulus: "
        #
        # # set freq for VS calculation
        # if self.stimulus == 'tone':
        #     f0 = self.f0
        #     print "Tone: f0=%.3f at %3.1f dbSPL, cell CF=%.3f" % (self.f0, self.dbspl, self.cf)
        # if self.stimulus == 'SAM':
        #     f0 = self.fMod
        #     print ("SAM Tone: f0=%.3f at %3.1f dbSPL, fMod=%3.1f  dMod=%5.2f, cell CF=%.3f" %
        #          (self.f0, self.dbspl, self.fMod, self.dMod, self.cf))
        # if self.stimulus == 'clicks':
        #     f0 = 1./self.click_rate
        #     print "Clicks: interval %.3f at %3.1f dbSPL, cell CF=%.3f " % (self.click_rate, self.dbspl, self.cf)
        # vs = PU.vector_strength(spikesinwin, f0)
        #
        # print 'AN Vector Strength: %7.3f, d=%.2f (us) Rayleigh: %7.3f  p = %.3e  n = %d' % (vs['r'], vs['d']*1e6, vs['R'], vs['p'], vs['n'])
        # (hist, binedges) = np.histogram(vs['ph'])
        # curve = p6.plot(binedges, hist, stepMode=True, fillBrush=(100, 100, 255, 150), fillLevel=0)
        # p6.setXRange(0., 2*np.pi)
        #
        # p7 = self.win.addPlot(title='%s phase' % self.cell, row=2, col=1)
        # spkin = bspk[np.where(bspk > phasewin[0]*1e3)]
        # spikesinwin = spkin[np.where(spkin <= phasewin[1]*1e3)]
        # vs = PU.vector_strength(spikesinwin, f0)
        # print '%s Vector Strength: %7.3f, d=%.2f (us) Rayleigh: %7.3f  p = %.3e  n = %d' % (self.cell, vs['r'], vs['d']*1e6, vs['R'], vs['p'], vs['n'])
        # (hist, binedges) = np.histogram(vs['ph'])
        # curve = p7.plot(binedges, hist, stepMode=True, fillBrush=(100, 100, 255, 150), fillLevel=0)
        # p7.setXRange(0., 2*np.pi)
        # p7.setXLink(p6)

        self.win.show()


if __name__ == '__main__':
    if len(sys.argv) > 1:
        cell = sys.argv[1]
    else:
        cell = 'bushy'
    if len(sys.argv) > 2:
        stimulus = sys.argv[2]
    else:
        stimulus = 'tone'
    print 'cell type: ', cell
    prot = SGCInputTestPSTH()
    prot.set_cell(cell)
    prot.run(stimulus=stimulus)
    prot.show()

    import sys
    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()
