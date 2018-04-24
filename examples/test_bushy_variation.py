"""
test_bushy_variation.py

Test inputs to bushy cells as we co-vary the KLT and IH conductances
Two tests: IV and spikes.
Simulation results are first written to disk (.p files); plotting is done separately

Usage:
python test_bushy_variation.py [a, b]
a runs IV curvers with variations of gKlt/gh
b runs PSTHs to AN input (CF tones) across the same variations.

"""


import sys
import numpy as np
import pyqtgraph as pg
import pyqtgraph.multiprocess as mp
import pickle
import time
from neuron import h
from cnmodel.protocols import Protocol
from cnmodel.protocols import iv_curve
from cnmodel import cells
from cnmodel.util import sound
from cnmodel.util import custom_init
from cnmodel.util import make_pulse
import cnmodel.util.pynrnutilities as PU
from cnmodel import data

import matplotlib.pyplot as mpl
import cnmodel.util.PlotHelpers as PH
import timeit




synapseType = 'multisite' # 'simple'
species = 'mouse'  # tables for other species do not yet exist

class RunTrial():
    def __init__(self, post_cell, info):
        """
        info is a dict
        """
        pre_cells = []
        synapses = []
        j = 0
        xmtr = {}
        for nsgc, sgc in enumerate(range(info['n_sgc'])):
            pre_cells.append(cells.DummySGC(cf=info['cf'], sr=info['sr']))
            if synapseType == 'simple':
                synapses.append(pre_cells[-1].connect(post_cell, type=synapseType))
                synapses[-1].terminal.netcon.weight[0] =info['gmax']*2.0
            elif synapseType == 'multisite':
                synapses.append(pre_cells[-1].connect(post_cell, post_opts={'AMPAScale': 2.0, 'NMDAScale': 2.0}, type=synapseType))
                for i in range(synapses[-1].terminal.n_rzones):
                    xmtr['xmtr%04d'%j] = h.Vector()
                    xmtr['xmtr%04d'%j].record(synapses[-1].terminal.relsite._ref_XMTR[i])
                    j = j + 1
                synapses[-1].terminal.relsite.Dep_Flag = False  # no depression in these simulations
            pre_cells[-1].set_sound_stim(info['stim'], seed = info['seed'] + nsgc, simulator=info['simulator'])
        Vm = h.Vector()
        Vm.record(post_cell.soma(0.5)._ref_v)
        rtime = h.Vector()
        rtime.record(h._ref_t)
        h.tstop = 1e3*info['run_duration'] # duration of a run
        h.celsius = info['temp']
        h.dt = info['dt']
        post_cell.cell_initialize()
        info['init']()
        h.t = 0.
        h.run()
        return {'time': np.array(rtime), 'vm': Vm.to_python(), 'xmtr': xmtr, 'pre_cells': pre_cells, 'post_cell': post_cell, 'synapses': synapses}

class SGCInputTestPSTH(Protocol):
    def set_cell(self, cell='bushy'):
        self.cell = cell
        self.parallelize = False
        
    def run(self, temp=34.0, dt=0.025, seed=575982035, reps=10, stimulus='tone', simulator='cochlea'):
        assert stimulus in ['tone', 'SAM', 'clicks']  # cases available
        assert self.cell in ['bushy', 'tstellate', 'octopus', 'dstellate']
        self.nrep = reps
        self.stimulus = stimulus
        self.run_duration = 0.20  # in seconds
        self.pip_duration = 0.05  # in seconds
        self.pip_start = [0.1]  # in seconds
        self.Fs = 100e3  # in Hz
        self.f0 = 4000.  # stimulus in Hz
        self.cf = 4000.  # SGCs in Hz
        self.fMod = 100.  # mod freq, Hz
        self.dMod = 0.  # % mod depth, Hz
        self.dbspl = 50.
        self.simulator = simulator
        self.sr = 1  # set SR group
        if self.stimulus == 'SAM':
            self.stim = sound.SAMTone(rate=self.Fs, duration=self.run_duration, f0=self.f0, 
                          fmod=self.fMod, dmod=self.dMod, dbspl=self.dbspl,
                          ramp_duration=2.5e-3, pip_duration=self.pip_duration,
                          pip_start=self.pip_start)
        if self.stimulus == 'tone':
            self.f0 = 4000.
            self.cf = 4000.
            self.stim = sound.TonePip(rate=self.Fs, duration=self.run_duration, f0=self.f0, dbspl=self.dbspl,
                          ramp_duration=2.5e-3, pip_duration=self.pip_duration,
                          pip_start=self.pip_start)

        if self.stimulus == 'clicks':
            self.click_rate = 0.020  # msec
            self.stim = sound.ClickTrain(rate=self.Fs, duration=self.run_duration,
                        f0=self.f0, dbspl=self.dbspl, click_start=0.010, click_duration=100.e-6,
                        click_interval=self.click_rate, nclicks=int((self.run_duration-0.01)/self.click_rate),
                        ramp_duration=2.5e-3)
        
        n_sgc = data.get('convergence', species=species, post_type=self.cell, pre_type='sgc')[0]
        self.n_sgc = int(np.round(n_sgc))
        # for simple synapses, need this value:
        self.AMPA_gmax = data.get('sgc_synapse', species=species,
                        post_type=self.cell, field='AMPA_gmax')[0]/1e3  # convert nS to uS for NEURON
        self.vms = [None for n in range(self.nrep)]
        self.synapses = [None for n in range(self.nrep)]
        self.xmtrs = [None for n in range(self.nrep)]
        self.pre_cells = [None for n in range(self.nrep)]
        self.time = [None for n in range(self.nrep)]
        info = {'n_sgc': self.n_sgc, 'gmax': self.AMPA_gmax, 'stim': self.stim,
                'simulator': self.simulator, 'cf': self.cf, 'sr': self.sr,
                'seed': seed, 'run_duration': self.run_duration,
                'temp': temp, 'dt': dt, 'init': custom_init}
        if not self.parallelize:
            for nr in range(self.nrep):
                info['seed'] = seed + 3*self.n_sgc*nr
                res = RunTrial(self.cell, info)
                # res contains: {'time': time, 'vm': Vm, 'xmtr': xmtr, 'pre_cells': pre_cells, 'post_cell': post_cell}
                self.pre_cells[nr] = res['pre_cells']
                self.time[nr] = res['time']
                self.xmtr = {k: v.to_python() for k, v in res['xmtr'].items()}
                self.vms[nr] = res['vm']
                self.synapses[nr] = res['synapses']
                self.xmtrs[nr] =self.xmtr

        if self.parallelize:
            ### Use parallelize with multiple workers
            tasks = range(len(self.nrep))
            results3 = results[:]
            start = time.time()
#            with mp.Parallelize(enumerate(tasks), results=results, progressDialog='processing in parallel..') as tasker:
            with mp.Parallelize(enumerate(tasks), results=results) as tasker:
                for i, x in tasker:
                    tot = 0
                    for j in xrange(size):
                        tot += j * x
                    tasker.results[i] = tot
            print( "\nParallel time, %d workers: %0.2f" % (mp.Parallelize.suggestedWorkerCount(), time.time() - start))
            print( "Results match serial:      %s" % str(results3 == results))

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
        for nr in range(self.nrep):
            xan = []
            yan = []
            for k in range(len(self.pre_cells[nr])):
                r = self.pre_cells[nr][k]._spiketrain
                xan.extend(r)
                yr = k + np.zeros_like(r) + 0.2
                yan.extend(yr)
            c = pg.PlotCurveItem()
            xp = np.repeat(np.array(xan), 2)
            yp = np.repeat(np.array(yan), 2)
            yp[1::2] = yp[::2] + 0.6
            c.setData(xp.flatten(), yp.flatten(), connect='pairs', pen=pg.mkPen(pg.intColor(nr, self.nrep), hues=self.nrep, width=1.0))
            p2.addItem(c)
        p2.setXLink(p1)

        p3 = self.win.addPlot(title='%s Spikes' % self.cell, row=2, col=0,
            labels={'bottom': 'T (ms)', 'left': 'Trial #'})
        xcn = []
        ycn = []
        xspks = []
        for k in range(self.nrep):
            bspk = PU.findspikes(self.time[k], self.vms[k], -35.)
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
        for nr in range(self.nrep):
            p4.plot(self.time[nr], self.vms[nr], pen=pg.mkPen(pg.intColor(nr, self.nrep), hues=self.nrep, width=1.0))
        p4.setXLink(p1)

        p5 = self.win.addPlot(title='xmtr', row=0, col=1,
            labels={'bottom': 'T (ms)', 'left': 'gSyn'})
        if synapseType == 'multisite':
            for nr in [0]:
                syn = self.synapses[nr]
                j = 0
                for k in range(self.n_sgc):
                    synapse = syn[k]
                    for i in range(synapse.terminal.n_rzones):
                        p5.plot(self.time[nr], self.xmtrs[nr]['xmtr%04d'%j], pen=pg.mkPen(pg.intColor(nr, self.nrep), hues=self.nrep, width=1.0))
                        j = j + 1
        p5.setXLink(p1)
        
        p6 = self.win.addPlot(title='AN PSTH', row=1, col=1,
            labels={'bottom': 'T (ms)', 'left': 'Sp/ms/trial'})
        bins = np.arange(0, 200, 1)
        (hist, binedges) = np.histogram(xan, bins)
        curve6 = p6.plot(binedges, hist, stepMode=True,
            fillBrush=(0, 0, 0, 255), brush=pg.mkBrush('k'), fillLevel=0)
        
        p7 = self.win.addPlot(title='%s PSTH' % self.cell, row=2, col=1,
            labels={'bottom': 'T (ms)', 'left': 'Sp/ms/trial'})
        bins = np.arange(0, 200, 1)
        (hist, binedges) = np.histogram(xcn, bins)
        curve7 = p7.plot(binedges, hist, stepMode=True,
            fillBrush=(0, 0, 0, 255), brush=pg.mkBrush('k'), fillLevel=0)

        self.win.show()

class Variations(Protocol):
    def __init__(self, runtype, runname):
        self.runtype = runtype
        self.runname = runname
        self.npost = 5  # number of post cells to test
        self.npre = 3  # number of presynaptic cells
        self.reset()
        
    def reset(self):
        super(Variations, self).reset()

    def make_cells(self, cf=16e3, temp=34.0, dt=0.025):
        self.pre_cells = []
        for n in range(self.npre):
            self.pre_cells.append(cells.DummySGC(cf=cf, sr=3))

        self.post_cells = []
        for n in range(self.npost):
            self.post_cell = cells.Bushy.create(species=species)
            self.post_cells.append(self.post_cell)
        for n in range(self.npre):
            for m in range(self.npost):
                synapse = self.pre_cells[n].connect(self.post_cells[m])
                self.synapse = synapse
                synapse.terminal.relsite.Dep_Flag = False
        
        # make variations in the postsynaptic cells
        varsg = [0.5, 0.75, 1.0, 1.5, 2.0]
        for i, m in enumerate(range(self.npost)):
            refgbar_klt = self.post_cells[m].soma().klt.gbar
            refgbar_ih = self.post_cells[m].soma().ihvcn.gbar
            self.post_cells[m].soma().klt.gbar = refgbar_klt * varsg[i]
            self.post_cells[m].soma().ihvcn.gbar = refgbar_ih * varsg[i]
            
        # self.stim = sound.TonePip(rate=100e3, duration=0.1, f0=4000, dbspl=80,
        #                           ramp_duration=2.5e-3, pip_duration=0.04,
        #                           pip_start=[0.02])
        #
        # preCell.set_sound_stim(self.stim, seed=seed)
        #
        # self['vm'] = postCell.soma(0.5)._ref_v
        # #self['prevm'] = preCell.soma(0.5)._ref_v
        # for i in range(30):
        #     self['xmtr%d'%i] = synapse.terminal.relsite._ref_XMTR[i]
        #     synapse.terminal.relsite.Dep_Flag = False

    def make_stimulus(self, stimulus='tone', cf=16000., f0=16000., simulator=None,
        rundur=0.2, pipdur=0.05, dbspl=50.,
        fmod=100., dmod=0.):
        self.stimulus = stimulus
        self.run_duration = rundur  # in seconds
        self.pip_duration = pipdur  # in seconds
        self.pip_start = [0.1]  # in seconds
        self.Fs = 100e3  # in Hz
        self.f0 = f0  # stimulus in Hz
        self.cf = cf  # SGCs in Hz
        self.fMod = fmod  # mod freq, Hz
        self.dMod = dmod  # % mod depth, Hz
        self.dbspl = dbspl
        self.simulator = simulator
        self.sr = 1  # set SR group
        if self.stimulus == 'SAM':
            self.stim = sound.SAMTone(rate=self.Fs, duration=self.run_duration, f0=self.f0, 
                          fmod=self.fMod, dmod=self.dMod, dbspl=self.dbspl,
                          ramp_duration=2.5e-3, pip_duration=self.pip_duration,
                          pip_start=self.pip_start)
        if self.stimulus == 'tone':
            self.stim = sound.TonePip(rate=self.Fs, duration=self.run_duration, f0=self.f0, dbspl=self.dbspl,
                          ramp_duration=2.5e-3, pip_duration=self.pip_duration,
                          pip_start=self.pip_start)

        if self.stimulus == 'clicks':
            self.click_rate = 0.020  # msec
            self.stim = sound.ClickTrain(rate=self.Fs, duration=self.run_duration,
                        f0=self.f0, dbspl=self.dbspl, click_start=0.010, click_duration=100.e-6,
                        click_interval=self.click_rate, nclicks=int((self.run_duration-0.01)/self.click_rate),
                        ramp_duration=2.5e-3)
        

    def run(self, mode='IV', cf=16e3, temp=34.0, dt=0.025, stimamp=0, iinj = [
        1.0]):
        self.dt = dt
        self.temp = temp
        
        self.make_cells(cf, temp, dt)
        print dir(self.pre_cells[0])
        seed = 0
        j = 0
        if mode == 'sound':
            self.make_stimulus(stimulus='tone')

            for np in range(len(self.pre_cells)):
                self.pre_cells[np].set_sound_stim(self.stim, seed=seed)
                seed += 1
                synapses.append(pre_cells[-1].connect(post_cell,
                    post_opts={'AMPAScale': 2.0, 'NMDAScale': 2.0}, type=synapseType))
                for i in range(synapses[-1].terminal.n_rzones):
                    xmtr['xmtr%04d'%j] = h.Vector()
                    xmtr['xmtr%04d'%j].record(synapses[-1].terminal.relsite._ref_XMTR[i])
                j = j + 1
                synapses[-1].terminal.relsite.Dep_Flag = False  # no depression in these simulations
                
        print 'setup to run'
        self.stim_params = []
        self.istim = []
        self.istims = []
        if mode == 'pulses':
            for i, pre_cell in enumerate(self.pre_cells):
                stim = {}
                stim['NP'] = 10
                stim['Sfreq'] = 100.0 # stimulus frequency
                stim['delay'] = 10.0
                stim['dur'] = 0.5
                stim['amp'] = stimamp
                stim['PT'] = 0.0
                stim['dt'] = dt
                (secmd, maxt, tstims) = make_pulse(stim)
                self.stim_params.append(stim)
                stim['amp'] = 0.
                (secmd2, maxt2, tstims2) = make_pulse(stim)
            
                # istim current pulse train
                istim = h.iStim(0.5, sec=pre_cell.soma)
                i_stim_vec = h.Vector(secmd)
                self.istim.append((istim, i_stim_vec))
                self['istim%02d' % i] = istim._ref_i
        
        self.postpars = []
        self.poststims = []
        if mode == 'IV':
            for i, post_cell in enumerate(self.post_cells):
                pstim = {}
                pstim['NP'] = 1
                pstim['Sfreq'] = 100.0 # stimulus frequency
                pstim['delay'] = 10.0
                pstim['dur'] = 100
                pstim['amp'] = iinj[0]
                pstim['PT'] = 0.0
                pstim['dt'] = dt
                (secmd, maxt, tstims) = make_pulse(pstim)
                self.postpars.append(pstim)
            
                # istim current pulse train
                pistim = h.iStim(0.5, sec=post_cell.soma)
                pi_stim_vec = h.Vector(secmd)
                self.poststims.append((pistim, pi_stim_vec))
                self['poststim%02d' % i] = pistim._ref_i            
        #
        # Run simulation
        #
        h.celsius = temp
        self.temp = temp
        
        # first find rmp for each cell
        for m in range(self.npost):
            self.post_cells[m].vm0 = None
            self.post_cells[m].cell_initialize()
            # set starting voltage...
            self.post_cells[m].soma(0.5).v = self.post_cells[m].vm0
        h.dt = 0.02
        h.t = 0  # run a bit to find true stable rmp
        h.tstop = 20.
        h.batch_save()
        h.batch_run(h.tstop, h.dt, 'v.dat')
        
        # order matters: don't set these up until we need to
        self['t'] = h._ref_t
        # set up recordings
        if mode == 'pulses':
            for i in range(self.npre):
                self.istim[i][1].play(self.istim[i][0]._ref_i, dt, 0)
                self['v_pre%02d'%i] = self.pre_cells[i].soma(0.5)._ref_v
        for m in range(self.npost):
            if mode == 'IV':
                self.poststims[m][1].play(self.poststims[m][0]._ref_i, dt, 0)
            self['v_post%02d' % m] = self.post_cells[m].soma(0.5)._ref_v
        h.finitialize()  # init and instantiate recordings
        print 'running'
        h.t = 0.
        h.tstop = 200.
        h.batch_run(h.tstop, h.dt, 'v.dat')
        # while h.t < h.tstop:  # get data (do not use h.run() - try it and see why!)
        #      h.fadvance()

    def runIV(self, parallelize):
        self.civ = {}
        self.iiv = []
        varsg = np.linspace(0.25, 2.0, int((2.0-0.25)/0.25)+1) #[0.5, 0.75, 1.0, 1.5, 2.0]  # covary Ih and gklt in constant ratio
        self.gklts = np.zeros(len(varsg))
        self.ghs = np.zeros(len(varsg))
        if not parallelize:
            for n in range(self.npost):
                self.civ[n] = []
    #            self.iiv[c] = []
            start = time.time()
            for inj in np.arange(-1.0, 1.51, 0.5):
                self.run(mode='IV', temp=34.0, dt=0.025, stimamp=10, iinj=[inj])
                print 'ran for current = ', inj
                for c in range(self.npost):
                    self.civ[c].append(self['v_post%02d' % c])
                    if c == 0:  # just the first
                        self.iiv.append(self['poststim%02d' %  c])
            print( "\nSerial time, %0.2f" % (time.time() - start))
            if runname is not None:
                f = open(runname, 'w')
                pickle.dump({'t': self['t'], 'v': self.civ, 'i': self.iiv}, f)
                f.close()
        else:
            # mp.parallelizer.multiprocessing.cpu_count()
            nworker = 16
            self.npost = len(varsg)
            tasks = range(self.npost)
            results = [None] * len(tasks)
            ivc = [None] * len(tasks)
            start = time.time()
#            with mp.Parallelize(enumerate(tasks), results=results, workers=nworker, progressDialog='processing in parallel..') as tasker:
            with mp.Parallelize(enumerate(tasks), results=results, workers=nworker) as tasker:
                for i, x in tasker:
                    post_cell = cells.Bushy.create(species=species)
                    refgbar_klt = post_cell.soma().klt.gbar
                    refgbar_ih = post_cell.soma().ihvcn.gbar
                    gklts = refgbar_klt * varsg[i]
                    ghs = refgbar_ih * varsg[i]
                    post_cell.soma().klt.gbar = gklts
                    post_cell.soma().ihvcn.gbar = ghs
                    post_cell.initial_mechanisms = None  # forget the mechanisms we set up initially
                    post_cell.save_all_mechs()  # and save new ones because we are explicitely varying them
                    ivc[i] = iv_curve.IVCurve()
                    ivc[i].run({'pulse': [(-1., 1.5, 0.25)]}, post_cell)
                    tasker.results[i] = {'v': ivc[i].voltage_traces, 'i': ivc[i].current_traces, 't': ivc[i].time_values, 'gklt': gklts, 'gh': ghs}
            print( "\nParallel time: %d workers,  %0.2f sec" % (nworker, time.time() - start))
            cell_info = {'varrange': varsg}
            print cell_info
            res = {'cells': cell_info, 'results': results}
            if runname is not None:
                f = open(runname, 'w')
                pickle.dump(res, f)
                f.close()
                

    def runSound(self, parallelize=False):
        self.civ = {}
        self.iiv = []
        if not parallelize:
            pass

        if parallelize:
            nworker = 16
            varsg = np.linspace(0.25, 2.0, int((2.0-0.25)/0.25)+1) #[0.5, 0.75, 1.0, 1.5, 2.0]  # covary Ih and gklt in constant ratio
            self.npost = len(varsg)
            nrep = 25
            tasks = range(self.npost)
            results = [None] * len(tasks)
            ivc = [None] * len(tasks)
            gklts = np.zeros(len(varsg))
            ghs = np.zeros(len(varsg))
            start = time.time()
            seed = 0
            cf = 16000.
            f0 = 16000.
            rundur = 0.25 # seconds
            pipdur = 0.1 # seconds
            dbspl = 50.
            fmod = 40.
            dmod = 0.
            stimulus = 'tone'
#            with mp.Parallelize(enumerate(tasks), results=results, workers=nworker, progressDialog='processing in parallel..') as tasker:
            with mp.Parallelize(enumerate(tasks), results=results, workers=nworker) as tasker:
                for i, x in tasker:
                    post_cell = cells.Bushy.create(species=species)
                    h.celsius = 34
                    self.temp = h.celsius
                    refgbar_klt = post_cell.soma().klt.gbar
                    refgbar_ih = post_cell.soma().ihvcn.gbar
                    gklts[i] = refgbar_klt * varsg[i]
                    ghs[i] = refgbar_ih * varsg[i]
                    post_cell.soma().klt.gbar = gklts[i]
                    post_cell.soma().ihvcn.gbar = ghs[i]
                    post_cell.initial_mechanisms = None  # forget the mechanisms we set up initially
                    post_cell.save_all_mechs()  # and save new ones because we are explicitely varying them
                    self.make_stimulus(stimulus=stimulus, cf=cf, f0=f0, rundur=rundur, pipdur=pipdur, 
                        dbspl=50., simulator=None, fmod=fmod, dmod=dmod)
                    
                    pre_cells = []
                    synapses = []
                    for n in range(self.npre):
                        pre_cells.append(cells.DummySGC(cf=cf, sr=2))
                        synapses.append(pre_cells[n].connect(post_cell, type=synapseType))
                    v_reps = []
                    i_reps = []
                    p_reps = [] # pre spike on 0'th sgc
                    for j in range(nrep):
                        for prec in range(len(pre_cells)):
                            pre_cells[prec].set_sound_stim(self.stim, seed=seed)
                            seed += 1
                            # for i in range(synapses[-1].terminal.n_rzones):
                            #     xmtr['xmtr%04d'%j] = h.Vector()
                            #     xmtr['xmtr%04d'%j].record(synapses[-1].terminal.relsite._ref_XMTR[i])
                            # j = j + 1
                            #synapses[-1].terminal.relsite.Dep_Flag = False  # no depression in these simulations
                        #
                        # Run simulation
                        #
                        post_cell.vm0 = None
                        post_cell.cell_initialize()
                        # set starting voltage...
                        post_cell.soma(0.5).v = post_cell.vm0
                        h.dt = 0.02
                        h.t = 0  # run a bit to find true stable rmp
                        h.tstop = 20.
                        h.batch_save()
                        h.batch_run(h.tstop, h.dt, 'v.dat')
                        self['t'] = h._ref_t
                        # set up recordings
                        self['v_post%02d' % j] = post_cell.soma(0.5)._ref_v
                        h.finitialize()  # init and instantiate recordings
                        print 'running %d' % i
                        h.t = 0.
                        h.tstop = rundur*1000.  # rundur is in seconds.
                        post_cell.check_all_mechs()  # make sure no further changes were introduced before run.
                        h.batch_run(h.tstop, h.dt, 'v.dat')
                        v_reps.append(self['v_post%02d' % j])
                        i_reps.append(0.*self['v_post%02d' % j])
                        p_reps.append(pre_cells[0]._stvec.to_python())
                    tasker.results[i] = {'v': v_reps, 'i': i_reps, 't': self['t'], 'pre': pre_cells[0]._stvec.to_python()}
            print( "\nParallel time: %d workers,  %0.2f sec" % (nworker, time.time() - start))
            cell_info = {'gklt': gklts, 'gh': ghs}
            stim_info = {'nreps': nrep, 'cf': cf, 'f0': f0, 'rundur': rundur, 'pipdur': pipdur, 'dbspl': dbspl, 'fmod': fmod, 'dmod': dmod}
            res = {'cells': cell_info, 'stim': stim_info, 'results': results}
            if runname is not None:
                f = open(runname, 'w')
                pickle.dump(res, f)
                f.close()

    #
    # def show(self):
    #
    #     self.win = pg.GraphicsWindow()
    #     self.win.resize(1000, 1000)
    #
    #     cmd_plot = self.win.addPlot(title='Stim')
    #     for i in range(len(self.pre_cells)):
    #         cmd_plot.plot(self['t'], self['istim%02d' % i], pen=pg.mkPen(pg.intColor(i, len(self.pre_cells)), hues=len(self.pre_cells), width=1.0))
    #
    #     self.win.nextRow()
    #     pre_plot = self.win.addPlot(title='SGC Vm')
    #     for i in range(len(self.pre_cells)):
    #         pre_plot.plot(self['t'], self['v_pre%02d'%i], pen=pg.mkPen(pg.intColor(i, len(self.pre_cells)), hues=len(self.pre_cells), width=1.0))
    #
    #     self.win.nextRow()
    #     post_plot = self.win.addPlot(title='Post Cell: %s' % self.post_cell.type)
    #     for m in range(self.npost):
    #         post_plot.plot(self['t'], self['v_post%02d' % m],
    #             pen=pg.mkPen(pg.intColor(m, len(self.post_cells)), hues=len(self.post_cells), width=1.0))

def showpicklediv(name):
    f = open(name, 'r')
    result = pickle.load(f)
    f.close()
    d = result['results']
    ncells = len(d)
    vr = result['cells']['varrange']
    fig, ax = mpl.subplots(ncells+1, 2, figsize=(8.5, 11.))
#    fig.set_size_inches(8.5, 11., forward=True)
    for ni in range(len(d[0]['i'])):
        ax[-1, 0].plot(d[0]['t'], d[0]['i'][ni], 'k', linewidth=0.5)
    ax[-1, 0].set_ylim([-2., 2.])
    for nc in range(ncells):
        for ni in range(len(d[nc]['v'])):
            ax[nc, 0].plot(d[nc]['t'], d[nc]['v'][ni], 'k', linewidth=0.5)
            ax[nc, 0].set_ylim([-180., 40.])
            if ni == 0:
                ax[nc, 0].annotate( '%.2f' % vr[nc], (180., 20.))
    PH.nice_plot(ax.ravel().tolist())
    PH.noaxes(ax.ravel()[:-1].tolist())
    PH.calbar(ax[0, 0], calbar=[120., -130., 25., 50.], axesoff=True, orient='left',
        unitNames={'x': 'ms', 'y': 'mV'}, fontsize=9, weight='normal', font='Arial')
    PH.calbar(ax[-1, 0], calbar=[120., 0.5, 25., 1.], axesoff=True, orient='left',
        unitNames={'x': 'ms', 'y': 'nA'}, fontsize=9, weight='normal', font='Arial')
            
    mpl.show()

def vector_plot(f, r, l, yp=None):
    ax2 = f.add_axes([yp.x1-0.06, yp.y1-0.06, 0.05, 0.05], polar=True)
    r = np.repeat(r, 3)
    l = np.repeat(l, 3)
    for i in range(2, len(l), 3):
        l[i] = 0.
        l[i-2] = 0.
    ax2.plot(r, l, lw=0.5)
    #ax2.arrow(0, 0, np.mean(l), np.mean(r), head_width=0.05, head_length=-0.1, fc='r', ec='r')

def phase_hist(f, spkphase, yp=None):
    ax2 = f.add_axes([yp.x1-0.06, yp.y1-0.06, 0.05, 0.05])
    n, bins = np.histogram(spkphase, np.linspace(0., 2*np.pi, 91.), density=False)
    ax2.bar(bins[:-1], n, width = bins[1], facecolor='k', alpha=0.75)
    
def clean_spiketimes(spikeTimes, mindT=0.7):
    """
    Clean up spike time array, removing all less than mindT
    spikeTimes is a 1-D list or array
    mindT is difference in time, same units as spikeTimes
    If 1 or 0 spikes in array, just return the array
    """
    if len(spikeTimes) > 1:
        dst = np.diff(spikeTimes)
        st = np.array(spikeTimes[0])  # get first spike
        sok = np.where(dst > mindT)
        st = np.append(st, [spikeTimes[s+1] for s in sok])
        # print st
        spikeTimes = st
    return spikeTimes

def showplots(name):
    """
    Show traces from sound stimulation - without current injection
    """
    f = open(name, 'r')
    d = pickle.load(f)
    f.close()
    ncells = len(d['results'])
    stiminfo = d['stim']
    dur = stiminfo['rundur']*1000.
    print 'dur: ', dur
    print 'stim info: '
    print '  fmod: ', stiminfo['fmod']
    print '  dmod: ', stiminfo['dmod']
    print '  f0:   ', stiminfo['f0']
    print '  cf:   ', stiminfo['cf']
    varsg = np.linspace(0.25, 2.0, int((2.0-0.25)/0.25)+1)  # was not stored... 
    fig, ax = mpl.subplots(ncells+1, 2, figsize=(8.5, 11.))
    spikelists = [[]]*ncells
    prespikes = [[]]*ncells
    xmin = 50.
    for i in range(ncells):
        vdat = d['results'][i]['v']
        idat = d['results'][i]['i']
        tdat = d['results'][i]['t']
        pdat = d['results'][i]['pre']
        PH.noaxes(ax[i, 0])
        # if i == 0:
        #     PH.calbar(ax[0, 0], calbar=[120., -120., 25., 20.], axesoff=True, orient='left',
        #         unitNames={'x': 'ms', 'y': 'mV'}, fontsize=9, weight='normal', font='Arial')
        for j in range(len(vdat)):
            if j == 2:
                ax[i, 0].plot(tdat-xmin, vdat[j], 'k', linewidth=0.5)
            if j == 0:
                ax[i, 0].annotate( '%.2f' % varsg[i], (180., 20.))
            
        ax[i, 0].set_xlim([0, dur-xmin])
        ax[i, 0].set_ylim([-75, 50])
        PH.referenceline(ax[i, 0], reference=-62.0, limits=None, color='0.33', linestyle='--' ,linewidth=0.5, dashes=[3, 3])

        for j in range(len(vdat)):
            detected = PU.findspikes(tdat, vdat[j], -20.)
            detected = clean_spiketimes(detected)
            spikelists[i].extend(detected)
            if j == 0:
                n, bins = np.histogram(detected, np.linspace(0., dur, 201), density=False)
            else:
                m, bins = np.histogram(detected, np.linspace(0., dur, 201), density=False)
                n += m
        prespikes[i].extend(pdat)
        if j == 0:
            n, bins = np.histogram(pdat, np.linspace(0., dur, 201), density=False)
        else:
            m, bins = np.histogram(pdat, np.linspace(0., dur, 201), density=False)
            n += m

        ax[i, 1].bar(bins[:-1]-xmin, n, width = bins[1], facecolor='k', alpha=0.75)
        ax[i, 1].set_xlim([0, dur-xmin])
        ax[i, 1].set_ylim([0, 30])
        vs = PU.vector_strength(spikelists[i], stiminfo['fmod'])
        pre_vs = PU.vector_strength(prespikes[i], stiminfo['fmod'])
        # print 'pre: ', pre_vs
        # print 'post: ', vs
#         apos = ax[i,1].get_position()
#         ax[i, 1].set_title('VS = %4.3f' % pre_vs['r'])
# #        vector_plot(fig, vs['ph'], np.ones(len(vs['ph'])), yp = apos)
#         phase_hist(fig, vs['ph'], yp=apos)
#        phase_hist(fig, pre_vs['ph'], yp=apos)
    prot = Variations(runtype, runname)
#   stim_info = {'nreps': nrep, 'cf': cf, 'f0': f0, 'rundur': rundur, 'pipdur': pipdur, 'dbspl': dbspl, 'fmod': fmod, 'dmod': dmod}
    if stiminfo['dmod'] > 0:
        stimulus = 'SAM'
    else:
        stimulus = 'tone'
    prot.make_stimulus(stimulus=stimulus, cf=stiminfo['cf'], f0=stiminfo['f0'], simulator=None,
            rundur=stiminfo['rundur'], pipdur = stiminfo['pipdur'], dbspl=stiminfo['dbspl'],
            fmod=stiminfo['fmod'], dmod=stiminfo['dmod'])
    ax[-1, 1].plot(prot.stim.time*1000.-xmin, prot.stim.sound, 'k-', linewidth=0.75)
    ax[-1, 1].set_xlim([0, (dur-xmin)])

    PH.noaxes(ax[-1, 0])
    ax[-1, 0].set_xlim([0, dur-xmin])
    ax[-1, 0].set_ylim([-75, 50])
#    PH.referenceline(ax[-1, 0], reference=-62.0, limits=None, color='0.33', linestyle='--' ,linewidth=0.5, dashes=[3, 3])
    PH.calbar(ax[-1, 0], calbar=[20., 0., 25., 20.], axesoff=True, orient='left',
        unitNames={'x': 'ms', 'y': 'mV'}, fontsize=9, weight='normal', font='Arial')

    PH.cleanAxes(ax.ravel().tolist())
    
    mpl.show()
    
if __name__ == '__main__':
    runname = None
    panel = None
    if len(sys.argv) == 2:
        panel = sys.argv[1]
    if panel == 'a':
        runtype = 'IV'
        runname = 'Figure6_IV'
    elif panel == 'b':
        runtype = 'sound'
        runname = 'Figure6_AN'
    else:
        runtype = panel
    if panel is None:
        raise ValueError("Must specify figure panel to generate: 'a', 'b'")
    if runtype in ['sound', 'IV']:
         prot = Variations(runtype, runname)
         if runtype == 'IV':
             start_time = timeit.default_timer()
             prot.runIV(parallelize = True)
             elapsed = timeit.default_timer() - start_time
             print ('Elapsed time for IV simulations: %f' % (elapsed))
             showpicklediv(runname)
         if runtype == 'sound':
             start_time = timeit.default_timer()
             prot.runSound(parallelize=True)
             elapsed = timeit.default_timer() - start_time
             print ('Elapsed time for AN simulations: %f' % (elapsed))
             showplots(runname)
#         pg.show()
         # if sys.flags.interactive == 0:
         #    pg.QtGui.QApplication.exec_()
            
    elif runtype in ['showiv']:
        showpicklediv(runname)

    elif runtype in ['plots']:
        showplots(runname)
    else:
        print 'run type should be one of sound, IV, showiv, plots'
        

