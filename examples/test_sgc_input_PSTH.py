import sys
import argparse
import numpy as np
import timeit
# import dask
import cloudpickle
import pyqtgraph as pg
import pyqtgraph.multiprocess as mp

from neuron import h
from cnmodel.protocols import Protocol
from cnmodel import cells
from cnmodel.util import sound
from cnmodel.util import custom_init
import cnmodel.util.pynrnutilities as PU
from cnmodel import data

species = "mouse"  # tables for other species do not yet exist



def buildCell(cell):
    """
    Make the postsynaptic cell
    
    Parameters
    ----------
    cell : str
        Name of a cell type to instantiate
    
    Return
    ------
    cnmodel postsynaptic cell object
    """
    if cell == "bushy":
        post_cell = cells.Bushy.create(species=species)
    elif cell == "tstellate":
        post_cell = cells.TStellate.create(species=species)
    elif cell == "octopus":
        post_cell = cells.Octopus.create(species=species)
    elif cell == "dstellate":
        post_cell = cells.DStellate.create(species=species)
    elif cell == "tuberculoventral":
        post_cell = cells.DStellate.create(species=species)
    elif cell == "pyramidal":
        post_cell = cells.DStellate.create(species=species)
    else:
        raise ValueError("cell %s is not yet implemented for PSTH testing" % self.cell)
    print('made postcell')
    return post_cell

def makeSynapses(post_cell, info):
    """
    Attach synapses to the postsynaptic cell by the rules in info
    
    Parameters
    post_cell : cnmodel cell object
    
    info: dict 
    
    Returns
    -------
    dict of:
        pre_cells (list of cnmodel dummy SGC cells)
        synapses (list of snmodel synapses) 
        and the xmtr variable (hoc)
    """
    pre_cells = []
    synapses = []
    j = 0
    xmtr = {}
    for nsgc, sgc in enumerate(range(info["n_sgc"])):
        pre_cells.append(cells.DummySGC(cf=info["cf"], sr=info["sr"]))
        if synapseType == "simple":
            synapses.append(pre_cells[-1].connect(post_cell, type=synapseType))
            synapses[-1].terminal.netcon.weight[0] = info["gmax"]
        elif synapseType == "multisite":
            synapses.append(
                pre_cells[-1].connect(
                    post_cell,
                    post_opts={"AMPAScale": 1.0, "NMDAScale": 1.0},
                    type=synapseType,
                )
            )
            for i in range(synapses[-1].terminal.n_rzones):
                xmtr["xmtr%04d" % j] = h.Vector()
                xmtr["xmtr%04d" % j].record(synapses[-1].terminal.relsite._ref_XMTR[i])
                j = j + 1
            synapses[
                -1
            ].terminal.relsite.Dep_Flag = False  # no depression in these simulations
        pre_cells[-1].set_sound_stim(
            info["stim"], seed=info["seed"] + nsgc, simulator=info["simulator"]
        )
    print('makesyn ok')
    return {'pre_cells':pre_cells, 'synapses': synapses, 'xmtr':xmtr}

# def setupNeuron(post_cell, info):
#     Vm = h.Vector()
#     Vm.record(post_cell.soma(0.5)._ref_v)
#     rtime = h.Vector()
#     rtime.record(h._ref_t)
#     h.tstop = 1e3 * info["run_duration"]  # duration of a run
#     h.celsius = info["temp"]
#     h.dt = info["dt"]
#     post_cell.cell_initialize()
#     # info["init"]()
#     print('setup ok')
#     return {'rtime': rtime, 'Vm': Vm}
    
def runNeuron(post_cell, info, psx):

    print('runNeuron entry')
    Vm = h.Vector()
    Vm.record(post_cell.soma(0.5)._ref_v)
    rtime = h.Vector()
    rtime.record(h._ref_t)
    h.tstop = 1e3 * info["run_duration"]  # duration of a run
    h.celsius = info["temp"]
    h.dt = info["dt"]
    post_cell.cell_initialize()
    info["init"]()
    h.t = 0.0
    h.run()
    print('runNeuron exit')
    return {
        "time": np.array(rtime),
        "vm": Vm.to_python(),
        "xmtr": psx['xmtr'],
        "pre_cells": psx['pre_cells'],
        "post_cell": post_cell,
        "synapses": psx['synapses'],
    }

def runTrial(cell, info):
    post_cell = buildCell(cell)
    psx = makeSynapses(post_cell, info )
    # setupNeuron(post_cell, info)
    res = runNeuron(post_cell, info, psx)
    return(res)

# # dask delayed routines for parallelism
# @dask.delayed
# def buildCell_dask(cell):
#     return buildCell(cell)
#
# @dask.delayed
# def makeSynapses_dask(post_cell, info):
#     return makeSynapses(post_cell, info)
#
# @dask.delayed
# def runNeuron_dask(post_cell, info,  psx): # , rtime, Vm):
#     return runNeuron(post_cell, info, psx)
    # Vm = h.Vector()
    # Vm.record(post_cell.soma(0.5)._ref_v)
    # rtime = h.Vector()
    # rtime.record(h._ref_t)
    # h.tstop = 1e3 * info["run_duration"]  # duration of a run
    # h.celsius = info["temp"]
    # h.dt = info["dt"]
    # post_cell.cell_initialize()
    # info["init"]()
    # h.t = 0.0
    # h.run()
    # return {
    #     "time": np.array(rtime),
    #     "vm": np.array(Vm),
    #     "xmtr": psx['xmtr'],
    #     "pre_cells": psx['pre_cells'],
    #     "post_cell": post_cell,
    #     "synapses": psx['synapses'],
    # }
    

# @dask.delayed
# def setupNeuron_dask(post_cell, info):
#     return setupNeuron(post_cell, info)

class SGCInputTestPSTH(Protocol):
    def __init__(self):
        self.nrep = 10
        self.stimulus = "tone"
        self.run_duration = 0.20  # in seconds
        self.pip_duration = 0.05  # in seconds
        self.pip_start = [0.1]  # in seconds
        self.Fs = 100e3  # in Hz
        self.f0 = 4000.0  # stimulus in Hz
        self.cf = 4000.0  # SGCs in Hz
        self.fMod = 100.0  # mod freq, Hz
        self.dMod = 0.0  # % mod depth, Hz
        self.dbspl = 50.0
        self.simulator = "cochlea"
        self.sr = 2  # set SR group

    def set_cell(self, cell="bushy"):
        self.cell = cell

    def set_freq(self, fr=4000.0):
        self.f0 = fr
    
    def set_cf(self, cf=4000.0):
        self.cf = fr
    
    def set_db(self, db=50.):
        self.dbspl = db
    
    def set_sr(self, srgroup=1):
        self.sr = srgroup
        
    def run(
        self,
        temp=34.0,
        dt=0.025,
        seed=575982035,
        reps=10,
        stimulus="tone",
        simulator="cochlea",
        parallelize=False,
    ):
        assert stimulus in ["tone", "SAM", "clicks"]  # cases available
        assert self.cell in [
            "bushy",
            "tstellate",
            "octopus",
            "dstellate",
            "tuberculoventral",
            "pyramidal",
        ]
        self.nrep = reps
        if self.stimulus == "SAM":
            self.dMod = 100.0
            self.stim = sound.SAMTone(
                rate=self.Fs,
                duration=self.run_duration,
                f0=self.f0,
                fmod=self.fMod,
                dmod=self.dMod,
                dbspl=self.dbspl,
                ramp_duration=2.5e-3,
                pip_duration=self.pip_duration,
                pip_start=self.pip_start,
            )
        if self.stimulus == "tone":
            self.f0 = 4000.0
            self.cf = 4000.0
            self.stim = sound.TonePip(
                rate=self.Fs,
                duration=self.run_duration,
                f0=self.f0,
                dbspl=self.dbspl,
                ramp_duration=2.5e-3,
                pip_duration=self.pip_duration,
                pip_start=self.pip_start,
            )

        if self.stimulus == "clicks":
            self.click_rate = 0.020  # msec
            self.stim = sound.ClickTrain(
                rate=self.Fs,
                duration=self.run_duration,
                f0=self.f0,
                dbspl=self.dbspl,
                click_start=0.010,
                click_duration=100.0e-6,
                click_interval=self.click_rate,
                nclicks=int((self.run_duration - 0.01) / self.click_rate),
                ramp_duration=2.5e-3,
            )

        n_sgc = data.get(
            "convergence", species=species, post_type=self.cell, pre_type="sgc"
        )[0]
        self.n_sgc = int(np.round(n_sgc))
        # for simple synapses, need this value:
        self.AMPA_gmax = (
            data.get(
                "sgc_synapse", species=species, post_type=self.cell, field="AMPA_gmax"
            )[0]
            / 1e3
        )  # convert nS to uS for NEURON
        self.vms = [None for n in range(self.nrep)]
        self.synapses = [None for n in range(self.nrep)]
        self.xmtrs = [None for n in range(self.nrep)]
        self.pre_cells = [None for n in range(self.nrep)]
        self.time = [None for n in range(self.nrep)]
        info = {
            "n_sgc": self.n_sgc,
            "gmax": self.AMPA_gmax,
            "stim": self.stim,
            "simulator": self.simulator,
            "cf": self.cf,
            "sr": self.sr,
            "seed": seed,
            "run_duration": self.run_duration,
            "temp": temp,
            "dt": dt,
            "init": custom_init,
        }
        if not parallelize:
            start_time = timeit.default_timer()
            for nr in range(self.nrep):
                info["seed"] = seed + 3 * self.n_sgc * nr
                res = runTrial(self.cell, info)
                # res contains: {'time': time, 'vm': Vm, 'xmtr': xmtr, 'pre_cells': pre_cells, 'post_cell': post_cell}
                self.pre_cells[nr] = res["pre_cells"]
                self.time[nr] = res["time"]
                self.xmtr = {k: v.to_python() for k, v in list(res["xmtr"].items())}
                self.vms[nr] = res["vm"]
                self.synapses[nr] = res["synapses"]
                self.xmtrs[nr] = self.xmtr
            elapsed = timeit.default_timer() - start_time
            print(f"Not Parallel Elapsed time for {self.nrep:d} stimuli: {elapsed:f} secs")

        if parallelize:

            results = {}
            workers = 1 if not parallelize else None
            tot_runs = self.nrep
            tasks = []
            for nr in range(self.nrep):
                tasks.append(nr)
            start_time = timeit.default_timer()
            with mp.Parallelize(enumerate(tasks), results=results, workers=workers) as tasker:
                for i, task in tasker:
                    iteration = task
                    repres = runTrial(self.cell, info)
                    tasker.results[i] = (stim, result)
                    print(f"  --- finished run {i+1:5d}/{tot_runs:5d} ---" )
        
            # get time of run before display
            elapsed = timeit.default_timer() - start_time
            print('Elapsed time for %d stimuli: %f  (%f sec per stim), synapses: %s' % (len(tasks), elapsed, elapsed/len(tasks), prot.bushy._synapsetype))


            # import dask.multiprocessing
            # dask.config.set(scheduler='processes')
            # # from dask.distributed import Client, LocalCluster
            # # cluster = LocalCluster()
            # # client = Client(cluster) # client = Client('127.0.0.1:8786')
            # start_time = timeit.default_timer()
            # def _tasks(nrep, nsgc, info, cell):
            #     res = []
            #     for nr in range(nrep):
            #         info["seed"] = seed + 3 * nsgc * nr
            #         post_cell = buildCell_dask(cell)
            #         psx = makeSynapses_dask(post_cell, info )
            #         # nvars = setupNeuron_dask(post_cell, info)
            #         repres = runNeuron_dask(post_cell, info, psx) # , nvars['rtime'], nvars['Vm'])
            #         res.append(repres)
            #     return(res)
            #
            # res = dask.compute(_tasks(self.nrep, self.n_sgc, info, self.cell)) #, scheduler='single-threaded')
            # elapsed = timeit.default_timer() - start_time
            # print(f"Parallel Elapsed time for {self.nrep:d} stimuli: {elapsed:f} secs")
            res = results
            for nr in range(self.nrep):
                self.pre_cells[nr] = res[nr]["pre_cells"]
                self.time[nr] = res[nr]["time"]
                self.xmtr = {k: v.to_python() for k, v in list(res[nr]["xmtr"].items())}
                self.vms[nr] = res[nr]["vm"]
                self.synapses[nr] = res[nr]["synapses"]
                self.xmtrs[nr] = self.xmtr


    def show(self):
        self.win = pg.GraphicsWindow()
        self.win.setBackground("w")
        Fs = self.Fs
        p1 = self.win.addPlot(
            title="Stimulus", row=0, col=0, labels={"bottom": "T (ms)", "left": "V"}
        )
        p1.plot(self.stim.time * 1000, self.stim.sound, pen=pg.mkPen("k", width=0.75))
        p1.setXLink(p1)

        p2 = self.win.addPlot(
            title="AN spikes",
            row=1,
            col=0,
            labels={"bottom": "T (ms)", "left": "AN spikes (first trial)"},
        )
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
            c.setData(
                xp.flatten(),
                yp.flatten(),
                connect="pairs",
                pen=pg.mkPen(pg.intColor(nr, self.nrep), hues=self.nrep, width=1.0),
            )
            p2.addItem(c)
        p2.setXLink(p1)

        p3 = self.win.addPlot(
            title="%s Spikes" % self.cell,
            row=2,
            col=0,
            labels={"bottom": "T (ms)", "left": "Trial #"},
        )
        xcn = []
        ycn = []
        xspks = []
        for k in range(self.nrep):
            bspk = PU.findspikes(self.time[k], self.vms[k], -35.0)
            xcn.extend(bspk)
            yr = k + np.zeros_like(bspk) + 0.2
            ycn.extend(yr)
        d = pg.PlotCurveItem()
        xp = np.repeat(np.array(xcn), 2)
        yp = np.repeat(np.array(ycn), 2)
        yp[1::2] = yp[::2] + 0.6
        d.setData(
            xp.flatten(), yp.flatten(), connect="pairs", pen=pg.mkPen("k", width=1.5)
        )
        p3.addItem(d)
        p3.setXLink(p1)

        p4 = self.win.addPlot(
            title="%s Vm" % self.cell,
            row=3,
            col=0,
            labels={"bottom": "T (ms)", "left": "Vm (mV)"},
        )
        for nr in range(self.nrep):
            p4.plot(
                self.time[nr],
                self.vms[nr],
                pen=pg.mkPen(pg.intColor(nr, self.nrep), hues=self.nrep, width=1.0),
            )
        p4.setXLink(p1)

        p5 = self.win.addPlot(
            title="xmtr", row=0, col=1, labels={"bottom": "T (ms)", "left": "gSyn"}
        )
        if synapseType == "multisite":
            for nr in [0]:
                syn = self.synapses[nr]
                j = 0
                for k in range(self.n_sgc):
                    synapse = syn[k]
                    for i in range(synapse.terminal.n_rzones):
                        p5.plot(
                            self.time[nr],
                            self.xmtrs[nr]["xmtr%04d" % j],
                            pen=pg.mkPen(
                                pg.intColor(nr, self.nrep), hues=self.nrep, width=1.0
                            ),
                        )
                        j = j + 1
        p5.setXLink(p1)

        p6 = self.win.addPlot(
            title="AN PSTH",
            row=1,
            col=1,
            labels={"bottom": "T (ms)", "left": "Sp/ms/trial"},
        )
        bins = np.arange(0, 200, 1)
        (hist, binedges) = np.histogram(xan, bins)
        curve6 = p6.plot(
            binedges,
            hist,
            stepMode=True,
            fillBrush=(0, 0, 0, 255),
            brush=pg.mkBrush("k"),
            fillLevel=0,
        )

        p7 = self.win.addPlot(
            title="%s PSTH" % self.cell,
            row=2,
            col=1,
            labels={"bottom": "T (ms)", "left": "Sp/ms/trial"},
        )
        bins = np.arange(0, 200, 1)
        (hist, binedges) = np.histogram(xcn, bins)
        curve7 = p7.plot(
            binedges,
            hist,
            stepMode=True,
            fillBrush=(0, 0, 0, 255),
            brush=pg.mkBrush("k"),
            fillLevel=0,
        )

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute AN only PSTH in postsynaptic cell"
    )
    parser.add_argument(
        type=str,
        dest="cell",
        default="bushy",
        choices=[
            "bushy",
            "tstellate",
            "dstellate",
            "octopus",
            "tuberculoventral",
            "pyramida",
        ],
        help="Select target cell",
    )
    parser.add_argument(
        "-s",
        "--stimulus",
        type=str,
        dest="stimulus",
        default="tone",
        choices=["tone", "SAM", "clicks"],
        help="Select stimulus from ['tone', 'SAM', 'clicks']",
    )
    parser.add_argument(
        "-t",
        "--type",
        type=str,
        dest="syntype",
        default="simple",
        choices=["simple", "multisite"],
        help="Set synapse type (simple, multisite)",
    )
    parser.add_argument(
        "-n",
        "--nrep",
        type=int,
        dest="nrep",
        default=10,
        help="Set number of repetitions",
    )

    parser.add_argument(
        "-d",
        "--dB",
        type=float,
        dest="dbspl",
        default=50.,
        help="Set sound intensity, dbSPL",
    )
    parser.add_argument(
        "--sr",
        type=int,
        dest="srgroup",
        default=2,
        help="SR group (1-high, 2-medium, 3-low)",
    )

    args = parser.parse_args()

    cell = args.cell
    synapseType = args.syntype

    print("cell type: ", cell)
    print('nrep: ', args.nrep)
    prot = SGCInputTestPSTH()
    prot.set_cell(cell)
    prot.set_db(args.dbspl)
    prot.set_sr(args.srgroup)
    
    prot.run(stimulus=args.stimulus, reps=args.nrep, simulator='cochlea', parallelize=True)
    prot.show()

    import sys

    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()
