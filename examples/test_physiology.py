"""
Test principal cell responses to tone pips of varying frequency and intensity.

This is an example of model construction from a very high level--we specify
only which populations of cells are present and which ones should be connected.
The population and cell classes take care of all the details of generating the
network needed to support a small number of output cells.

Note: run time for this example can be very long. To speed things up, reduce
n_frequencies or n_levels, or reduce the number of selected output cells (see
cells_per_band).
 
"""

import os, sys, time
from collections import OrderedDict
import numpy as np
import scipy.stats
from neuron import h
import pyqtgraph as pg
import pyqtgraph.multiprocess as mp
from pyqtgraph.Qt import QtGui, QtCore
from cnmodel import populations
from cnmodel.util import sound, random_seed
from cnmodel.protocols import Protocol
import timeit


class CNSoundStim(Protocol):
    def __init__(self, seed, temp=34.0, dt=0.025, synapsetype='simple'):
        Protocol.__init__(self)
        
        self.seed = seed
        self.temp = temp
        self.dt = dt
#        self.synapsetype = synapsetype  # simple or multisite
        
        # Seed now to ensure network generation is stable
        random_seed.set_seed(seed)
        # Create cell populations.
        # This creates a complete set of _virtual_ cells for each population. No 
        # cells are instantiated at this point.
        self.sgc = populations.SGC(model='dummy')
        self.bushy = populations.Bushy()
        self.dstellate = populations.DStellate()
        self.tstellate = populations.TStellate()
        self.tuberculoventral = populations.Tuberculoventral()
        
        pops = [self.sgc, self.dstellate, self.tuberculoventral, self.tstellate, self.bushy]
        self.populations = OrderedDict([(pop.type,pop) for pop in pops])
        
        # set synapse type to use in the sgc population - simple is fast, multisite is slower
        # (eventually, we could do this for all synapse types..)
        self.sgc._synapsetype = synapsetype
            
        # Connect populations. 
        # This only defines the connections between populations; no synapses are 
        # created at this stage.
        self.sgc.connect(self.bushy, self.dstellate, self.tuberculoventral, self.tstellate)
        self.dstellate.connect(self.bushy, self.tstellate)  # should connect to dstellate as well?
        self.tuberculoventral.connect(self.bushy, self.tstellate)
        self.tstellate.connect(self.bushy)

        # Select cells to record from.
        # At this time, we actually instantiate the selected cells.
        
        # Pick a single bushy cell near 16kHz, with medium-SR inputs
        bc = self.bushy.cells
        msr_cells = bc[bc['sgc_sr'] == 1]  # filter for msr cells
        ind = np.argmin(np.abs(msr_cells['cf'] - 16e3))  # find the one closest to 16kHz
        cell_id = msr_cells[ind]['id']
        self.bushy.create_cells([cell_id])  # instantiate just one cell

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

    def run(self, stim, seed):
        """Run the network simulation with *stim* as the sound source and a unique
        *seed* used to configure the random number generators.
        """
        self.reset()
        
        # Generate 2 new seeds for the SGC spike generator and for the NEURON simulation
        rs = np.random.RandomState()
        rs.seed(self.seed ^ seed)
        seed1, seed2 = rs.randint(0, 2**32, 2)
        random_seed.set_seed(seed1)
        self.sgc.set_seed(seed2)
        
        self.sgc.set_sound_stim(stim, parallel=False)
        
        # set up recording vectors
        for pop in self.bushy, self.dstellate, self.tstellate, self.tuberculoventral:
            for ind in pop.real_cells():
                cell = pop.get_cell(ind)
                self[cell] = cell.soma(0.5)._ref_v
        self['t'] = h._ref_t
            
        h.tstop = stim.duration * 1000
        h.celsius = self.temp
        h.dt = self.dt
        
        print("init..")
        self.custom_init()
        print("start..")
        last_update = time.time()
        while h.t < h.tstop:
            h.fadvance()
            now = time.time()
            if now - last_update > 1.0:
                print("%0.2f / %0.2f" % (h.t, h.tstop))
                last_update = now
        
        # record vsoma and spike times for all cells
        vec = {}
        for k in self._vectors:
            v = self[k].copy()
            if k == 't':
                vec[k] = v
                continue
            spike_inds = np.argwhere((v[1:]>-20) & (v[:-1]<=-20))[:,0]
            spikes = self['t'][spike_inds]
            pop = k.type
            cell_ind = getattr(self, pop).get_cell_index(k)
            vec[(pop, cell_ind)] = [v, spikes]
        
        # record SGC spike trains
        for ind in self.sgc.real_cells():
            cell = self.sgc.get_cell(ind)
            vec[('sgc', ind)] = [None, cell._spiketrain]
        
        return vec


class NetworkSimDisplay(pg.QtGui.QSplitter):
    def __init__(self, prot, results, baseline, response):
        pg.QtGui.QSplitter.__init__(self, QtCore.Qt.Horizontal)
        self.selected_cell = None
        
        self.prot = prot
        self.baseline = baseline  # (start, stop)
        self.response = response  # (start, stop)

        self.ctrl = QtGui.QWidget()
        self.layout = pg.QtGui.QVBoxLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.ctrl.setLayout(self.layout)
        self.addWidget(self.ctrl)

        self.nv = NetworkVisualizer(prot.populations)
        self.layout.addWidget(self.nv)
        self.nv.cell_selected.connect(self.nv_cell_selected)
                
        self.stim_combo = pg.QtGui.QComboBox()
        self.layout.addWidget(self.stim_combo)
        self.trial_combo = pg.QtGui.QComboBox()
        self.layout.addWidget(self.trial_combo)
        self.results = OrderedDict()
        self.stim_order = []
        freqs = set()
        levels = set()
        max_iter = 0
        for k,v in list(results.items()):
            f0, dbspl, iteration = k
            max_iter = max(max_iter, iteration)
            stim, result = v
            key = 'f0: %0.0f  dBspl: %0.0f' % (f0, dbspl)
            self.results.setdefault(key, [stim, {}])
            self.results[key][1][iteration] = result
            self.stim_order.append((f0, dbspl))
            freqs.add(f0)
            levels.add(dbspl)
            self.stim_combo.addItem(key)
        self.freqs = sorted(list(freqs))
        self.levels = sorted(list(levels))
        self.iterations = max_iter + 1
        self.trial_combo.addItem("all trials")
        for i in range(self.iterations):
            self.trial_combo.addItem(str(i))
            
        self.stim_combo.currentIndexChanged.connect(self.stim_selected)
        self.trial_combo.currentIndexChanged.connect(self.trial_selected)

        self.tuning_plot = pg.PlotWidget()
        self.tuning_plot.setLogMode(x=True, y=False)
        self.tuning_plot.scene().sigMouseClicked.connect(self.tuning_plot_clicked)
        self.layout.addWidget(self.tuning_plot)

        self.tuning_img = pg.ImageItem()
        self.tuning_plot.addItem(self.tuning_img)
        
        df = np.log10(self.freqs[1]) - np.log10(self.freqs[0])
        dl = self.levels[1] - self.levels[0]
        self.stim_rect = QtGui.QGraphicsRectItem(QtCore.QRectF(0, 0, df, dl))
        self.stim_rect.setPen(pg.mkPen('c'))
        self.stim_rect.setZValue(20)
        self.tuning_plot.addItem(self.stim_rect)

        #self.network_tree = NetworkTree(self.prot)
        #self.layout.addWidget(self.network_tree)
        
        self.pw = pg.GraphicsLayoutWidget()
        self.addWidget(self.pw)
        
        self.stim_plot = self.pw.addPlot()
        self.pw.ci.layout.setRowFixedHeight(0, 100)
        
        self.pw.nextRow()
        self.cell_plot = self.pw.addPlot(labels={'left': 'Vm'})

        self.pw.nextRow()
        self.input_plot = self.pw.addPlot(labels={'left': 'input #', 'bottom': 'time'}, title="Input spike times")
        self.input_plot.setXLink(self.cell_plot)
        self.stim_plot.setXLink(self.cell_plot)
        
        self.stim_selected()
        
    def update_stim_plot(self):
        stim = self.selected_stim
        self.stim_plot.plot(stim.time*1000, stim.sound, clear=True, antialias=True)
        
    def update_raster_plot(self):
        self.input_plot.clear()
        if self.selected_cell is None:
            return
        pop, ind = self.selected_cell
        
        rec = pop._cells[ind]
        i = 0
        plots = []
        # plot spike times for all presynaptic cells
        labels = []
        if rec['connections'] == 0:
            return
        
        pop_colors = {'dstellate': 'y', 'tuberculoventral': 'r', 'sgc': 'g', 'tstellate': 'b'}
        pop_symbols = {'dstellate': 'x', 'tuberculoventral': '+', 'sgc': 't', 'tstellate': 'o'}
        pop_order = [self.prot.sgc, self.prot.dstellate, self.prot.tuberculoventral]
        trials = self.selected_trials()
        for pop in pop_order:
            pre_inds = rec['connections'].get(pop, [])
            for preind in pre_inds:
                # iterate over all trials
                for j in trials:
                    result = self.selected_results[j]
                    spikes = result[(pop.type, preind)][1]
                    y = np.ones(len(spikes)) * i + j / (2. * len(self.selected_results))
                    self.input_plot.plot(spikes, y, pen=None, symbolBrush=pop_colors[pop.type], symbol='+', symbolPen=None)
                i += 1
                labels.append(pop.type + " " + str(preind))
        self.input_plot.getAxis('left').setTicks([list(enumerate(labels))])

    def update_cell_plot(self):
        self.cell_plot.clear()
        if self.selected_cell is None:
            return
        pop, cell_ind = self.selected_cell
        
        self.cell_plot.setTitle("%s %d   %s" % (pop.type, cell_ind, str(self.stim_combo.currentText())))
        trials = self.selected_trials()
        for i in trials:
            result = self.selected_results[i]
            y = result[(pop.type, cell_ind)][0]
            if y is not None:
                p = self.cell_plot.plot(self.selected_results[0]['t'], y, 
                    name='%s-%d' % self.selected_cell, antialias=True, 
                    pen=(i, len(self.selected_results)*1.5))
        #p.curve.setClickable(True)
        #p.sigClicked.connect(self.cell_curve_clicked)
        #p.cell_ind = ind

    def tuning_plot_clicked(self, event):
        spos = event.scenePos()
        stimpos = self.tuning_plot.plotItem.vb.mapSceneToView(spos)
        x = 10**stimpos.x()
        y = stimpos.y()
        
        best = None
        for stim, result in list(self.results.values()):
            f0 = stim.opts['f0']
            dbspl = stim.opts['dbspl']
            if x < f0 or y < dbspl:
                continue
            if best is None:
                best = stim
                continue
            if f0 > best.opts['f0'] or dbspl > best.opts['dbspl']:
                best = stim
                continue
        
        if best is None:
            return
        self.select_stim(best.opts['f0'], best.opts['dbspl'])

    def nv_cell_selected(self, nv, cell):
        self.select_cell(*cell)

    def stim_selected(self):
        key = str(self.stim_combo.currentText())
        results = self.results[key]
        self.selected_results = results[1]
        self.selected_stim = results[0]
        self.update_stim_plot()
        self.update_raster_plot()
        self.update_cell_plot()
        
        self.stim_rect.setPos(np.log10(results[0].opts['f0']), results[0].opts['dbspl'])

    def trial_selected(self):
        self.update_raster_plot()
        self.update_cell_plot()
        self.update_tuning()

    def selected_trials(self):
        if self.trial_combo.currentIndex() == 0:
            return list(range(self.iterations))
        else:
            return [self.trial_combo.currentIndex() - 1]

    def select_stim(self, f0, dbspl):
        i = self.stim_order.index((f0, dbspl))
        self.stim_combo.setCurrentIndex(i)
        
    def select_cell(self, pop, cell_id):
        self.selected_cell = pop, cell_id
        self.update_tuning()
        self.update_cell_plot()
        self.update_raster_plot()
        
    #def cell_curve_clicked(self, c):
        #if self.selected is not None:
            #pen = self.selected.curve.opts['pen']
            #pen.setWidth(1)
            #self.selected.setPen(pen)
            
        #pen = c.curve.opts['pen']
        #pen.setWidth(3)
        #c.setPen(pen)
        #self.selected = c

        #self.show_cell(c.cell_ind)

    def update_tuning(self):
        # update matrix image
        if self.selected_cell is None:
            return
        
        pop, ind = self.selected_cell
        fvals = set()
        lvals = set()
        
        # first get lists of all frequencies and levels in the matrix
        for stim, vec in list(self.results.values()):
            fvals.add(stim.key()['f0'])
            lvals.add(stim.key()['dbspl'])
        fvals = sorted(list(fvals))
        lvals = sorted(list(lvals))

        # Get spontaneous rate statistics
        spont_spikes = 0
        spont_time = 0
        for stim, iterations in list(self.results.values()):
            for vec in list(iterations.values()):
                spikes = vec[(pop.type, ind)][1]
                spont_spikes += ((spikes >= self.baseline[0]) & (spikes < self.baseline[1])).sum()
                spont_time += self.baseline[1] - self.baseline[0]
        spont_rate = spont_spikes / spont_time
        
        # next count the number of spikes for the selected cell at each point in the matrix
        matrix = np.zeros((len(fvals), len(lvals)))
        trials = self.selected_trials()
        for stim, iteration in list(self.results.values()):
            for i in trials:
                vec = iteration[i]
                spikes = vec[(pop.type, ind)][1]
                n_spikes = ((spikes >= self.response[0]) & (spikes < self.response[1])).sum()
                i = fvals.index(stim.key()['f0'])
                j = lvals.index(stim.key()['dbspl'])
                matrix[i, j] += n_spikes - spont_rate * (self.response[1]-self.response[0])
        matrix /= self.iterations
        
        # plot and scale the matrix image 
        # note that the origin (lower left) of each image pixel indicates its actual test freq/level. 
        self.tuning_img.setImage(matrix)
        self.tuning_img.resetTransform()
        self.tuning_img.setPos(np.log10(min(fvals)), min(lvals))
        self.tuning_img.scale((np.log10(max(fvals)) - np.log10(min(fvals))) / (len(fvals) - 1), 
                              (max(lvals) - min(lvals)) / (len(lvals) - 1))
        


class NetworkTree(QtGui.QTreeWidget):
    def __init__(self, prot):
        self.prot = prot
        QtGui.QTreeWidget.__init__(self)
        self.setColumnCount(2)
        
        self.update_tree()
        
    def update_tree(self):
        for pop_name in ['bushy', 'tstellate', 'dstellate', 'tuberculoventral', 'sgc']:
            if not hasattr(self.prot, pop_name):
                continue
            pop = getattr(self.prot, pop_name)
            grp = QtGui.QTreeWidgetItem([pop_name])
            self.addTopLevelItem(grp)
            for cell in pop.real_cells():
                self.add_cell(grp, pop, cell)
                
    def add_cell(self, grp_item, pop, cell):
        item = QtGui.QTreeWidgetItem([str(cell)])
        grp_item.addChild(item)
        all_conns = pop.cell_connections(cell)
        if all_conns == 0:
            return
        for cpop, conns in list(all_conns.items()):
            pop_grp = QtGui.QTreeWidgetItem([cpop.type, str(conns)])
            item.addChild(pop_grp)


class NetworkVisualizer(pg.PlotWidget):
    
    cell_selected = pg.QtCore.Signal(object, object)
    
    def __init__(self, populations):
        self.pops = populations
        pg.PlotWidget.__init__(self)
        self.setLogMode(x=True, y=False)
        
        self.cells = pg.ScatterPlotItem(clickable=True)
        self.cells.setZValue(10)
        self.addItem(self.cells)
        self.cells.sigClicked.connect(self.cells_clicked)
        
        self.selected = pg.ScatterPlotItem()
        self.selected.setZValue(20)
        self.addItem(self.selected)
        
        self.connections = pg.PlotCurveItem()
        self.addItem(self.connections)
        
        # first assign positions of all cells
        cells = []
        for y,pop in enumerate(self.pops.values()):
            pop.cell_spots = []
            pop.fwd_connections = {}
            for i,cell in enumerate(pop._cells):
                pos = (np.log10(cell['cf']), y)
                real = cell['cell'] != 0
                if not real:
                    pop.cell_spots.append(None)
                    continue
                brush = pg.mkBrush('b') if real else pg.mkBrush(255, 255, 255, 30)
                spot = {'x': pos[0], 'y': pos[1], 'symbol': 'o' if real else 'x', 'brush': brush, 'pen': None, 'data': (pop, i)}
                cells.append(spot)
                pop.cell_spots.append(spot)
        
        self.cells.setData(cells)
        
        self.getAxis('left').setTicks([list(enumerate(self.pops.keys()))])
        
        # now assign connection lines and record forward connectivity
        con_x = []
        con_y = []
        for pop in list(self.pops.values()):
            for i,cell in enumerate(pop._cells):
                conns = cell['connections']
                if conns == 0:
                    continue
                for prepop, precells in list(conns.items()):
                    spot = pop.cell_spots[i]
                    if spot is None:
                        continue
                    p1 = spot['x'], spot['y']
                    for j in precells:
                        prepop.fwd_connections.setdefault(j, [])
                        prepop.fwd_connections[j].append((pop, i))
                        spot2 = prepop.cell_spots[j]
                        if spot2 is None:
                            return
                        p2 = spot2['x'], spot2['y']
                        con_x.extend([p1[0], p2[0]])
                        con_y.extend([p1[1], p2[1]])
        self.connections.setData(x=con_x, y=con_y, connect='pairs', pen=(255, 255, 255, 60))
        
    def cells_clicked(self, *args):
        selected = None
        for spot in args[1]:
            # find the first real cell
            pop, i = spot.data()
            if pop._cells[i]['cell'] != 0:
                selected = spot
                break
        if selected is None:
            self.selected.hide()
            return
        
        rec = pop._cells[i]
        pos = selected.pos()
        spots = [{'x': pos.x(), 'y': pos.y(), 'size': 15, 'symbol': 'o', 'pen': 'y', 'brush': 'b'}]

        # display presynaptic cells
        if rec['connections'] != 0:
            for prepop, preinds in list(rec['connections'].items()):
                for preind in preinds:
                    spot = prepop.cell_spots[preind].copy()
                    spot['size'] = 15
                    spot['brush'] = 'r'
                    spots.append(spot)
                
        # display postsynaptic cells
        for postpop, postind in pop.fwd_connections.get(i, []):
            spot = postpop.cell_spots[postind].copy()
            spot['size'] = 15
            spot['brush'] = 'g'
            spots.append(spot)
        
        self.selected.setData(spots)
        self.selected.show()
        
        self.cell_selected.emit(self, selected.data())


if __name__ == '__main__':
    import pickle, os, sys
    app = pg.mkQApp()
    pg.dbg()
    
    # Create a sound stimulus and use it to generate spike trains for the SGC
    # population
    stims = []
    parallel = True
    
    nreps = 5
    fmin = 4e3
    fmax = 32e3
    octavespacing = 1/8.
    #octavespacing = 1.
    n_frequencies = int(np.log2(fmax/fmin) / octavespacing) + 1
    fvals = np.logspace(np.log2(fmin/1000.), np.log2(fmax/1000.), num=n_frequencies, endpoint=True, base=2)*1000.
    
    n_levels = 11
    #n_levels = 3
    levels = np.linspace(20, 100, n_levels)
    
    print(("Frequencies:", fvals/1000.))
    print(("Levels:", levels))

    syntype = 'multisite'
    path = os.path.dirname(__file__)
    cachepath = os.path.join(path, 'cache')
    if not os.path.isdir(cachepath):
        os.mkdir(cachepath)

    seed = 34657845
    prot = CNSoundStim(seed=seed, synapsetype=syntype)
    i = 0
    
    start_time = timeit.default_timer()
    
    #stimpar = {'dur': 0.06, 'pip': 0.025, 'start': [0.02], 'baseline': [10, 20], 'response': [20, 45]}
    stimpar = {'dur': 0.2, 'pip': 0.04, 'start': [0.1], 'baseline': [50, 100], 'response': [100, 140]}
    tasks = []
    for f in fvals:
        for db in levels:
            for i in range(nreps):
                tasks.append((f, db, i))
            
    results = {}
    workers = 1 if not parallel else None
    tot_runs = len(fvals) * len(levels) * nreps
    with mp.Parallelize(enumerate(tasks), results=results, progressDialog='Running parallel simulation..', workers=workers) as tasker:
        for i, task in tasker:
            f, db, iteration = task
            stim = sound.TonePip(rate=100e3, duration=stimpar['dur'], f0=f, dbspl=db,  # dura 0.2, pip_start 0.1 pipdur 0.04
                                    ramp_duration=2.5e-3, pip_duration=stimpar['pip'], 
                                    pip_start=stimpar['start'])
    
            print(("=== Start run %d/%d ===" % (i+1, tot_runs)))
            cachefile = os.path.join(cachepath, 'seed=%d_f0=%f_dbspl=%f_syntype=%s_iter=%d.pk' % (seed, f, db, syntype, iteration))
            if '--ignore-cache' in sys.argv or not os.path.isfile(cachefile):
                result = prot.run(stim, seed=i)
                pickle.dump(result, open(cachefile, 'wb'))
            else:
                print("  (Loading cached results)")
                result = pickle.load(open(cachefile, 'rb'))
            tasker.results[(f, db, iteration)] = (stim, result)
            print(('--- finished run %d/%d ---' % (i+1, tot_runs)))
        
    # get time of run before display
    elapsed = timeit.default_timer() - start_time
    print('Elapsed time for %d stimuli: %f  (%f sec per stim), synapses: %s' % (len(tasks), elapsed, elapsed/len(tasks), prot.bushy._synapsetype))
    
    nd = NetworkSimDisplay(prot, results, baseline=stimpar['baseline'], response=stimpar['response'])
    nd.show()
    
    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()
