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
from neuron import h
import pyqtgraph as pg
import pyqtgraph.multiprocess as mp
from pyqtgraph.Qt import QtGui, QtCore
from cnmodel import populations
from cnmodel.util import sound, random
from cnmodel.protocols import Protocol


class CNSoundStim(Protocol):
    def __init__(self, seed, temp=34.0, dt=0.025):
        Protocol.__init__(self)
        
        random.set_seed(seed)
        self.temp = temp
        self.dt = dt

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

        # Connect populations. 
        # This only defines the connections between populations; no synapses are 
        # created at this stage.
        self.sgc.connect(self.bushy, self.dstellate, self.tuberculoventral, self.tstellate)
        self.dstellate.connect(self.bushy, self.tstellate)  # should connect to dstellate as well?
        self.tuberculoventral.connect(self.bushy, self.tstellate)
        self.tstellate.connect(self.bushy)

        # Select cells to record from.
        # At this time, we actually instantiate the selected cells.
        # Select 4 cells centered around 16kHz
        frequencies = [16e3]  #[6e3, 16e3, 32e3]
        cells_per_band = 4
        for f in frequencies:
            bushy_cell_ids = self.bushy.select(cells_per_band, cf=f, create=True)
            #tstel_cell_ids = self.tstellate.select(cells_per_band, cf=f, create=True)  

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

        self.sgc.set_seed(seed)

    def run(self, stim):
        self.reset()
        self.sgc.set_sound_stim(stim, parallel=True)
        
        # set up recording vectors
        for pop in self.bushy, self.dstellate, self.tstellate, self.tuberculoventral:
            for ind in pop.real_cells():
                cell = pop.get_cell(ind)
                self[cell] = cell.soma(0.5)._ref_v
        self['t'] = h._ref_t
            
        h.tstop = stim.duration * 1000
        h.celsius = self.temp
        h.dt = self.dt
        
        print "init.."
        self.custom_init()
        print "start.."
        last_update = time.time()
        while h.t < h.tstop:
            h.fadvance()
            now = time.time()
            if now - last_update > 1.0:
                print "%0.2f / %0.2f" % (h.t, h.tstop)
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
    def __init__(self, prot, results):
        pg.QtGui.QSplitter.__init__(self, QtCore.Qt.Horizontal)
        
        self.selected = None
        
        self.prot = prot
        
        self.ctrl = QtGui.QWidget()
        self.layout = pg.QtGui.QVBoxLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.ctrl.setLayout(self.layout)
        self.addWidget(self.ctrl)
        
        self.stim_combo = pg.QtGui.QComboBox()
        self.layout.addWidget(self.stim_combo)
        self.results = {}
        for stim, result in results:
            key = 'f0: %0.0f  dBspl: %0.0f' % (stim.key()['f0'], stim.key()['dbspl'])
            self.results[key] = (stim, result)
            self.stim_combo.addItem(key)
        self.stim_combo.currentIndexChanged.connect(self.load_stim)

        self.network_tree = NetworkTree(self.prot)
        self.layout.addWidget(self.network_tree)
        
        self.pw = pg.GraphicsLayoutWidget()
        self.addWidget(self.pw)
        
        self.matrix_plot = self.pw.addPlot()
        self.matrix_plot.setLogMode(x=True, y=False)
        self.level_line = self.matrix_plot.addLine(y=50, movable=True)
        self.freq_line = self.matrix_plot.addLine(x=2, movable=True)
        
        self.pw.nextRow()
        self.cell_plot = self.pw.addPlot(labels={'left': 'Vm'}, title='Bushy cell Vm')

        self.pw.nextRow()
        self.input_plot = self.pw.addPlot(labels={'left': 'input #', 'bottom': 'time'}, title="Input spike times")
        self.input_plot.setXLink(self.cell_plot)
        
    def load_stim(self):
        key = str(self.stim_combo.currentText())
        results = self.results[key]
        self.selected_results = results[1]
        self.selected_stim = results[0]
        self.cell_plot.clear()
        real = self.prot.bushy.real_cells()
        for i, ind in enumerate(real):
            p = self.cell_plot.plot(self.selected_results['t'], 
                                    self.selected_results[('bushy', ind)][0], 
                                    pen=(i, len(real)*1.5), name='bushy-%d' % ind)
            p.curve.setClickable(True)
            p.sigClicked.connect(self.cell_curve_clicked)
            p.cell_ind = ind
        
    def cell_curve_clicked(self, c):
        if self.selected is not None:
            pen = self.selected.curve.opts['pen']
            pen.setWidth(1)
            self.selected.setPen(pen)
            
        pen = c.curve.opts['pen']
        pen.setWidth(3)
        c.setPen(pen)
        self.selected = c

        self.show_cell(c.cell_ind)
        
    def show_cell(self, ind):
        """Show spike trains of inputs to selected cell.
        """
        self.input_plot.clear()
        rec = self.prot.bushy._cells[ind]
        i = 0
        plots = []
        # plot spike times for all presynaptic cells
        for j, pop in enumerate((self.prot.dstellate, self.prot.sgc)):
            if pop not in rec['connections']:
                continue
            pre_inds = rec['connections'][pop]
            for preind in pre_inds:
                spikes = self.selected_results[(pop.type, preind)][1]
                y = np.ones(len(spikes)) * i
                self.input_plot.plot(spikes, y, pen=None, symbol=['o', 't'][j], symbolBrush=(i, 30))
                i += 1
                    
        # update matrix image
        self.matrix_plot.clear()
        fvals = set()
        lvals = set()
        
        # first get lists of all frequencies and levels in the matrix
        for stim, vec in self.results.values():
            fvals.add(stim.key()['f0'])
            lvals.add(stim.key()['dbspl'])
        fvals = sorted(list(fvals))
        lvals = sorted(list(lvals))
        
        # next count the number of spikes for the selected cell at each point in the matrix
        matrix = np.zeros((len(fvals), len(lvals)))
        for stim, vec in self.results.values():
            spikes = vec[('bushy', ind)][1]
            i = fvals.index(stim.key()['f0'])
            j = lvals.index(stim.key()['dbspl'])
            matrix[i, j] = len(spikes)
            
        # plot and scale the matrix image 
        # note that the origin (lower left) of each image pixel indicates its actual test freq/level. 
        self.matrix_img = pg.ImageItem(matrix)
        self.matrix_plot.addItem(self.matrix_img)
        self.matrix_img.setPos(np.log10(min(fvals)), min(lvals))
        self.matrix_img.scale((np.log10(max(fvals)) - np.log10(min(fvals))) / (len(fvals) - 1), 
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
        for cpop, conns in all_conns.items():
            pop_grp = QtGui.QTreeWidgetItem([cpop.type, str(conns)])
            item.addChild(pop_grp)


class NetworkVisualizer(pg.PlotWidget):
    def __init__(self, populations):
        self.pops = populations
        pg.PlotWidget.__init__(self)
        self.setLogMode(x=True, y=False)
        
        self.cells = pg.ScatterPlotItem()
        self.cells.setZValue(10)
        self.addItem(self.cells)
        
        self.connections = pg.PlotCurveItem()
        self.addItem(self.connections)
        
        # first assign positions of all cells
        cells = []
        colors = {'sgc': 'g', 'bushy': 'b', 'dstellate': 'y', 'tuberculoventral': 'r', 'tstellate': 'm'}
        for y,pop in enumerate(self.pops.values()):
            symbol = {'sgc': '+', 'bushy': 'o', 'tuberculoventral': 's', 'dstellate': 't', 'tstellate': 'x'}[pop.type]
            pop.cell_positions = []
            for cell in pop._cells:
                pos = (np.log10(cell['cf']), y)
                real = cell['cell'] != 0
                brush = pg.mkBrush(colors[pop.type]) if real else pg.mkBrush(255, 255, 255, 30)
                cells.append({'x': pos[0], 'y': pos[1], 'symbol': 'o' if real else 'x', 'brush': brush, 'pen': None})
                pop.cell_positions.append(pos)
        
        self.cells.setData(cells)
        
        self.getAxis('left').setTicks([list(enumerate(self.pops.keys()))])
        
        # now assign connection lines
        con_x = []
        con_y = []
        for pop in self.pops.values():
            for i,cell in enumerate(pop._cells):
                conns = cell['connections']
                if conns == 0:
                    continue
                for prepop, precells in conns.items():
                    p1 = pop.cell_positions[i]
                    for j in precells:
                        p2 = prepop.cell_positions[j]
                        con_x.extend([p1[0], p2[0]])
                        con_y.extend([p1[1], p2[1]])
        self.connections.setData(x=con_x, y=con_y, connect='pairs', pen=(255, 255, 255, 100))
        



if __name__ == '__main__':
    import pickle, os, sys
    app = pg.mkQApp()
    pg.dbg()
    
    # Create a sound stimulus and use it to generate spike trains for the SGC
    # population
    stims = []
    fmin = 4e3
    fmax = 32e3
    octavespacing = 0.3
    n_frequencies = int((1./octavespacing)*np.log2(fmax/1000.)/np.log2(fmin/1000.) - 1)
    #n_frequencies = 10
    n_levels = 5
    #fvals = fmin * (fmax/fmin)**(np.arange(n_frequencies) / (n_frequencies-1.))
    fvals = np.logspace(np.log2(fmin/1000.), np.log2(fmax/1000.), num=n_frequencies, endpoint=True, base=2)*1000.
    levels = np.linspace(20, 100, n_levels)
    print 'n frequencies: ', n_frequencies
    print("Frequencies:", fvals/1000.)
    print("Levels:", levels)

    path = os.path.dirname(__file__)
    cachepath = os.path.join(path, 'cache')
    if not os.path.isdir(cachepath):
        os.mkdir(cachepath)

    seed = 34657845
    prot = CNSoundStim(seed=seed)
    
    nv = NetworkVisualizer(prot.populations)
    nv.show()
    raise Exception()

    i = 0
    results = []
    for f in fvals:
        for db in levels:
            stim = sound.TonePip(rate=100e3, duration=0.1, f0=f, dbspl=db,
                                 ramp_duration=2.5e-3, pip_duration=0.04, 
                                 pip_start=[0.02])
        
            print("=== Start run %d/%d ===" % (i+1, len(fvals)*len(levels)))
            cachefile = os.path.join(cachepath, 'seed=%d_f0=%f_dbspl=%f.pk' % (seed, f, db))
            if '--ignore-cache' in sys.argv or not os.path.isfile(cachefile):
                result = prot.run(stim)
                pickle.dump(result, open(cachefile, 'wb'))
            else:
                print("  (Loading cached results)")
                result = pickle.load(open(cachefile, 'rb'))
            results.append((stim, result))
            i += 1

    nd = NetworkSimDisplay(prot, results)
    nd.show()
    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()
