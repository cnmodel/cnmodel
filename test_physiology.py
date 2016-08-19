"""
Test physiological response properties.
"""

import os, time
import numpy as np
from neuron import h
import pyqtgraph as pg
import pyqtgraph.multiprocess as mp
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
        #self.tstellate = populations.TStellate()

        # Connect populations. 
        # This only defines the connections between populations; no synapses are 
        # created at this stage.
        self.sgc.connect(self.bushy, self.dstellate)#, self.tstellate)
        #self.dstellate.connect(self.bushy)#, self.tstellate)
        #self.tstellate.connect(self.bushy)

        # Select cells to record from.
        # At this time, we actually instantiate the selected cells.
        # select 10 bushy cells closest to 16kHz
        bushy_cell_ids = self.bushy.select(4, cf=16e3, create=True)  
        # select 10 stellate cells closest to 16kHz
        #tstel_cell_ids = self.tstellate.select(10, cf=16e3, create=True)  

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
        self.sgc.set_sound_stim(stim)
        
        # set up recording vectors
        for pop in self.bushy, self.dstellate:
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


class NetworkSimDisplay(pg.QtGui.QWidget):
    def __init__(self, prot, results):
        pg.QtGui.QWidget.__init__(self)
        
        self.selected = None
        
        self.prot = prot
        
        self.layout = pg.QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.stim_combo = pg.QtGui.QComboBox()
        self.layout.addWidget(self.stim_combo, 0, 0)
        self.results = {}
        for stim, result in results:
            self.results[str(stim.key())] = (stim, result)
            self.stim_combo.addItem(str(stim.key()))
        self.stim_combo.currentIndexChanged.connect(self.load_stim)
        
        self.pw = pg.GraphicsLayoutWidget()
        self.layout.addWidget(self.pw, 1, 0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        
        self.matrix_plot = self.pw.addPlot()
        self.matrix_plot.setLogMode(x=True, y=False)
        self.level_line = self.matrix_plot.addLine(y=50, movable=True)
        self.freq_line = self.matrix_plot.addLine(x=2, movable=True)
        
        self.pw.nextRow()
        self.cell_plot = self.pw.addPlot()

        self.pw.nextRow()
        self.input_plot = self.pw.addPlot()
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
        
        
        
if __name__ == '__main__':
    import pickle, os, sys
    app = pg.mkQApp()
    
    # Create a sound stimulus and use it to generate spike trains for the SGC
    # population
    stims = []
    fmin = 4e3
    fmax = 40e3
    fn = 10
    fvals = fmin * (fmax/fmin)**(np.arange(fn) / (fn-1.))
    levels = np.linspace(0, 100, 11)
    #levels = map(int, sys.argv[1:])
    print("Frequencies:", fvals/1000.)
    print("Levels:", levels)
    
    path = os.path.dirname(__file__)
    cachepath = os.path.join(path, 'cache')
    if not os.path.isdir(cachepath):
        os.mkdir(cachepath)
    
    seed = 34657845
    prot = CNSoundStim(seed=seed)
    i = 0
    results = []
    for f in fvals:
        for db in levels:
            stim = sound.TonePip(rate=100e3, duration=0.1, f0=f, dbspl=db,
                                 ramp_duration=2.5e-3, pip_duration=0.04, 
                                 pip_start=[0.02])
        
            print("=== Start run %d/%d ===" % (i+1, len(fvals)*len(levels)))
            cachefile = os.path.join(cachepath, 'seed=%d_f0=%f_dbspl=%f.pk' % (seed, f, db))
            if not os.path.isfile(cachefile):
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
