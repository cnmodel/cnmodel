#!/usr/bin/python
"""
Basic test of initialization of multiple cells in the model, and running multiple cells at one time.
Plots the resposnes to a series of current injections for most implemented baseic cell types in
in cnmodel.

Usage:
    python examples/toy_model.py  (no arguments)

"""

from __future__ import print_function


import sys
from neuron import h
import numpy as np
import cnmodel.cells as cells
from cnmodel.protocols import Protocol
from cnmodel.util import custom_init
from collections import OrderedDict
import re
import pyqtgraph.exporters
from cnmodel.util import pyqtgraphPlotHelpers as PH
from cnmodel.protocols import IVCurve


try:  # check for pyqtgraph install
    import pyqtgraph as pg
except ImportError:
    raise ImportError("This model requires pyqtgraph")

from cnmodel.util.stim import make_pulse

def autorowcol(n):
    """
    return a reasonable layout (cols, rows) for n plots on a page
    up to 16.
    Otherwise return floor(sqrt(n)) + 1 for both r and c.
    """
    nmap = {1: (1,1), 2: (2,1), 3: (3,1), 4: (2,2), 5: (3,2), 6: (3,2), 7: (3,3),
            8: (3,3), 9: (3,3), 10: (3,4), 11: (3,4), 12: (3,4), 13: (4,4), 14: (4,4),
            15: (4,4), 16: (4,4)}
    if n <= 16:
        return nmap[n][0], nmap[n][1]
    else:
        nx = np.floor(np.sqrt(n)) + 1
        return nx, nx


def makeLayout(cols=1, rows=1, letters=True, margins=4, spacing=4, nmax=None):
    """
    Create a multipanel plot, returning the various pyptgraph elements.
    The layout is always a rectangular grid with shape (cols, rows)
    if letters is true, then the plot is labeled "A, B, C..."
    margins sets the margins around the outside of the plot
    spacing sets the spacing between the elements of the grid
    """
    import string
    letters = string.ascii_uppercase
    widget = pg.QtGui.QWidget()
    gridLayout = pg.QtGui.QGridLayout()
    widget.setLayout(gridLayout)
    gridLayout.setContentsMargins(margins, margins, margins, margins)
    gridLayout.setSpacing(spacing)
    plots = [[0 for x in range(cols)] for x in range(rows)]
    i = 0
    for c in range(cols):
        for r in range(rows):
            plots[r][c] = pg.PlotWidget()
            gridLayout.addWidget(plots[r][c], r, c)
            #labelUp(plots[r][c], 'T(s)', 'Y', title = letters[i])
            i += 1
            if i > 25:
                i = 0
            if nmax is not None and i >= nmax:
                break  # that's all - leave out empty plots

    return(plots, widget, gridLayout)

def getnextrowcol(plx, row, col, cols):
    col += 1
    if col >= cols:
        col = 0
        row += 1
    return (plx[row][col], row, col)
    

class Toy(Protocol):
    """
    Calls to encapsulate the model runs
    Run a set of cells with defined parameters to show excitability patterns. 
    Note that cells from Rothman and Manis are run at 22C; others at various
    temperatures depending on how they were initially measured and defined.
    
    """
    def __init__(self):
        super(Toy, self).__init__()

    def current_name(self, name, n):
        """
        From the name of the current model, get the current injection information
        
        Parameters
        ---------
        name : str (no default)
            name of the cell type
        n : int (no default)
        """
        if len(self.celltypes[name][3]) > 2:
            injs = self.celltypes[name][3]
            injarr = np.linspace(injs[0], injs[1], injs[2], endpoint=True)
            return '%.3f' % injarr[n]
        else:
            return '%.3f' % self.celltypes[name][3][n]

    def getname(self, cell, ninj):
        name = self.make_name(cell)
        iname = self.current_name(name, ninj)
        nname = name + ' ' + iname
        return name, nname

    def make_name(self, cell):
        return cell + ', ' + self.celltypes[cell][1] + ':'
        
    def run(self):
        sre = re.compile('(?P<cell>\w+)(?:[, ]*)(?P<type>[\w-]*)(?:[, ]*)(?P<species>[\w-]*)')  # regex for keys in cell types
        self.celltypes = OrderedDict([('Bushy, II', (cells.Bushy,      "II",   'guineapig', (-0.5, 0.5, 11), 22)),
                                ('Bushy, II-I',     (cells.Bushy,      "II-I", 'guineapig', (-0.5,0.5, 11), 22)),
                                ('DStellate, I-II', (cells.DStellate,  'I-II', 'guineapig', (-0.3, 0.3, 9), 22)),
                                ('TStellate, I-c',  (cells.TStellate,  "I-c",  'guineapig', (-0.15, 0.15, 9), 22)), 
                                ('TStellate, I-t',  (cells.TStellate,  "I-t",  'guineapig', (-0.15, 0.15, 9), 22)),
                                ('Octopus, II-o',   (cells.Octopus,    'II-o', 'guineapig', (-2.5, 2.5, 11), 22)),
                                ('Bushy, II, Mouse',       (cells.Bushy,      "II",   'mouse', (-1, 1.2, 13), 34)),
                                ('TStellate, I-c, Mouse',  (cells.TStellate,  "I-c",  'mouse', (-1, 1, 9), 34)), 
                                ('DStellate, I-II, Mouse', (cells.DStellate,  'I-II', 'mouse', (-0.5, 0.5, 9), 34)),
                                ('Pyramidal, I, Rat',    (cells.Pyramidal,  'I',    'rat', (-0.3, 0.4, 11), 34)),
                                ('Cartwheel, I, Mouse',    (cells.Cartwheel,  'I',    'mouse', (-0.5, 0.5, 9), 34)),
                                ('Tuberculoventral, TVmouse, Mouse', (cells.Tuberculoventral, 'TVmouse', 'mouse', (-0.35, 1, 11), 34)),
                                ('SGC, sgc-bm, Mouse',         (cells.SGC,        "sgc-bm",   'mouse', (-0.2, 0.6, 9), 34)),
                                ('SGC, sgc-a, Mouse',          (cells.SGC,         "sgc-a",   'mouse', (-0.2, 0.6, 9), 34)),
                                ])

        dt = 0.025
        h.dt = dt
        h.celsius = 22

        stim = {
            'NP': 1,
            'delay': 10,
            'dur': 100,
            'amp': 0.,
            'dt': h.dt,
            }
        tend = stim['delay'] + stim['dur'] + 20.

        netcells = {}
        for c in list(self.celltypes.keys()):
            g = sre.match(c)
            cellname = g.group('cell')
            modelType = g.group('type')
            species = self.celltypes[c][2]
            if g.group('type') == '':
                netcells[c] = self.celltypes[c][0].create()
            else:
                netcells[c] = self.celltypes[c][0].create(modelType=modelType, species=species, debug=False)
        # dicts to hold data
        pl = OrderedDict([])
        pl2 = OrderedDict([])
        rvec = OrderedDict([])
        vec = OrderedDict([])
        istim = OrderedDict([])
        ncells = len(list(self.celltypes.keys()))
        #
        # build plotting area
        #
        app = pg.mkQApp()
        self.win = pg.GraphicsWindow()
        self.win.setBackground('w')
        self.win.resize(800, 600)
        cols, rows = autorowcol(ncells)

        row = 0
        col = 0
        labelStyle = {'color': '#000', 'font-size': '9pt', 'weight': 'normal'}
        tickStyle = pg.QtGui.QFont('Arial', 9, pg.QtGui.QFont.Light)
        self.iv = IVCurve()  # use standard IVCurve here...
        for n, name in enumerate(self.celltypes.keys()):
            nrn_cell = netcells[name]  # get the Neuron object we are using for this cell class
            injcmds = list(self.celltypes[name][3])  # list of injections
            injcmds[2] = (injcmds[1]-injcmds[0])/(float(injcmds[2]-1))  # convert to pulse format for IVCurve
            temperature = self.celltypes[name][4]
            nrn_cell.set_temperature(float(temperature))
            ninjs = len(injcmds)
            print( 'cell: ', name)
            # print( 'injs: ', injcmds)
            pl[name] = self.win.addPlot(labels={'left': 'V (mV)', 'bottom': 'Time (ms)'})
            PH.nice_plot(pl[name])
            pl[name].setTitle(title=name, font=pg.QtGui.QFont('Arial', 10) )
            col += 1
            if col >= cols:
                col = 0
                self.win.nextRow()
                row += 1
            self.iv.reset()
            self.iv.run({'pulse': [injcmds]}, nrn_cell, durs=(stim['delay'], stim['dur'], 20.),
                   sites=None, reppulse=None, temp=float(temperature))
            for k in range(len(self.iv.voltage_traces)):
                pl[name].plot(self.iv.time_values, self.iv.voltage_traces[k], pen=pg.mkPen('k', width=0.75))
            pl[name].setRange(xRange=(0., 130.), yRange=(-160., 40.))
            PH.noaxes(pl[name])
            PH.calbar(pl[list(self.celltypes.keys())[0]], calbar=[0, -120., 10., 20.], unitNames={'x': 'ms', 'y': 'mV'})

            text = (u"{0:2d}\u00b0C {1:.2f}-{2:.2f} nA".format(int(temperature), np.min(self.iv.current_cmd), np.max(self.iv.current_cmd)))
            ti = pg.TextItem(text, anchor=(1, 0))
            ti.setFont(pg.QtGui.QFont('Arial', 9))
            ti.setPos(120., -120.)
            pl[name].addItem(ti)
        # get overall Rin, etc; need to initialize all cells
        nrn_cell.cell_initialize()
        for n, name in enumerate(self.celltypes.keys()):
            nrn_cell = netcells[name]
            nrn_cell.vm0 = nrn_cell.soma.v
            pars = nrn_cell.compute_rmrintau(auto_initialize=False)
            print('{0:>14s} [{1:>24s}]   *** Rin = {2:6.1f} M\ohm  Tau = {3:6.1f} ms   Vm = {4:6.1f} mV'.
                format(nrn_cell.status['name'], name, pars['Rin'], pars['tau'], pars['v']))


if __name__ == "__main__":
    t=Toy()
    t.run()
    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()
    # exporter = pg.exporters.ImageExporter(t.win.scene())
    # exporter.export('~/Desktop/Model_Figure2.svg')
