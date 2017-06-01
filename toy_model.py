#!/usr/bin/python
from __future__ import print_function
__author__ = 'pbmanis'

"""
Basic test of initialization of multiple cells in the model, and running multiple cells at one time.

"""

import sys
from neuron import h
import numpy as np
import cnmodel.cells as cells
from cnmodel.protocols import Protocol
from collections import OrderedDict
import re

try:  # check for pyqtgraph install
    import pyqtgraph as pg
    HAVE_PG = True
except ImportError:
    HAVE_PG = False

if HAVE_PG:
    from PyQt4 import QtGui

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
    widget = QtGui.QWidget()
    gridLayout = QtGui.QGridLayout()
    widget.setLayout(gridLayout)
    gridLayout.setContentsMargins(margins, margins, margins, margins)
    gridLayout.setSpacing(spacing)
    plots = [[0 for x in xrange(cols)] for x in xrange(rows)]
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
        if len(self.celltypes[name][2]) > 2:
            injs = self.celltypes[name][2]
            injarr = np.linspace(injs[0], injs[1], injs[2], endpoint=True)
            return '%.3f' % injarr[n]
        else:
            return '%.3f' % self.celltypes[name][2][n]

    def getname(self, cell, ninj):
        name = self.make_name(cell)
        iname = self.current_name(name, ninj)
        nname = name + ' ' + iname
        return name, nname

    def make_name(self, cell):
        return cell + ', ' + self.celltypes[cell][1] + ':'
        
    def run(self):
        sre = re.compile('(?P<cell>\w+)(?:[, ]*)(?P<type>[\w-]*)')  # regex for keys in cell types
        self.celltypes = OrderedDict([('Bushy, II', (cells.Bushy, "II", (-0.5, 0.5, 11), 22)),
                                ('Bushy, II-I', (cells.Bushy, "II-I", (-0.5,0.5, 11), 22)),
                                ('Octopus, II-o', (cells.Octopus, 'II-o', (-2.5, 2.5, 11), 22)),
                                ('TStellate, I-c', (cells.TStellate, "I-c", (-0.15, 0.15, 9), 22)), 
                                ('TStellate, I-t', (cells.TStellate, "I-t", (-0.15, 0.15, 9), 22)),
                                ('DStellate, I-II', (cells.DStellate, 'I-II', (-0.2, 0.2, 9), 22)),
                                ('Pyramidal, I', (cells.Pyramidal, 'I', (-0.3, 0.4, 11), 34)),
                                ('Cartwheel, I', (cells.Cartwheel, 'I', (-0.2, 0.2, 9), 34)),
                                ('Tubercuoventral, I', (cells.Tuberculoventral, 'I', (-0.35, 1, 11), 34)),
                                ('SGC, bm', (cells.SGC, "bm", (-0.2, 0.2, 5), 22)),
                                ('SGC, a', (cells.SGC, "a", (-0.2, 0.2, 5), 22)),
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
        for c in self.celltypes.keys():
            g = sre.match(c)
            cellname = g.group('cell')
            modelType = g.group('type')
            if g.group('type') == '':
                netcells[c] = self.celltypes[c][0].create()
            else:
                netcells[c] = self.celltypes[c][0].create(modelType=modelType, debug=False)
        # dicts to hold data
        pl = OrderedDict([])
        pl2 = OrderedDict([])
        rvec = OrderedDict([])
        vec = OrderedDict([])
        istim = OrderedDict([])
        ncells = len(self.celltypes.keys())
        #
        # build plotting area
        #
        app = pg.mkQApp()
        win = pg.GraphicsWindow()
        self.win=win
        win.resize(800, 600)
        cols, rows = autorowcol(ncells)
        (plx, widget, gridlayout) = makeLayout(cols=cols, rows=rows, nmax=ncells)
        win.setLayout(gridlayout)
        win.show()
        row = 0
        col = 0
        
        for n, name in enumerate(self.celltypes.keys()):
            pl[name] = plx[row][col]
            pl[name].setLabels(left='Vm (mV)', bottom='Time (ms)')
            pl[name].setTitle(name)
            col += 1
            if col >= cols:
                col = 0
                row += 1
        # initialize some recording vectors
        for n, name in enumerate(self.celltypes.keys()):
            nrn_cell = netcells[name]  # get the Neuron object we are using for this cell class
            injcmds = self.celltypes[name][2]  # list of injections
            temperature = self.celltypes[name][3]
            ninjs = len(injcmds)
            if ninjs > 2:  # 2 values or a range?
                injcmds = np.linspace(injcmds[0], injcmds[1], num=injcmds[2], endpoint=True)
                ninjs = len(injcmds)
            print( 'cell: ', name)
            print( 'injs: ', injcmds)
            for ninj in range(ninjs):  # for each current level
                iname = self.current_name(name, ninj)
                runname = name + ' ' + iname 
                rvec[runname] = {'v_soma': h.Vector(), 'i_inj': h.Vector(), 'time': h.Vector()}
                stim['amp'] = injcmds[ninj]  # currents[name]
                (secmd, maxt, tstims) = make_pulse(stim)
                if ninj == 0: # install stimulus electronde first run only
                    istim[name] = h.iStim(0.5, sec=nrn_cell.soma)
                    istim[name].delay = 5.
                    istim[name].dur = 1e9 # these actually do not matter...
                vec[runname] = {'i_stim': h.Vector(secmd)}
        
                rvec[runname]['v_soma'].record(nrn_cell.soma(0.5)._ref_v)
                rvec[runname]['i_inj'].record(istim[name]._ref_i)
                rvec[runname]['time'].record(h._ref_t)
                # connect current command vector
                h.dt = dt
                h.celsius = temperature
                vec[runname]['i_stim'].play(istim[name]._ref_i, h.dt, 0, sec=nrn_cell.soma)

                nrn_cell.cell_initialize()
                self.custom_init()
                h.t = 0.
                h.tstop = tend
                while h.t < h.tstop:
                    h.fadvance()
                # h.t = 0.
                # h.tstop = tend
                # h.batch_save() # save nothing
                # h.batch_run(h.tstop, h.dt, "v.dat")
                
                pl[name].plot(np.array(rvec[runname]['time']), np.array(rvec[runname]['v_soma']))

        # get overall Rin, etc; need to initialize all cells
        nrn_cell.cell_initialize()
        for n, name in enumerate(self.celltypes.keys()):
            nrn_cell = netcells[name]
            nrn_cell.vm0 = nrn_cell.soma.v
            pars = nrn_cell.compute_rmrintau(auto_initialize=False)
            print(u'{0:>14s} [{1:>14s}]   *** Rin = {2:6.1f} M\u03A9  \u03C4 = {3:6.1f} ms   Vm = {4:6.1f} mV'.
                format(nrn_cell.status['name'], name, pars['Rin'], pars['tau'], pars['v']))


if __name__ == "__main__":
    t=Toy()
    t.run()
    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()
