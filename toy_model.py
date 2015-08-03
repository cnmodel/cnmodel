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
try:
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
    def __init__(self):
        super(Toy, self).__init__()
        print 'made Toy'

    def current_name(self, name, n):
        return '%.3f' % self.celltypes[name][2][n]

    def getname(self, cell, ninj):
        name = cell.status['name'] + ', ' + cell.status['type']
        iname = self.current_name(name, ninj)
        nname = name + iname
        return name, nname

    def run(self):
        print 'running'
        sre = re.compile('(?P<cell>\w+)(?:[, ]*)(?P<type>[\w-]*)')  # regex for keys in cell types
        self.celltypes = OrderedDict([('Bushy, II', (cells.Bushy, "II", (-0.5,0.5))), ('Bushy, II-I', (cells.Bushy, "II-I", (-0.5,0.5))),
                                ('Octopus, II-o', (cells.Octopus, 'II-o', (-2.5, 2.5))),
                                ('TStellate, I-c', (cells.TStellate, "I-c", (-0.15, 0.15))), ('TStellate, I-t', (cells.TStellate, "I-t", (-0.15, 0.15))),
                                ('DStellate, I-II', (cells.DStellate, 'I-II', (-0.2, 0.2))),
                                ('Pyramidal, I', (cells.Pyramidal, 'I', (-0.1, 0.06))), ('Cartwheel, I', (cells.Cartwheel, 'I', (-0.12, 0.1))),
                                ('SGC, bm', (cells.SGC, "bm", (-0.2, 0.2))),
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
        tend = stim['delay'] + stim['dur'] + 20

        netcells = []
        for c in self.celltypes.keys():
            g = sre.match(c)
            cellname = g.group('cell')
            print cellname
            type = g.group('type')
            if g.group('type') == '':
                netcells.append(self.celltypes[c][0].create())
            else:
                netcells.append(self.celltypes[c][0].create(type=type))
        istim = {}

        #
        # build plotting area
        #
        app = pg.mkQApp()
        win = pg.GraphicsWindow()
        self.win=win
        win.resize(800, 600)
        pl = {}
        rvec = {}
        vec = {}
        ncells = len(self.celltypes.keys())
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
        # initialized some recording vectors
        for n, name in enumerate(self.celltypes.keys()):
            for ninj in range(len(self.celltypes[name][2])):
                iname = self.current_name(name, ninj)
                rvec[name+iname] = {'v_soma': h.Vector(), 'i_inj': h.Vector(), 'time': h.Vector()}
                vec[name+iname] = {'istim': h.Vector()}
       #
       # Build injection arrays for each stimulus level
       #
        for cell in netcells:
            name = cell.status['name'] + ', ' + cell.status['type']
            print 'building for: ', name
            for ninj in range(len(self.celltypes[name][2])):
                name, nname  = self.getname(cell, ninj)
                stim['amp'] =  self.celltypes[name][2][ninj]  # currents[name]
                (secmd, maxt, tstims) = make_pulse(stim)
                istim[name] = h.iStim(0.5, sec=cell.soma)
                istim[name].delay = 5.
                istim[name].dur = 1e9 # these actually do not matter...
                istim[name].iMax = 0.
                vec[nname]['i_stim'] = h.Vector(secmd)

        for ninj in range(2):  # all have 2 currents ...
        #
        # initialize all cells
            for cell in netcells:
                name, nname  = self.getname(cell, ninj)
                cell.cell_initialize()
                rvec[nname]['v_soma'].record(cell.soma(0.5)._ref_v)
                rvec[nname]['i_inj'].record(istim[name]._ref_i)
                rvec[nname]['time'].record(h._ref_t)
                # connect current command vector
                vec[nname]['i_stim'].play(istim[name]._ref_i, h.dt, 0, sec=cell.soma)

            h.dt = dt
            h.tstop = tend
            self.custom_init()
            # Now run all cells at one time for the selected individual current levels
            while h.t < h.tstop:
                h.fadvance()
            for cell in netcells:  # plot data from each cell from this run
                (name, nname) = self.getname(cell, ninj)
                pl[name].plot(np.array(rvec[nname]['time']), np.array(rvec[nname]['v_soma']))

        # get overall Rin, etc; need to initialize all cells first
        self.custom_init()
        for cell in netcells:
            cell.vm0 = cell.soma.v
            Rin, tau, v = cell.measure_rintau(auto_initialize=False)
            print(u'{0:>14s}   *** Rin = {1:6.1f} M\u03A9  \u03C4 = {2:6.1f} ms   Vm = {3:6.1f} mV'.format(cell.status['name'], Rin, tau, v))


if __name__ == "__main__":
    t=Toy()
    t.run()
    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()
