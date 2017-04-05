#!/usr/bin/python

"""
Brute force parameter space exploration

Two steps:
1. Rin, taum and Vm exploration
2. FI curve exploration (once resting parameters are set)
"""
import os
import numpy as np

from neuron import h

import cnmodel
import cnmodel.cells as cells
from cnmodel.protocols import IVCurve, VCCurve
from cnmodel.util import nstomho
path = os.path.dirname(cnmodel.__file__)
h.load_file("stdrun.hoc")
h.load_file(os.path.join(path, "custom_init.hoc"))

cell = cells.Tuberculoventral.create(debug=False, ttx=False, modelType='I', species='mouse')
print 'setup'
somas = np.arange(12., 50.1, 2.)
leaks = np.arange(0., 30.1, 2)
ihgbar = np.arange(.2, 20., 0.2)
nsomas = len(somas)
nleaks = len(leaks)
nih = len(ihgbar)
slmap = np.zeros((nsomas, nleaks, nih))
cell.vrange=[-100., 100.]
for somasize in somas:
    cell.set_soma_size_from_Cm(somasize)
    for leakbar in leaks:
        cell.soma().leak.gbar = nstomho(leakbar, cell.somaarea)
        for ihbar in ihgbar:
            cell.soma().ihvcn.gbar = nstomho(ihbar, cell.somaarea)
            cell.vm0 = None
            cell.get_mechs(cell.soma)
            cell.vrange=[-100., 100.]
            cell.cell_initialize(vrange=cell.vrange)
            resting = cell.compute_rmrintau(auto_initialize=True)
            if (resting['tau'] > 8.) and (resting['tau'] < 12.) and (-75. < resting['v'] < -50.):
                print ('somasize: %.1f leakbar: %.3f, ihbar: %.3f v: %.2f, Rin: %.1f  tau: %.1f' % 
                (somasize, cell.soma().leak.gbar*1000., cell.soma().ihvcn.gbar*1000., resting['v'], resting['Rin'], resting['tau']))
            #print '   L, h: ', cell.soma().cm, cell.soma().area(), cell.soma().diam, cell.soma().sec.L
        # except:
        #     pass
