"""
Test synaptic connections between two different cell types. 

Usage:  python compare_simple_multisynapses.py <pre_celltype> <post_celltype>

This script:

1. Creates single pre- and postsynaptic cells
2. Creates a single synaptic terminal between the two cells, using the multisite synapse method.
3. Stimulates the presynaptic cell by current injection.
4. Records and analyzes the resulting post-synaptic events.
5. Repeats 3, 4 100 times to get an average postsynaptic event.
6. stores the resulting waveform in a pandas database

This is used mainly to check that the strength, kinetics, and dynamics of 
each synapse type is working as expected. 
"""

import pyqtgraph as pg
from cnmodel.protocols import SynapseTest
from cnmodel import cells
from cnmodel.synapses import Synapse
import pandas as pd

import sys

args = sys.argv[1:]

convergence = {
    'sgc': {'bushy': 1, 'tstellate': 1, 'dstellate': 1,  'octopus': 1, 'pyramidal': 1, 'tuberculoventral': 1},
    'dstellate': {'bushy': 1, 'tstellate': 1, 'dstellate': 1, 'pyramidal': 1, 'tuberculoventral': 1},
    'tuberculoventral': {'bushy': 1, 'tstellate': 1, 'pyramidal': 1, 'tuberculoventral': 1},
    }

def runtest(args, synapsetype='multisite'):
    if len(args) < 2:
        print("Usage:  python compare_simple_multisynapses.py <pre_celltype> <post_celltype>")
        print("   Supported cell types: sgc, bushy, tstellate, dstellate, octopus, tuberculoventral, pyramidal, cartwheel")
        sys.exit(1)



    celltypes = args[:2]
    
    
    c = []
    for cellType in celltypes:
        if cellType == 'sgc':
            cell = cells.SGC.create()
        elif cellType == 'tstellate':
            cell = cells.TStellate.create(debug=True, ttx=False)
        elif cellType == 'dstellate': # Type I-II Rothman model, similiar excitability (Xie/Manis, unpublished)
            cell = cells.DStellate.create(model='RM03', debug=True, ttx=False)
        # elif cellType == 'dstellate_eager': # From Eager et al.
        #     cell = cells.DStellate.create(model='Eager', debug=True, ttx=False)
        elif cellType == 'bushy':
            cell = cells.Bushy.create(debug=True, ttx=True)
        elif cellType == 'octopus':
            cell = cells.Octopus.create(debug=True, ttx=True)
        elif cellType == 'tuberculoventral':
            cell = cells.Tuberculoventral.create(debug=True, ttx=True)
        elif cellType == 'pyramidal':
            cell = cells.Pyramidal.create(debug=True, ttx=True)
        elif cellType == 'cartwheel':
            cell = cells.Cartwheel.create(debug=True, ttx=True)
        else:
            raise ValueError("Unknown cell type '%s'" % cellType)
        c.append(cell)

    preCell, postCell = c
    
    print('celltypes: ', celltypes)
    nTerminals = convergence.get(args[0], {}).get(args[1], None)
    if nTerminals is None:
        nTerminals = 1
        print("Warning: Unknown convergence for %s => %s, assuming %d" % (args[0], args[1], nTerminals))

    if args[1:3] == ['sgc', 'bushy']:
        niter = 5
    else:
        niter = 20
    assert(synapsetype in ['simple', 'multisite'])
    st = SynapseTest()
    
    st.run(preCell.soma, postCell.soma, nTerminals, vclamp=-65., iterations=niter, synapsetype=synapsetype)
    st.show_result()
    st.plots['VPre'].setYRange(-70., 10.)
    st.plots['EPSC'].setYRange(-2.0, 0.5)
    st.plots['latency2080'].setYRange(0., 1.0)
    st.plots['halfwidth'].setYRange(0., 1.0)
    st.plots['RT'].setYRange(0., 0.2)
    st.plots['latency'].setYRange(0., 1.0)
    st.plots['latency_distribution'].setYRange(0., 1.0)
    return st  # need to keep st alive in memory

def runall():
    st = []
    for pre in convergence.keys():
        for post in convergence[pre]:
#            if pre == 'dstellate' and post == 'dstellate':
            sti = runtest([pre, post], 'simple')
            st.append(sti)
    return(st)

if __name__ == '__main__':
    #st = runtest()
    st = runall()
    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()
