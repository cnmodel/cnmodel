"""
Test synaptic connections between two different cell types. 

Usage:  python test_synapses.py <pre_celltype> <post_celltype>

This script:

1. Creates single pre- and postsynaptic cells
2. Creates multiple synaptic terminals between the two cells.
   (the convergence is hard-coded below).
3. Stimulates the presynaptic cell by current injection.
4. Records and analyzes the resulting post-synaptic events.

This is used mainly to check that the strength, kinetics, and dynamics of 
each synapse type is working as expected. A higher-level approach is
demonstrated in test_populations.py, in which the presynaptic cells are 
automatically generated using expected patterns of connectivity.
"""

import pyqtgraph as pg
from cnmodel.protocols import SynapseTest
from cnmodel import cells
from cnmodel.synapses import Synapse


import sys

def runtest():
    if len(sys.argv) < 3:
        print "Usage:  python test_synapses.py <pre_celltype> <post_celltype>"
        print "   Supported cell types: sgc, bushy, tstellate, dstellate"
        sys.exit(1)


    convergence = {
        'sgc': {'bushy': 1, 'tstellate': 1, 'dstellate': 1, 'dstellate_eager': 10, 'octopus': 10},
        'dstellate': {'bushy': 10, 'tstellate': 15, 'dstellate': 5},
        }



    c = []
    for cellType in sys.argv[1:3]:
        if cellType == 'sgc':
            cell = cells.SGC.create()
        elif cellType == 'tstellate':
            cell = cells.TStellate.create(debug=True, ttx=False)
        elif cellType == 'dstellate': # Type I-II Rothman model, similiar excitability (Xie/Manis, unpublished)
            cell = cells.DStellate.create(model='RM03', debug=True, ttx=False)
        elif cellType == 'dstellate_eager': # From Eager et al.
            cell = cells.DStellate.create(model='Eager', debug=True, ttx=False)
        elif cellType == 'bushy':
            cell = cells.Bushy.create(debug=True, ttx=True)
        elif cellType == 'octopus':
            cell = cells.Octopus.create(debug=True, ttx=True)
        else:
            raise ValueError("Unknown cell type '%s'" % cellType)
        c.append(cell)

    preCell, postCell = c
    
    nTerminals = convergence.get(sys.argv[1], {}).get(sys.argv[2], None)
    if nTerminals is None:
        nTerminals = 1
        print "Warning: Unknown convergence for %s => %s, assuming %d" % (sys.argv[1], sys.argv[2], nTerminals)

    if sys.argv[1:3] == ['sgc', 'bushy']:
        niter = 5
    else:
        niter = 20
    

    st = SynapseTest()
    st.run(preCell.soma, postCell.soma, nTerminals, vclamp=-65., iterations=niter)
    st.show_result()
    st.plots['VPre'].setYRange(-70., 10.)
    st.plots['EPSC'].setYRange(-2.0, 0.5)
    st.plots['latency2080'].setYRange(0., 1.0)
    st.plots['halfwidth'].setYRange(0., 1.0)
    st.plots['RT'].setYRange(0., 0.2)
    st.plots['latency'].setYRange(0., 1.0)
    st.plots['latency_distribution'].setYRange(0., 1.0)
    return st  # need to keep st alive in memory

if __name__ == '__main__':
    st = runtest()
    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()
