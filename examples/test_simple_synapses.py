import pyqtgraph as pg
from cnmodel.protocols.simple_synapse_test import SimpleSynapseTest
from cnmodel import cells
from cnmodel.synapses import Synapse


import sys

def runtest():
    if len(sys.argv) < 3:
        print "Usage:  python test_synapses.py <pre_celltype> <post_celltype>"
        print "   Supported cell types: sgc, bushy, tstellate, dstellate"
        sys.exit(1)

    convergence = {
        'sgc': {'bushy': 1, 'tstellate': 6, 'dstellate': 10, 'dstellate_eager': 10},
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
        else:
            raise ValueError("Unknown cell type '%s'" % cellType)
        c.append(cell)

    preCell, postCell = c
    
    nTerminals = convergence.get(sys.argv[1], {}).get(sys.argv[2], None)
    if nTerminals is None:
        nTerminals = 1
        print "Warning: Unknown convergence for %s => %s, assuming %d" % (sys.argv[1], sys.argv[2], nTerminals)

    if sys.argv[1:3] == ['sgc', 'bushy']:
        niter = 1
    else:
        niter = 20

    st = SimpleSynapseTest()
    st.run(preCell.soma, postCell.soma, iterations=niter)
    st.show()


if __name__ == '__main__':
    runtest()
    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()
