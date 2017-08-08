"""
Test connection between two cell populations.

Usage: python test_populations.py <pre_celltype> <post_celltype>

This script: 

1. Creates two cell populations (pop1 and pop2)
2. Connects pop1 => pop2
3. Instantiates a single cell in pop2
4. Automatically generates presynaptic cells and synapses from pop1
5. Stimulates presynaptic cells and records postsynaptically

This is a high-level approach to generating networks in that the supporting
cells (those in pop1) are created automatically based on expected patterns
of connectivity in the cochlear nucleus. A lower-level approach is demonstrated
in test_synapses.py, in which the individual pre- and postsynaptic cells are
manually created and connected.
"""
from cnmodel import populations
from cnmodel.protocols import PopulationTest
import pyqtgraph as pg
import sys

def testpopulation():
    if len(sys.argv) < 3:
        print "Usage:  python test_populations.py <pre_celltype> <post_celltype>"
        sys.exit(1)

    pop_types = {
        'sgc': populations.SGC,
        'bushy': populations.Bushy,
        'tstellate': populations.TStellate,
        'dstellate': populations.DStellate,
        'pyramidal': populations.Pyramidal,
        'tuberculoventral': populations.Tuberculoventral,
        }

    pops = []
    for cell_type in sys.argv[1:3]:
        if cell_type not in pop_types:
            print '\nUnsupported cell type: "%s". Options are %s' % (cell_type, pop_types.keys())
            sys.exit(-1)
        pops.append(pop_types[cell_type]())


    pt = PopulationTest()
    pt.run(pops)
    pt.show()

    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()

if __name__ == '__main__':
    testpopulation()
