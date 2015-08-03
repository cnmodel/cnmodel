"""
1. Create two cell populations (pop1 and pop2)
2. Connect pop1 => pop2
3. Instantiate a single cell in pop2
4. Automatically generate presynaptic cells and synapses from pop1
5. Stimulate presynaptic cells and record postsynaptically
"""
from cnmodel import populations
from cnmodel.protocols import PopulationTest

import sys
if len(sys.argv) < 3:
    print "Usage:  python test_circuits.py <pre_celltype> <post_celltype>"
    sys.exit(1)

pop_types = {
    'sgc': populations.SGC,
    'bushy': populations.Bushy,
    'tstellate': populations.TStellate,
    'dstellate': populations.DStellate,
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

