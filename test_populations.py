"""
1. Create two cell populations (pop1 and pop2)
2. Connect pop1 => pop2
3. Instantiate a single cell in pop2
4. Automatically generate presynaptic cells and synapses from pop1
5. Stimulate presynaptic cells and record postsynaptically
"""
from nrnlibrary import populations

import sys
if len(sys.argv) < 3:
    print "Usage:  python test_circuits.py <pre_celltype> <post_celltype>"
    sys.exit(1)

pop_types = {
    'sgc': populations.SGC,
    'bushy': populations.Bushy,
    #'tstellate': tstellate,
    #'dstellate': dstellate,
    }

pops = []
for cell_type in sys.argv[1:3]:
    if cell_type not in pop_types:
        print 'Unsupported cell type: "%s". Options are %s' % (cell_type, pops.keys())
        sys.exit(-1)
    pops.append(pop_types[cell_type]())

pre_pop, post_pop = pops
pre_pop.connect(post_pop)

# start with one cell, selected from the user-selected population, that has
# a cf close to 4kHz
cell = post_pop.select(1, cf=4000, create=True)
post_pop.resolve_inputs(depth=1)


pt = PopulationTest()
pt.run(pops)
pt.show()

