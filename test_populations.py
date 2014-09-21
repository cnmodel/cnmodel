"""
Create a single neuron of specified type, then generate and connect the 
supporting circuitry needed to completely simulate the behavior of that 
neuron.
"""
from nrnlibrary import populations



import sys
if len(sys.argv) < 2:
    print "Usage:  python test_circuits.py <celltype>"
    sys.exit(1)

sgc = populations.SGC()
bushy = populations.Bushy()
tstellate = populations.TStellate()
dstellate = populations.DStellate()

sgc.connect(bushy, tstellate, dstellate)
dstellate.connect(bushy, tstellate, dstellate)
tstellate.connect(bushy, tstellate, dstellate)

pops = {
    'sgc': sgc,
    'bushy': bushy,
    'tstellate': tstellate,
    'dstellate': dstellate,
    }

cellType = sys.argv[1]

# start with one cell
cell = pops[cellType].add_cell()
pops[cellType].resolve_inputs(depth=1)


ct = CircuitTest()
ct.run(pops)
ct.show()

