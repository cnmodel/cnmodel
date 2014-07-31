from nrnlibrary.protocols import SynapseTest
from nrnlibrary import cells
from nrnlibrary.synapses import Synapse

import sys
if len(sys.argv) < 2:
    print "Usage:  python test_synapses.py <celltype>"
    print "   Supported cell types: bushy, tstellate, dstellate"
    sys.exit(1)

cellType = sys.argv[1]


if cellType == 'tstellate':
    postCell = cells.TStellate.create(debug=True, ttx=False) # make a postsynaptic cell
    nANTerminals = 6
elif cellType == 'dstellate': # Type I-II Rothman model, similiar excitability (Xie/Manis, unpublished)
    postCell = cells.DStellate.create(debug=True, ttx=False) # make a postsynaptic cell
    nANTerminals = 10
elif cellType == 'dstellate_eager': # From Eager et al.
    postCell = cells.DStellateEager.create(debug=True, ttx=False) # make a postsynaptic cell
    nANTerminals = 10
elif cellType == 'bushy':
    postCell = cells.Bushy.create(debug=True, ttx=True) # make a postsynaptic cell
    nANTerminals = 3
else:
    raise ValueError("Unknown cell type '%s'" % cellType)

preCell = cells.DStellate.create()
#preCell = cells.SGC.create()
#synapses = [Synapse() for i in range(nANTerminals)]
st = SynapseTest()
st.run(preCell, postCell, nANTerminals)
st.analyze()

