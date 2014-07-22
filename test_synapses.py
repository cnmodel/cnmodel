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
    TargetCell = cells.TStellate(debug=True, ttx=False) # make a postsynaptic cell
    nANTerminals = 6
elif cellType == 'dstellate': # Type I-II Rothman model, similiar excitability (Xie/Manis, unpublished)
    TargetCell = cells.DStellate(debug=True, ttx=False) # make a postsynaptic cell
    nANTerminals = 12
elif cellType == 'dstellate_eager': # From Eager et al.
    TargetCell = cells.DStellateEager(debug=True, ttx=False) # make a postsynaptic cell
    nANTerminals = 12
elif cellType == 'bushy':
    TargetCell = cells.Bushy(debug=True, ttx=True) # make a postsynaptic cell
    nANTerminals = 3
else:
    raise ValueError("Unknown cell type '%s'" % cellType)

synapses = [Synapse() for i in range(nANTerminals)]
st = SynapseTest()
st.run(TargetCell, synapses)
st.analyze()

