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
elif cellType == 'dstellate': # similar to t-stellate in abasence of other data
    TargetCell = cells.DStellate(debug=True, ttx=False) # make a postsynaptic cell
elif cellType == 'bushy':
    TargetCell = cells.Bushy(debug=True, ttx=True) # make a postsynaptic cell
else:
    raise ValueError("Unknown cell type '%s'" % cellType)

synapse = Synapse()
st = SynapseTest()
st.run(TargetCell, synapse)
st.analyze()

