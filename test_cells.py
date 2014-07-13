import argparse
import os, sys
from neuron import h

import nrnlibrary
import nrnlibrary.cells as cells
from nrnlibrary.protocols import IVCurve, VCCurve

debugFlag = True
parser = argparse.ArgumentParser(description=('test_cells.py:',
' Biophysical representations of neurons (mostly auditory), test file'))

cclamp = False
cellinfo = {'types': ['bushy', 'stellate', 'steldend', 'dstellate', 'sgc',
                        'cartwheel', 'pyramidal', 'octopus'],
            'configs': ['std, ''waxon', 'dendrite'],
            'nav': ['std', 'jsrna', 'nav11'],
            'species': ['guineapig', 'guineapig-bushy-II-I',
                                'guineapig-bushy-II', 'rat', 'mouse'],
            'pulse': ['step', 'pulse']}
ccivrange = {'bushy': (-0.5, 0.5, 0.05),
            'stellate': (-.2, 0.2, 0.05),
            'steldend': (-1.0, 1.0, 0.1),
            'dstellate': (-0.25, 1.0, 0.05),
            'sgc:': (-0.5, 0.5, 0.05),
            'cartwheel': (-0.5, 0.5, 0.05),
            'pyramidal': (-1., 1., 0.1),
            'octopus': (-2., 2., 0.2)}
# scales holds some default scalint to use in the cciv plots
# argument is {cellname: (xmin, xmax, IVymin, IVymax, FIspikemax,
# offset(for spikes), crossing (for IV) )}
## the "offset" refers to setting the axes back a bit
scale = {'bushy': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
                'crossing', [0, -60]),
            'stellate': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
                'crossing', [0, -60]),
            'steldend': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
                'crossing', [0, -60]),
            'dstellate': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
                'crossing', [0, -60]),
            'sgc:': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
                'crossing', [0, -60]),
            'cartwheel': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
                'crossing', [0, -60]),
            'pyramidal': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
                'crossing', [0, -60]),
            'octopus': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
                'crossing', [0, -60])}
ax = None
h.celsius = 32
parser.add_argument('celltype', action='store')
parser.add_argument('species', action='store', default='guineapig')
    # species is an optional option....
parser.add_argument('-c', action="store", dest="configuration",
    default='std', help=("Set axon config: %s " %
        [cfg for cfg in cellinfo['configs']]))
parser.add_argument('--nav', action="store", dest="nav", default="std",
    help=("Choose sodium channel: %s " % [ch for ch in cellinfo['nav']]))
parser.add_argument('-p', action="store", dest="pulsetype", default="step",
    help=("Set CCIV pulse to step or repeated pulse"))
clampgroup = parser.add_mutually_exclusive_group()
clampgroup.add_argument('--vc', action='store_true',
    help="Run in voltage clamp mode")
clampgroup.add_argument('--cc', action='store_true',
    help="Run in current clamp mode")
clampgroup.add_argument('--demo', action='store_true',
    help="Run in  voltage clamp demo")
args = parser.parse_args()
print args.celltype
if args.celltype in cellinfo['types']:
    print 'cell: %s is ok' % args.celltype
else:
    print 'cell: %s is not in our list of cell types' % (args.celltype)
    print 'celltypes: ', cellinfo['types']
    sys.exit(1)

path = os.path.dirname(nrnlibrary.__file__)
#h.nrn_load_dll(os.path.join(path, 'i386/special'))
h.load_file("stdrun.hoc")
h.load_file(os.path.join(path, "custom_init.hoc"))
# replace init with one that gets closer to steady state

print 'configuration: ', args.configuration
sites = None
if args.pulsetype == 'step':
    ptype = None
else:
    ptype = 'pulses'
if args.configuration in cellinfo['configs']:
    print 'Configuration %s is ok' % args.configuration


if args.celltype == 'sgc':
    (cell, sgcaxon) = cells.SGC(debug=debugFlag, species='mouse',
    nach = 'nav11', chlist = ['ih'])

#
# T-stellate tests
#
elif (args.celltype == 'stellate' and args.nav == 'nav11'
        and args.species == 'guineapig'):
    cell = cells.TStellateNav11(debug=debugFlag)
elif (args.celltype == 'steldend'):
    cell = cells.TStellateNav11(debug=debugFlag, dend=True,
                                ttx=False, cs=False)
elif (args.celltype == 'stellate' and args.nav == 'nav11'
        and args.species == 'mouse'):
    cell = cells.TStellate(species=args.species, nav11=True, 
                            debug=debugFlag)
elif args.celltype == 'stellate' and args.nav == 'std':
    cell = cells.TStellateFast(debug=debugFlag)

#
# Bushy tests
#
elif (args.celltype == 'bushy' and args.configuration == 'waxon'):
    cell = cells.BushyWithAxon(debug=debugFlag)
elif args.celltype == 'bushy' and args.configuration == 'std':
    cell = cells.Bushy(debug=debugFlag)
    sites = [cell.soma, None, None, None]
elif args.celltype == 'bushy' and args.configuration == 'dendrite':
    dendriteFlag = True
    cell = cells.Bushy(debug=debugFlag, dendrite=dendriteFlag)
    sites = [cell.soma, cell.maindend, cell.secdend[0]]
    
#
# D-stellate tests
#
elif args.celltype == 'dstellate':
    cell = cells.DStellate(debug=debugFlag)
    
else:
    print ("Cell Type %s and configurations nav=%s or config=%s are not available" % (args.celltype, args.nav, args.configuration))
    sys.exit(1)
#    seg = cell()

print("Cell model: %s" % cell.__class__.__name__)
print(cell.__doc__)

import pyqtgraph as pg
app = pg.mkQApp()

#
# define the current clamp electrode and default settings
#
if args.cc is True:
    #run_iv(ccivrange[args.celltype], cell,
        #sites=sites, reppulse=ptype)
    iv = IVCurve()
    iv.run(ccivrange[args.celltype], cell, sites=sites, reppulse=ptype)
    iv.show()
elif args.vc is True:
    vc = VCCurve()
    vc.run((-120, -40, 5), cell)
    vc.show()
elif args.demo is True:
    run_democlamp(cell, dendrites)
else:
    print("Nothing to run. Specify one of --cc, --vc, --democlamp.")
    sys.exit(1)
#-----------------------------------------------------------------------------
#
# If we call this directly, provide a test with the IV function
# - see below to switch cells
#


if sys.flags.interactive == 0:
    pg.QtGui.QApplication.exec_()
