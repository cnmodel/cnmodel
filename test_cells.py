import argparse
import os, sys
from neuron import h

import cnmodel
import cnmodel.cells as cells
from cnmodel.protocols import IVCurve, VCCurve

debugFlag = True
parser = argparse.ArgumentParser(description=('test_cells.py:',
' Biophysical representations of neurons (mostly auditory), test file'))

cclamp = False
cellinfo = {'types': ['bushy', 'stellate', 'steldend', 'dstellate', 'dstellateeager', 'sgc',
                        'cartwheel', 'pyramidal', 'octopus'],
            'configs': ['std, ''waxon', 'dendrite'],
            'nav': ['std', 'jsrna', 'nav11'],
            'species': ['guineapig', 'cat', 'rat', 'mouse'],
            'pulse': ['step', 'pulse']}
# Format for ivranges is list of tuples. This allows finer increments in selected ranges, such as close to rest
ccivrange = {'bushy': [(-0.5, 0.5, 0.025)],
            'stellate': [(-0.2, 0.2, 0.02), (-0.015, 0, 0.005)],
            'steldend': [(-1.0, 1.0, 0.1)],
            'dstellate': [(-0.2, 0.2, 0.0125)],
            'dstellateeager': [(-0.6, 1.0, 0.025)],
            'sgc': [(-0.3, 0.3, 0.01)],
            'cartwheel': [(-0.2, 0.1, 0.02)],
            'pyramidal': [(-0.3, 0.3, 0.025), (-0.040, 0.025, 0.005)],
            'octopus': [(-3., 3., 0.2)],
            }
# scales holds some default scaling to use in the cciv plots
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
            'dstellateeager': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
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
h.celsius = 22
parser.add_argument('celltype', action='store')
parser.add_argument('species', action='store', default='guineapig')
parser.add_argument('--type', action='store', default=None)
parser.add_argument('--temp', action='store', default=22.0,
                    help=("Temp DegC (22 default)"))
    # species is an optional option....
parser.add_argument('-c', action="store", dest="configuration",
    default='std', help=("Set axon config: %s " %
        [cfg for cfg in cellinfo['configs']]))
parser.add_argument('--nav', action="store", dest="nav", default="na",
    help=("Choose sodium channel: %s " % [ch for ch in cellinfo['nav']]))
parser.add_argument('--ttx', action="store_true", dest="ttx", default=False,
    help=("Use TTX (no sodium current"))
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

path = os.path.dirname(cnmodel.__file__)
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

#
# Spiral Ganglion cell tests
#

if args.celltype == 'sgc':
    cell = cells.SGC.create(debug=debugFlag, species=args.species, nach=args.nav, ttx=args.ttx, type=args.type)
#
# T-stellate tests
#
elif args.celltype == 'stellate':
    cell = cells.TStellate.create(debug=debugFlag, species=args.species, nach=args.nav, type=args.type, ttx=args.ttx)
#
# Bushy tests
#
elif args.celltype == 'bushy' and args.configuration == 'waxon':
    cell = cells.Bushy.create(debug=debugFlag, species=args.species, nach=args.nav, type=args.type, ttx=args.ttx)
    cell.add_axon()

elif args.celltype == 'bushy' and args.configuration == 'std':
    cell = cells.Bushy.create(debug=debugFlag, species=args.species, nach=args.nav, type=args.type, ttx=args.ttx)
#
# Ocotpus tests
#
elif args.celltype == 'octopus' and args.configuration == 'std':
    cell = cells.Octopus.create(debug=debugFlag, species=args.species, nach='jsrna', type=args.type, ttx=args.ttx)
#
# D-stellate tests
#
elif args.celltype == 'dstellate':
    cell = cells.DStellate.create(debug=debugFlag, ttx=args.ttx, type=args.type)

elif args.celltype == 'dstellateeager':
    cell = cells.DStellateEager.create(debug=debugFlag, ttx=args.ttx, type=args.type)

#
# DCN pyramidal cell tests
#
elif args.celltype == 'pyramidal':
    cell = cells.Pyramidal.create(debug=debugFlag, ttx=args.ttx, type=args.type)

#
# DCN cartwheel cell tests
#
elif args.celltype == 'cartwheel':
    cell = cells.Cartwheel.create(debug=debugFlag, ttx=args.ttx, type=args.type)

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
    iv.run(ccivrange[args.celltype],  cell, durs=[10., 100., 20.],
           sites=sites, reppulse=ptype, temp=float(args.temp))
    iv.show(cell=cell)
elif args.vc is True:
    vc = VCCurve()
    vc.run((-120, 40, 5), cell)
    vc.show(cell=cell)
elif args.demo is True:
    run_democlamp(cell, dendrites)
else:
    print("Nothing to run. Specify one of --cc, --vc, --democlamp.")
    sys.exit(1)
#-----------------------------------------------------------------------------
#
# If we call this directly, provide a test with the IV function
#


if sys.flags.interactive == 0:
    pg.QtGui.QApplication.exec_()
