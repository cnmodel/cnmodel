"""
Test synaptic connections between two different cell types. 

usage: test_synapses.py [-h] [-t {simple,multisite}] [-c]
                        {sgc,tstellate,dstellate,tuberculoventral}
                        {bushy,tstellate,dstellate,octopus,tuberculoventral,pyramidal}

Compute AN only PSTH in postsynaptic cell

positional arguments:
  {sgc,tstellate,dstellate,tuberculoventral}
                        Select presynaptic cell type
  {bushy,tstellate,dstellate,octopus,tuberculoventral,pyramidal}
                        Select postsynaptic cell type

optional arguments:
  -h, --help            show this help message and exit
  -t {simple,multisite}, --type {simple,multisite}
                        Set synapse type (simple, multisite)
  -c, --convergence     Use convergence = 1 for comparision between simple and
                        multi, instead of default table


This script:

1. Creates single pre- and postsynaptic cells
2. Creates multiple synaptic terminals between the two cells.
   (the convergence is hard-coded below).
3. Stimulates the presynaptic cell by current injection.
4. Records and analyzes the resulting post-synaptic events.

This is used mainly to check that the strength, kinetics, and dynamics of 
each synapse type is working as expected. A higher-level approach is
demonstrated in test_populations.py, in which the presynaptic cells are 
automatically generated using expected patterns of connectivity.
"""

import argparse
import pyqtgraph as pg
from cnmodel.protocols import SynapseTest
from cnmodel import cells
from cnmodel.synapses import Synapse


import sys

def runtest():
    parser = argparse.ArgumentParser(description='Compute AN only PSTH in postsynaptic cell')
    parser.add_argument(type=str, dest='precell', default='sgc',
                        choices = ['sgc', 'tstellate', 'dstellate', 'tuberculoventral'],
                        help='Select presynaptic cell type')
    parser.add_argument(type=str, dest='postcell', default='bushy',
                        choices = ['bushy', 'tstellate', 'dstellate', 'octopus',
                            'tuberculoventral', 'pyramidal'],
                        help='Select postsynaptic cell type')
    parser.add_argument('-t', '--type', type=str, dest='syntype', default='multisite',
                        choices=['simple', 'multisite'],
                        help='Set synapse type (simple, multisite)')
    parser.add_argument('-c', '--convergence', action='store_true', dest='convergence', 
                        help='Use convergence = 1 for comparision between simple and multi, instead of default table')
                        
    args = parser.parse_args()

    precell = args.precell
    postcell = args.postcell
    synapseType = args.syntype
    use_conv_table = args.convergence

# These must be se3t to 1 to match data in original tables. Otherwise, it would be better
# to use the original tables.

    if not use_conv_table:
       convergence = {
            'sgc': {'bushy': 1, 'tstellate': 1, 'dstellate': 1, 'dstellate_eager': 10, 'octopus': 10, 
                'tuberculoventral': 1, 'pyramidal': 1, 'cartwheel': 0},
            'dstellate': {'bushy': 10, 'tstellate': 15, 'dstellate': 5},
            }
    else:
       convergence = {
            'sgc': {'bushy': 1, 'tstellate': 1, 'dstellate': 1, 'dstellate_eager': 1, 'octopus': 1,
                    'tuberculoventral': 1, 'pyramidal': 1, 'cartwheel': 0},
            'dstellate': {'bushy': 1, 'tstellate': 1, 'dstellate': 0, 'tuberculoventral': 1, 'pyramidal': 1, 'cartwheel': 0},
            'tuberculoventral': {'bushy': 1, 'tstellate': 1, 'dstellate': 0,
                    'tuberculoventral': 1, 'pyramidal': 1, 'cartwheel': 0},
            }        

    c = []
    for cellType in [precell, postcell]:
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
        elif cellType == 'tuberculoventral':
            cell = cells.Tuberculoventral.create(debug=True, ttx=True)
        elif cellType == 'pyramidal':
            cell = cells.Pyramidal.create(debug=True, ttx=True)
        elif cellType == 'octopus':
            cell = cells.Octopus.create(debug=True, ttx=True)
        else:
            raise ValueError("Unknown cell type '%s'" % cellType)
        c.append(cell)

    preCell, postCell = c
    
    if not use_conv_table:
        nTerminals = convergence.get(precell, {}).get(postcell, None)
    else:
        nTerminals = 1
        # print("Warning: Unknown convergence for %s => %s, assuming %d" % (precell, postcell, nTerminals))

    if [precell, postcell] == ['sgc', 'bushy']:
        niter = 5
    else:
        niter = 20
    # syntype = 'multisite'
    # if len(sys.argv) > 3:
    #     syntype = sys.argv[3]
    # assert(syntype in ['simple', 'multisite'])
    if synapseType == 'simple':
        niter = 1

    st = SynapseTest()
    st.run(preCell.soma, postCell.soma, nTerminals, vclamp=-65., iterations=niter, synapsetype=synapseType)
    st.show_result()
    st.plots['VPre'].setYRange(-70., 10.)
    st.plots['EPSC'].setYRange(-2.0, 0.5)
    st.plots['latency2080'].setYRange(0., 1.0)
    st.plots['halfwidth'].setYRange(0., 1.0)
    st.plots['RT'].setYRange(0., 0.2)
    st.plots['latency'].setYRange(0., 1.0)
    st.plots['latency_distribution'].setYRange(0., 1.0)
    return st  # need to keep st alive in memory

def main():
    st = runtest()
    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()
    
if __name__ == '__main__':
    main()
