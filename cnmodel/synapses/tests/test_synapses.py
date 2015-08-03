"""
Create presynaptic and postsynaptic neurons, automatically connect them with 
a synapse, stimulate the presynaptic cell, and analyze the resulting PSCs
in the postsynaptic cell.
"""
import faulthandler
faulthandler.enable()
import os, pickle, pprint
import numpy as np
import neuron

import cnmodel
import cnmodel.cells as cells
from cnmodel.util import UserTester
from cnmodel.protocols import SynapseTest
from cnmodel.util import random, reset

#
# Synapse tests
#
def test_sgc_bushy():
    SynapseTester('sgc', 'bushy')

def test_sgc_tstellate():
    SynapseTester('sgc', 'tstellate')

def test_sgc_tstellate2():  # again to test RNG stability
    SynapseTester('sgc', 'tstellate')

def test_sgc_dstellate():
    SynapseTester('sgc', 'dstellate')

def test_dstellate_bushy():
    SynapseTester('dstellate', 'bushy')

def test_dstellate_tstellate():
    SynapseTester('dstellate', 'tstellate')

def test_dstellate_dstellate():
    SynapseTester('dstellate', 'dstellate')



#
# Supporting functions
#
convergence = {
    'sgc': {'bushy': 3, 'tstellate': 6, 'dstellate': 10, 'dstellate_eager': 10},
    'dstellate': {'bushy': 10, 'tstellate': 15, 'dstellate': 5},
    }


def make_cell(typ):
    if typ == 'sgc':
        cell = cells.SGC.create()
    elif typ == 'tstellate':
        cell = cells.TStellate.create(debug=True, ttx=False)
    elif typ == 'dstellate': # Type I-II Rothman model, similiar excitability (Xie/Manis, unpublished)
        cell = cells.DStellate.create(model='RM03', debug=True, ttx=False)
    elif typ == 'dstellate_eager': # From Eager et al.
        cell = cells.DStellate.create(model='Eager', debug=True, ttx=False)
    elif typ == 'bushy':
        cell = cells.Bushy.create(debug=True, ttx=True)
    else:
        raise ValueError("Unknown cell type '%s'" % typ)
    return cell


class SynapseTester(UserTester):
    def __init__(self, pre, post):
        self.st = None
        UserTester.__init__(self, "%s_%s" % (pre, post), pre, post)
    
    def run_test(self, pre, post):
        # Make sure no objects are left over from previous tests
        reset()
        
        # seed random generator using the name of this test
        seed = "%s_%s" % (pre, post)
        
        pre_cell = make_cell(pre)
        post_cell = make_cell(post)
        
        n_term = convergence.get(pre, {}).get(post, None)
        if n_term is None:
            n_term = 1
        st = SynapseTest()
        st.run(pre_cell.soma, post_cell.soma, n_term, seed=seed)
        if self.audit:
            st.show()
        
        info = dict(
            rel_events=st.release_events(),
            rel_timings=st.release_timings(),
            open_prob=st.open_probability(),
            event_analysis=st.analyze_events(),
            )
        self.st = st
        
        #import weakref
        #global last_syn
        #last_syn = weakref.ref(st.synapses[0].terminal.relsi)
        
        
        
        return info
    
    def assert_test_info(self, *args, **kwds):
        try:
            super(SynapseTester, self).assert_test_info(*args, **kwds)
        finally:
            if self.st is not None:
                self.st.hide()
    
