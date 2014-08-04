import os, pickle, pprint
import numpy as np
import neuron

import nrnlibrary
import nrnlibrary.cells as cells
from nrnlibrary.protocols import IVCurve

#
# Consistency tests for synaptic connections between cell types. 
# Basic metrics for psg kinetics, synaptic dynamics / release probability,
# and latency.
#
def test_sgc_bushy():
    pre = cells.SGC.create()
    post = cells.Bushy.create(species='guineapig', type='II')
    CellSynapseTester('sgc_bushy', pre, post)

def test_sgc_tstellate():
    pre = cells.SGC.create()
    post = cells.TStellate.create(species='guineapig', type='I-c')
    assert_test_info('sgc_tstellate', pre, post)



#
# Supporting functions
#


class CellSynapseTester(UserTester):
    def run_test(self, pre_cell, post_cell):
        # run I/V test on cell
        iv = IVCurve()
        iv.run(cell.i_test_range, cell)
        iv.show(cell)
        
        info = dict(
            icmd=iv.current_cmd,
            spikes=iv.spike_times(),
            rmp=iv.rest_vm(),
            rm=iv.input_resistance(),
            vpeak=iv.peak_vm(),
            vss=iv.steady_vm(),
            )
        self.iv = iv
        return info
    
    def assert_test_info(self, *args, **kwds):
        try:
            super(CellTester, self).assert_test_info(*args, **kwds)
        finally:
            self.iv.win.hide()
    

