import os, pickle, pprint
import numpy as np
import neuron

import cnmodel
import cnmodel.cells as cells
from cnmodel.util import UserTester, reset
from cnmodel.protocols import IVCurve

#
# Cell-type tests
#

def test_bushy():
    reset()
    cell = cells.Bushy.create(species='guineapig', type='II')
    CellTester('bushy_guineapig-typeII', cell)

def test_bushy21():
    reset()
    cell = cells.Bushy.create(species='guineapig', type='II-I')
    CellTester('bushy_guineapig-typeII-I', cell)

def test_tstellate():
    reset()
    cell = cells.TStellate.create(species='guineapig', type='I-c')
    CellTester('tstellate_guineapig-typeI-c', cell)

def test_tstellatet():
    reset()
    cell = cells.TStellate.create(species='guineapig', type='I-t')
    CellTester('tstellate_guineapig-typeI-t', cell)

def test_dstellate():
    reset()
    cell = cells.DStellate.create(species='guineapig', type='I-II')
    CellTester('dstellate_guineapig-typeI-II', cell)

def test_octopus():
    reset()
    cell = cells.Octopus.create(species='guineapig', type='II-o')
    CellTester('octopus_guineapig-typeII-o', cell)

def test_pyramidal():
    reset()
    cell = cells.Pyramidal.create(species='rat', type='I')
    CellTester('pyramidal_rat_I', cell)

def test_cartwheel():
    reset()
    cell = cells.Cartwheel.create(species='rat', type='I')
    CellTester('cartwheel_rat_I', cell)

def test_sgc_basal_middle():
    reset()
    cell = cells.SGC.create(species='mouse', type='bm')
    CellTester('SGC_rat_bm', cell)

def test_sgc_apical():
    reset()
    cell = cells.SGC.create(species='mouse', type='a')
    CellTester('SGC_rat_a', cell)



#
# Supporting functions
#

class CellTester(UserTester):
    data_dir = 'cell_data'
    
    def run_test(self, cell):
        # run I/V test on cell
        iv = IVCurve()
        self.iv = iv
        iv.run(cell.i_test_range, cell)
        if self.audit:
            iv.show(cell)
        
        info = dict(
            icmd=iv.current_cmd,
            spikes=iv.spike_times(),
            rmp=iv.rest_vm(),
            rm_taum=iv.input_resistance_tau(),
            vpeak=iv.peak_vm(),
            vss=iv.steady_vm(),
            rmrintau=cell.compute_rmrintau(),
            )
        return info
    
    def assert_test_info(self, *args, **kwds):
        try:
            super(CellTester, self).assert_test_info(*args, **kwds)
        finally:
            if hasattr(self, 'iv') and hasattr(self.iv, 'win'):
                self.iv.win.hide()
    

#def result_file(key):
    #"""
    #Return a file name to be used for storing / retrieving test results
    #given *key*.
    #"""
    #path = os.path.dirname(__file__)
    #return os.path.join(path, 'cell_data', key + '.pk')

#def load_cell_info(key):
    #"""
    #Load prior test results for *key*.
    #If there are no prior results, return None.
    #"""
    #fn = result_file(key)
    #if os.path.isfile(fn):
        #return pickle.load(open(fn, 'rb'))
    #return None

#def save_cell_info(info, key):
    #"""
    #Store test results for *key*.
    #"""
    #fn = result_file(key)
    #dirname = os.path.dirname(fn)
    #if not os.path.isdir(dirname):
        #os.mkdir(dirname)
    #pickle.dump(info, open(fn, 'wb'))
    
    
#def CellTester(key):
    #"""
    #Test *cell* and raise exception if the results do not match prior
    #data.
    #"""
    #audit = cnmodel.AUDIT_TESTS
    
    ## run I/V test on cell
    #iv = IVCurve()
    #iv.run(cell.i_test_range, cell)
    #iv.show(cell)
    
    #try:
        #info = dict(
            #icmd=iv.current_cmd,
            #spikes=iv.spike_times(),
            #rmp=iv.rest_vm(),
            #rm=iv.input_resistance(),
            #vpeak=iv.peak_vm(),
            #vss=iv.steady_vm(),
            #)

        #expect = load_cell_info(key)
        
        #if expect is not None:
            
            ## Check test structures are the same
            #assert len(info) == len(expect)
            #for k in info:
                #assert k in expect
                
            ## Check data matches
            #for k in info:
                #if isinstance(info[k], list):
                    #assert len(info[k]) == len(expect[k])
                    #for i in range(len(info[k])):
                        #assert np.allclose(info[k][i], expect[k][i])
                #else:
                    #assert np.allclose(info[k], expect[k])
        #else:
            #if not audit:
                #raise Exception("No prior test results for cell type '%s'. "
                                #"Run test.py --audit store new test data." % key)
            
            #print "\n=== New test results for %s: ===\n" % key
            #pprint.pprint(info)
            #print "Store new test results? [y/n]",
            #yn = raw_input()
            #if yn.lower().startswith('y'):
                #save_cell_info(info, key)
            #else:
                #raise Exception("Rejected test results for '%s'" % key)
    #finally:
        #iv.win.hide()
    
