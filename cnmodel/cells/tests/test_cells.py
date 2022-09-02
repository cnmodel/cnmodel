import os, pickle, pprint
import numpy as np
import neuron

import cnmodel
import cnmodel.cells as cells
from cnmodel.util import UserTester, reset
from cnmodel.protocols import IVCurve

"""
Cell-type tests
"""

def test_bushy():
    reset(raiseError=False)
    cell = cells.Bushy.create(species='guineapig', modelType='II')
    CellTester('bushy_guineapig-typeII', cell)

def test_bushy21():
    reset(raiseError=False)
    cell = cells.Bushy.create(species='guineapig', modelType='II-I')
    CellTester('bushy_guineapig-typeII-I', cell)

def test_bushy_mouse():
    reset(raiseError=False)
    cell = cells.Bushy.create(species='mouse', modelType='II', nach='na')
    CellTester('bushy-mouse-typeII', cell)
    
def test_tstellate():
    reset(raiseError=False)
    cell = cells.TStellate.create(species='guineapig', modelType='I-c')
    CellTester('tstellate_guineapig-typeI-c', cell)

def test_tstellate_mouse():
    reset(raiseError=False)
    cell = cells.TStellate.create(species='mouse', modelType='I-c')
    CellTester('tstellate_mouse-typeI-c', cell)
    
def test_tstellatet():
    reset(raiseError=False)
    cell = cells.TStellate.create(species='guineapig', modelType='I-t')
    CellTester('tstellate_guineapig-typeI-t', cell)

def test_dstellate():
    reset(raiseError=False)
    cell = cells.DStellate.create(species='guineapig', modelType='I-II')
    CellTester('dstellate_guineapig-typeI-II', cell)

def test_dstellate_mouse():
    reset(raiseError=False)
    cell = cells.DStellate.create(species='mouse', modelType='I-II')
    CellTester('dstellate_mouse-typeI-II', cell)
    
def test_octopus():
    reset(raiseError=False)
    cell = cells.Octopus.create(species='guineapig', modelType='II-o')
    CellTester('octopus_guineapig-typeII-o', cell)

def test_octopus_mouse():
    reset(raiseError=False)
    cell = cells.Octopus.create(species='mouse', modelType='II-o')
    CellTester('octopus_mouse-typeII-o', cell)

def test_pyramidal():
    reset(raiseError=False)
    cell = cells.Pyramidal.create(species='rat', model='POK', modelType='I')
    CellTester('pyramidal_rat_I', cell)

def test_pyramidal_ceballos():
    reset(raiseError=False)
    cell = cells.PyramidalCeballos.create(species='mouse', model='Ceballos', modelType='I')
    CellTester('pyramidal_mouse_I', cell)
    
def test_tuberculoventral():
    reset(raiseError=False)
    cell = cells.Tuberculoventral.create(species='mouse', modelType='TVmouse')
    CellTester('tuberculoventral_mouse_I', cell)

def test_cartwheel():
    reset(raiseError=False)
    cell = cells.Cartwheel.create(species='mouse', modelType='I')
    CellTester('cartwheel_rat_I', cell)

def test_sgc_basal_middle():
    reset(raiseError=False)
    cell = cells.SGC.create(species='mouse', modelType='sgc-bm')
    CellTester('SGC_rat_bm', cell)

def test_sgc_apical():
    reset(raiseError=False)
    cell = cells.SGC.create(species='mouse', modelType='sgc-a')
    CellTester('SGC_rat_a', cell)

def test_mso():
    reset(raiseError=False)
    cell = cells.MSO.create(species='guineapig', modelType='principal')
    CellTester('mso_guineapig-principal', cell)

#
# Supporting functions
#

class CellTester(UserTester):
    data_dir = 'cell_data'
    
    def run_test(self, cell):
        # run I/V test on cell
        V0 = cell.find_i0(showinfo=True)
        rmrintau = cell.compute_rmrintau(auto_initialize=False, vrange=None)
        iv = IVCurve()
        self.iv = iv
        iv.run(cell.i_test_range, cell)
        if self.audit:
            iv.show(cell)
        
        info = dict(
            temp=iv.temp,
            icmd=iv.current_cmd,
            spikes=iv.spike_times(),
            rmp=iv.rest_vm(),
            rm_taum=iv.input_resistance_tau(),
            vpeak=iv.peak_vm(),
            vss=iv.steady_vm(),
            rmrintau=rmrintau,
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
    

# The following is superseeded by the built in unit tests.
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
    
