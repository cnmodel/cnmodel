# -*- encoding: utf-8 -*-
from __future__ import print_function
import sys
import numpy as np
import pyqtgraph as pg
import neuron

import cnmodel
import cnmodel.cells as cells
from cnmodel.protocols import SynapseTest
from cnmodel.util import random_seed, reset
from cnmodel import data

"""
Check that sgc PSDs have correct AMPA and NMDA peak conductances / CV.
"""

def test_sgc_bushy_psd(plot=False):
    sgc_psd_test(cells.Bushy, seed=23572385, tstop=4.0, plot=plot)

    
def test_sgc_tstellate_psd(plot=False):
    sgc_psd_test(cells.TStellate, seed=34754398, plot=plot)


def test_sgc_dstellate_psd(plot=False):
    sgc_psd_test(cells.DStellate, seed=54743998, plot=plot, n_syn=50)


def test_sgc_octopus_psd(plot=False):
    sgc_psd_test(cells.Octopus, seed=54743998, plot=plot, n_syn=50)


def sgc_psd_test(cell_class, seed, plot=False, tstop=5.0, n_syn=20):
    """
    Tests a multisite synapse from the SGC to a target cell.
    The values returned from an actual set of runs of the synapse are compared
    to the expected values in the synapses.py table. This is needed because
    the maximal open probability of the receptor models is not 1, so the maximal
    conductance per receptor needs to be adjusted empirically. If the measured current
    does not match the expected current, then we print a message with the expected value,
    and fail with an assert statment in the test. 
    The measurement itself is made in measure_gmax().
    
    Parameters
    ----------
    cell_class : an instance of the cell class
    seed : int
        random number seed for the call
    plot : boolean (default False)
        plot request, passed to measure_gmax
    tstop : float (default 5.0 ms)
        duration of run for measurement of gmax. Needs to be long enough to find the
        maximum of the EPSC/IPSC.
    n_syn : int (default 20)
        number of synapses to instantiate for testing (to get an average value)
    
    """
    celltyp = cell_class.__name__.lower()
    
    random_seed.set_seed(seed)
    reset(raiseError=False)  # avoid failure because we cannot release NEURON objects completely.
    tsc = cell_class.create(ttx=True)
    (ampa_gmax, nmda_gmax, epsc_cv) = measure_gmax(tsc, n_syn=n_syn, tstop=tstop, plot=plot)
    exp_ampa_gmax = data.get('sgc_synapse', species='mouse', post_type=celltyp, field='AMPA_gmax')[0]
    exp_nmda_gmax = data.get('sgc_synapse', species='mouse', post_type=celltyp, field='NMDA_gmax')[0]
    exp_epsc_cv = data.get('sgc_synapse', species='mouse', post_type=celltyp, field='EPSC_cv')
    ampa_correct = np.allclose(exp_ampa_gmax, ampa_gmax)
    if not ampa_correct:
        AMPAR_gmax = data.get('sgc_synapse', species='mouse', post_type=celltyp, field='AMPAR_gmax')
        ratio = exp_ampa_gmax/ampa_gmax
        print('AMPA Receptor conductance in model should be %.16f (table is %.16f)'
                % (AMPAR_gmax * ratio, AMPAR_gmax))
    nmda_correct = np.allclose(exp_nmda_gmax, nmda_gmax)
    if not nmda_correct:
        NMDAR_gmax = data.get('sgc_synapse', species='mouse', post_type=celltyp, field='NMDAR_gmax')
        ratio = exp_nmda_gmax/nmda_gmax
        print('ratio: ', ratio, exp_nmda_gmax, nmda_gmax)
        print('NMDA Receptor conductance in model should be %.16f (table is %.16f)'
                % (NMDAR_gmax * ratio, NMDAR_gmax))
    cv_correct = (abs(exp_epsc_cv / epsc_cv - 1.0) < 0.1)
    print ('cv_correct: ', cv_correct)
    if not cv_correct:
        ratio = exp_epsc_cv/epsc_cv
        print('CV Receptor in synapses.py model should be %.6f (measured = %.6f; table = %.6f)'
                % (epsc_cv * ratio, epsc_cv, exp_epsc_cv))
        print ((abs(exp_epsc_cv / (epsc_cv * ratio) - 1.0) < 0.1))
    assert cv_correct
    assert ampa_correct and nmda_correct


def measure_gmax(cell, n_syn=20, tstop=5.0, plot=False):
    sgc = cells.SGC.create()
    prot = SynapseTest()
    
    # Connect 20 synapses and stimulate once each
    # Temp is 33 C and vm=+40 to match Cao & Oertel 2010
    prot.run(sgc.soma, cell.soma, n_synapses=n_syn, temp=33.0, dt=0.025, vclamp=40.0,
             tstop=tstop, stim_params={'NP': 1, 'delay': 0.1})
    
    # For each synapse:
    #    * Add up ampa and nmda conductances across all sites
    #      (although SGC-TS synapses currently have only 1 site)
    #    * Keep track of the maximum ampa and nmda conductance
    if plot:
        global plt
        plt = pg.plot()
    ampa_gmax = []
    nmda_gmax = []
    epsc_gmax = []
    for syn in prot.synapses:
        ampa = np.zeros_like(syn.psd.get_vector('ampa', 'g'))
        nmda = ampa.copy()
        ampa_po = ampa.copy()
        for i in range(syn.psd.n_psd):
            
            ampa += syn.psd.get_vector('ampa', 'g', i)*1e-3  # convert pS from mechanism to nS
            nmda += syn.psd.get_vector('nmda', 'g', i)*1e-3
            if nmda[-1] - nmda[-2] > 0.001:
                raise Exception("Did not reach nmda gmax; need longer run.")
        amax = ampa.max()
        ampa_gmax.append(amax)
        nmda_gmax.append(nmda.max())
        epsc_gmax.append((ampa+nmda).max())
        tb = np.linspace(0., len(ampa)*prot.dt, len(ampa))
        if plot:
            plt.plot(tb, ampa, pen='g')
            plt.plot(tb, nmda, pen='y')
    return (np.mean(ampa_gmax), np.mean(nmda_gmax), 
            np.std(epsc_gmax) / np.mean(epsc_gmax))


    
if __name__ == '__main__':
    if len(sys.argv[0]) > 1:
        testcell = sys.argv[1]
    if testcell not in ['bushy', 'tstellate', 'dstellate', 'octopus', 'all']:
        print ('PSD test for cell type %s is not yet supported.' % testcell)
        exit(1)
    else:
        if testcell in ['bushy']:
            test_sgc_bushy_psd(plot=True)
        if testcell in ['tstellate']:
            test_sgc_tstellate_psd(plot=True)
        if testcell in ['dstellate']:
            test_sgc_dstellate_psd(plot=True)
        if testcell in ['octopus']:
            test_sgc_octopus_psd(plot=True)
        if testcell in ['all']:
            test_sgc_bushy_psd(plot=True)
            test_sgc_tstellate_psd(plot=True)
            test_sgc_dstellate_psd(plot=True)
            test_sgc_octopus_psd(plot=True)
            
#    pg.show()
    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()
    