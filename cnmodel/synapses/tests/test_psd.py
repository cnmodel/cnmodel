import numpy as np
import pyqtgraph as pg
import neuron

import cnmodel
import cnmodel.cells as cells
from cnmodel.protocols import SynapseTest
from cnmodel.util import random, reset
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


def sgc_psd_test(cell_class, seed, plot=False, tstop=5.0, n_syn=20):
    celltyp = cell_class.__name__.lower()
    
    random.set_seed(seed)
    reset()
    tsc = cell_class.create(ttx=True)
    (ampa_gmax, nmda_gmax, epsc_cv) = measure_gmax(tsc, n_syn=n_syn, tstop=tstop, plot=plot)
    
    exp_ampa_gmax = data.get('sgc_synapse', species='mouse', post_type=celltyp, field='AMPA_gmax')[0]
    exp_nmda_gmax = data.get('sgc_synapse', species='mouse', post_type=celltyp, field='NMDA_gmax')[0]
    exp_epsc_cv = data.get('sgc_synapse', species='mouse', post_type=celltyp, field='EPSC_cv')
    
    assert abs(exp_epsc_cv / epsc_cv - 1) < 0.1
    assert np.allclose((exp_ampa_gmax, exp_nmda_gmax), (ampa_gmax, nmda_gmax))


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
        for i in range(syn.psd.n_psd):
            ampa += syn.psd.get_vector('ampa', 'g', i)
            nmda += syn.psd.get_vector('nmda', 'g', i)
            if nmda[-1] - nmda[-2] > 0.001:
                raise Exception("Did not reach nmda gmax; need longer run.")
        amax = ampa.max()
        ampa_gmax.append(amax)
        nmda_gmax.append(nmda.max())
        epsc_gmax.append((ampa+nmda).max())
        if plot:
            plt.plot(ampa, pen='g')
            plt.plot(nmda, pen='y')
    
    return (np.mean(ampa_gmax), np.mean(nmda_gmax), 
            np.std(epsc_gmax) / np.mean(epsc_gmax))


    
if __name__ == '__main__':
    test_sgc_bushy_psd(plot=True)
    #test_sgc_tstellate_psd(plot=True)