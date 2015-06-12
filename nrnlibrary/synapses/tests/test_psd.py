import numpy as np
import pyqtgraph as pg
import neuron

import nrnlibrary
import nrnlibrary.cells as cells
from nrnlibrary.protocols import SynapseTest
from nrnlibrary.util import random, reset
from nrnlibrary import data


def test_sgc_bushy_psd(plot=False):
    """Check that sgc-bushy synapses have correct AMPA and NMDA peak conductances.
    """
    random.set_seed(23572385)
    reset()
    bc = cells.Bushy.create(ttx=True)
        
    (ampa_gmax, nmda_gmax, epsc_cv) = measure_gmax(bc, n_syn=20, tstop=4.0, plot=plot)

    exp_ampa_gmax = data.get('sgc_synapse', species='mouse', post_type='bushy', field='AMPA_gmax')[0]
    exp_nmda_gmax = data.get('sgc_synapse', species='mouse', post_type='bushy', field='NMDA_gmax')[0] 
    exp_epsc_cv = data.get('sgc_synapse', species='mouse', post_type='bushy', field='EPSC_cv')
    
    assert np.allclose((exp_ampa_gmax, exp_nmda_gmax), (ampa_gmax, nmda_gmax))
    assert abs(exp_epsc_cv / epsc_cv - 1) < 0.1

    
def test_sgc_tstellate_psd(plot=False):
    """Check that sgc-bushy synapses have correct AMPA and NMDA peak conductances.
    """
    random.set_seed(34754398)
    reset()
    tsc = cells.TStellate.create(ttx=True)
    (ampa_gmax, nmda_gmax, epsc_cv) = measure_gmax(tsc, n_syn=20, tstop=5.0, plot=plot)
    
    exp_ampa_gmax = data.get('sgc_synapse', species='mouse', post_type='tstellate', field='AMPA_gmax')[0]
    exp_nmda_gmax = data.get('sgc_synapse', species='mouse', post_type='tstellate', field='NMDA_gmax')[0]
    exp_epsc_cv = data.get('sgc_synapse', species='mouse', post_type='tstellate', field='EPSC_cv')
    
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