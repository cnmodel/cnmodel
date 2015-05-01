import numpy as np
import pyqtgraph as pg
import neuron

import nrnlibrary
import nrnlibrary.cells as cells
from nrnlibrary.protocols import SynapseTest
from nrnlibrary.util import random, reset
from nrnlibrary import data


def test_sgc_bushy_psd():
    """Check that sgc-bushy synapses have correct AMPA and NMDA peak conductances.
    """
    prot = SynapseTest()
    reset()
    
    sgc = cells.SGC.create()
    bc = cells.Bushy.create(ttx=True)
    
    # Connect 20 synapses and stimulate once each
    prot.run(sgc.soma, bc.soma, n_synapses=20, temp=33.0, dt=0.025, vclamp=40.0,
             tstop=4.0, stim_params={'NP': 1, 'delay': 0.1})
    
    # For each synapse:
    #    * Add up ampa and nmda conductances across all sites
    #    * Keep track of the maximum ampz and nmda conductance
    #global plt
    #plt = pg.plot()
    ampa_gmax = []
    nmda_gmax = []
    for syn in prot.synapses:
        ampa = np.zeros_like(syn.psd.get_vector('ampa', 'g'))
        nmda = ampa.copy()
        for i in range(syn.psd.n_psd):
            ampa += syn.psd.get_vector('ampa', 'g', i)
            nmda += syn.psd.get_vector('nmda', 'g', i)
            if nmda[-1] - nmda[-2] > 0.001:
                print nmda
                raise Exception("Did not reach nmda gmax; need longer run.")
        ampa_gmax.append(ampa.max())
        nmda_gmax.append(nmda.max())
        #plt.plot(ampa, pen='g')
        #plt.plot(nmda, pen='y')

    exp_ampa_gmax = data.get('sgc_synapse', species='mouse', post_type='bushy', field='AMPA_gmax')[0]
    exp_nmda_gmax = data.get('sgc_synapse', species='mouse', post_type='bushy', field='NMDA_gmax')[0] 
    assert np.allclose((exp_ampa_gmax, exp_nmda_gmax), (np.mean(ampa_gmax), np.mean(nmda_gmax)))

    
def test_sgc_tstellate_psd():
    """Check that sgc-bushy synapses have correct AMPA and NMDA peak conductances.
    """
    random.set_seed(34783)
    reset()
    sgc = cells.SGC.create()
    bc = cells.Bushy.create(ttx=True)
    tsc = cells.TStellate.create(ttx=True)
    
    bsyn = sgc.connect(bc)
    tsyn = sgc.connect(tsc)
    
    # Force terminals to release when the presynaptic cell fires
    bsyn.terminal.relsite.F = 1.0
    bsyn.terminal.relsite.Dep_Flag = 0
    tsyn.terminal.relsite.F = 1.0
    tsyn.terminal.relsite.Dep_Flag = 0
    
    neuron.h.celsius = 37.0
    neuron.h.dt = 0.025
    neuron.h.finitialize()
    
    
    N = 1000
    isyn = np.empty((N, 2, 2))  # (time, bushy/stell, ampa/nmda)
    t = np.empty(N)
    for i in range(N):
        # voltage clamp postsynaptic cells at +40
        bc.soma.v = 40
        tsc.soma.v = 40
        
        # Generate a spike
        if i == 100:
            sgc.soma.v = 100
        #if 100< i < 104:
            #for j in range(bsyn.terminal.n_rzones):
                #bsyn.terminal.relsite.XMTR[j] = 3
            #tsyn.terminal.relsite.XMTR[0] = 3
        #else:
            #for j in range(bsyn.terminal.n_rzones):
                #bsyn.terminal.relsite.XMTR[j] = 0
            #tsyn.terminal.relsite.XMTR[0] = 0
        
        neuron.h.fadvance()
        
        isyn[i, 0, 0] = bsyn.psd.ampa_psd[0].i
        isyn[i, 0, 1] = bsyn.psd.nmda_psd[0].i
        isyn[i, 1, 0] = tsyn.psd.ampa_psd[0].i
        isyn[i, 1, 1] = tsyn.psd.nmda_psd[0].i
        t[i] = neuron.h.t
        
    import pyqtgraph as pg
    plt = pg.plot()
    plt.plot(t, isyn[:, 0, 0], pen='r')
    plt.plot(t, isyn[:, 0, 1], pen='y')
    plt.plot(t, isyn[:, 1, 0], pen='b')
    plt.plot(t, isyn[:, 1, 1], pen='c')
    
    # TODO: replace with correct values
    bushy_ampa_max = 12345.678
    bushy_nmda_max = 12345.678
    ts_ampa_max = 12345.678
    ts_nmda_max = 12345.678
    assert np.close(isyn.max(axis=0), np.array([[bushy_ampa_max, bushy_nmda_max], 
                                                [ts_ampa_max, ts_nmda_max]]))
    
    
if __name__ == '__main__':
    test_sgc_bushy_psd()