import numpy as np
from neuron import h

from .psd import PSD


class GluPSD(PSD):
    """
    Glutamatergic PSD with ionotropic AMPA / NMDA receptors

    This creates a set of postsynaptoc NMDA and AMPA receptors, one pair
    per terminal release site. Receptors are connected to the XMTR range
    variable of the terminal release mechanisms.
    
    Parameters
    ==========
    section : Section instance
        The postsynaptic section into which the receptor mechanisms should be
        attached
    terminal : Terminal instance
        The presynaptic terminal that provides input to the receptor XMTR
        variables. 
    ampa_gmax : float
        Maximum conductance of AMPARs
    nmda_gmax : float
        Maximum conductance of NMDARs
    gvar : float
        Coefficient of variation for randomly adjusting ampa_gmax and nmda_gmax.
        Note that ampa and nmda maximum conductances will be scaled together,
        but the scale values will be selected randomly for each pair of 
        receptor mechanisms.
    eRev : float
        Reversal potential to use for both receptor types.
        
    Notes
    =====
    
    *ampa_gmax* and *nmda_gmax* should be provided as the maximum *measured*
    conductances; these will be automatically corrected for the maximum open
    probability of the receptor mechanisms.
    
    GluPSD does not include a cleft mechanism because AMPATRUSSELL implements
    its own cleft and NMDA_Kampa is slow enough that a cleft would have 
    insignificant effect.
    """
    def __init__(self, section, terminal,
                 ampa_gmax, nmda_gmax
                 gvar=0, eRev=0,
                 identifier=0):
        PSD.__init__(self, section, terminal)
        
        self.pre_sec = terminal.section
        self.post_sec = section
        
        from .. import cells
        
        self.pre_cell = cells.cell_from_section(self.pre_sec)
        self.post_cell = cells.cell_from_section(self.post_sec)

        relzone = terminal.relsite
        n_rzones = terminal.n_rzones
        
        # and then make a set of postsynaptic receptor mechanisms
        (ampa_psd, nmda_psd, par, parn) = self.template_iGluR_PSD()
        
        for k in range(0, n_rzones):
            # Connect terminal to psd
            relzone.setpointer(relzone._ref_XMTR[k], 'XMTR', ampa_psd[k])
            relzone.setpointer(relzone._ref_XMTR[k], 'XMTR', nmda_psd[k])
            
            # add a little variability - gvar is CV of amplitudes
            v = 1.0 + gvar * np.random.standard_normal()
            
            # set gmax and eRev for each postsynaptic receptor mechanism
            ampa_psd[k].gmax = ampa_gmax * v 
            ampa_psd[k].Erev = eRev
            nmda_psd[k].gmax = nmda_gmax * v
            nmda_psd[k].Erev = eRev
            nmda_psd[k].vshift = 0
        
        par = list(par)
        par.extend(parn)
                
        self.ampa_psd = ampa_psd
        self.nmda_psd = nmda_psd
        self.all_psd = nmda_psd + ampa_psd
        self.par = par

    def template_iGluR_PSD(self, cellname=None):
        """
        Create an ionotropic Glutamate receptor "PSD"
        Each PSD has receptors for each active zone, which must be matched (connected) to presynaptic
        terminals. Each PSD recetpor consists of an AMPATRUSSELL and an NMDA_KAMPA receptor
        Inputs:
            sec: The template requires a segment to insert the receptors into
            nReceptors: The number of receptor sites to insert
            cellname: Bushy/MNTB/stellate: determines ampa receptor kinetics
        Outputs:
            (psd, psdn, par, parn)
            psd is the list of PSDs that were created (AMPA)
            psdn is the list of NMDA PSDs (same number as psd, just the NMDARs)
            par: dictionary of AMPAR kinetics as inserted
            parn: NMDA ratio
        Side Effecdts: None
        """
        sec = self.post_sec
        nReceptors = self.terminal.n_rzones
        
        psd = []
        psdn = []
        sec.push()
        for k in range(0, nReceptors):
            psd.append(h.AMPATRUSSELL(0.5, sec)) # raman/trussell AMPA with rectification
            psdn.append(h.NMDA_Kampa(0.5, sec)) # Kampa state model NMDA receptors

            if cellname in ['bushy', 'MNTB']:
                psd[-1].Ro1 = 107.85
                psd[-1].Ro2 = 0.6193
                psd[-1].Rc1 = 3.678
                psd[-1].Rc2 = 0.3212
            if cellname == 'stellate':
                psd[-1].Ro1 = 39.25
                psd[-1].Ro2 = 4.40
                psd[-1].Rc1 = 0.667
                psd[-1].Rc2 = 0.237
                psd[-1].PA = 0.1

        h.pop_section()
        par = {'Ro1': ('r', psd[0].Ro1),
            'Ro2': ('r', psd[0].Ro2),
            'Rc1': ('r', psd[0].Rc1),
            'Rc2': ('r', psd[0].Rc2), }
        
        return (psd, psdn, par, {})
