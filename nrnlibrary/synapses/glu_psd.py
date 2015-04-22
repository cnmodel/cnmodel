import numpy as np
from neuron import h

from .psd import PSD


class GluPSD(PSD):
    """
    Glutamatergic PSD with ionotropic AMPA / NMDA receptors
    """
    def __init__(self, section, terminal,
                 ampa_gmax,
                 nmda_ampa_ratio,
                 message=None, debug=False,
                 gvar=0, eRev=0,
                 nmda_ratio=1.0, identifier=0):
        """ This routine generates the synaptic connections from one presynaptic
            input onto a postsynaptic cell.
            Each connection is a stochastic presynaptic synapse ("presynaptic") with
            depression and facilitation.
            Each synapse can have  multiple releasesites "releasesites" on the
            target cell, as set by "NRZones".
            Each release site releases transmitter using "cleftXmtr"
            Each release site's transmitter is in turn attached to a PSD at each ending ("psd")
            Each psd can have a different conductance centered about the mean of
            gmax, according to a gaussian distribution set by gvar.
            
        Notes:
        
        *ampa_gmax* should be provided as the maximum *measured* AMPA conductance;
        this will be automatically corrected for the maximum open probability of
        the AMPA mechanism.
        
        *nmda_ampa_ratio* should be the ratio nmda/ampa Po measured at +40 mV.
        """
        PSD.__init__(self, section, terminal)
        
        pre_sec = terminal.section
        post_sec = section
        
        from .. import cells
        self.AN_Po_Ratio = 23.2917 # ratio of open probabilities for AMPA and NMDAR's at peak currents
        self.AMPA_Max_Po = 0.44727
        self.NMDARatio = 0.0
        
        self.pre_cell = cells.cell_from_section(pre_sec)
        self.post_cell = cells.cell_from_section(post_sec)

        # get AMPA gmax corrected for max open probability
        gmax = ampa_gmax / self.AMPA_Max_Po
        
        relzone = terminal.relsite
        n_rzones = terminal.n_rzones
        
        #
        # Create cleft mechanisms
        # 
        clefts = []
        for k in range(0, n_rzones):
            cl = h.cleftXmtr(0.5, sec=post_sec)
            clefts.append(cl) # cleft
        
        # and then make a set of postsynaptic receptor mechanisms
        #        print 'PSDTYPE: ', psdtype
        (ampa_psd, nmda_psd, par, parn) = template_iGluR_PSD(post_sec, nReceptors=n_rzones,
                                                    nmda_ratio=nmda_ratio)
        
        # Connect terminal to psd (or cleft)
        for k in range(0, n_rzones):
            # Note: cleft kinetics is implemented in the AMPA mechanism
            relzone.setpointer(relzone._ref_XMTR[k], 'XMTR', ampa_psd[k])
            
            # Note: NMDA has no cleft mechanism, but it has a slow response that
            # would not be strongly affected by the relatively fast cleft kinetics.
            relzone.setpointer(relzone._ref_XMTR[k], 'XMTR', nmda_psd[k]) # include NMDAR's as well at same release site
            
            v = 1.0 + gvar * np.random.standard_normal()
            ampa_psd[k].gmax = gmax * v # add a little variability - gvar is CV of amplitudes
            ampa_psd[k].Erev = eRev # set the reversal potential
            
            # also adjust the nmda receptors at the same synapse
            gNAR = nmda_ampa_ratio * self.AN_Po_Ratio * self.NMDARatio
            nmda_psd[k].gmax = gmax * v * gNAR
            nmda_psd[k].Erev = eRev
            nmda_psd[k].vshift = 0
        
        par = list(par)
        par.extend(parn)
        if message is not None:
            print message
                
        self.ampa_psd = ampa_psd
        self.nmda_psd = nmda_psd
        self.all_psd = nmda_psd + ampa_psd
        self.clefts = clefts
        self.par = par


def template_iGluR_PSD(sec, nReceptors=1, debug=False, cellname=None, message=None, nmda_ratio=1):
    """
    Create an ionotropic Glutamate receptor "PSD"
    Each PSD has receptors for each active zone, which must be matched (connected) to presynaptic
    terminals. Each PSD recetpor consists of an AMPATRUSSELL and an NMDA_KAMPA receptor
    Inputs:
        sec: The template requires a segment to insert the receptors into
        nReceptors: The number of receptor sites to insert
        debug: flag for debugging (prints extra information)
        cellname: Bushy/MNTB/stellate: determines ampa receptor kinetics
        message: Not used.
        nmda_ratio: The relative conductance of the open NMDA receptors to the open AMPA receptors.
    Outputs:
        (psd, psdn, par, parn)
        psd is the list of PSDs that were created (AMPA)
        psdn is the list of NMDA PSDs (same number as psd, just the NMDARs)
        par: dictionary of AMPAR kinetics as inserted
        parn: NMDA ratio
    Side Effecdts: None
    """
    psd = []
    psdn = []
    sec.push()
    AN_Po_Ratio = 23.2917 # ratio of open probabilities for AMPA and NMDAR's at peak currents
    for k in range(0, nReceptors):
        psd.append(h.AMPATRUSSELL(0.5, sec)) # raman/trussell AMPA with rectification
        psdn.append(h.NMDA_Kampa(0.5, sec)) # Kampa state model NMDA receptors

        if cellname in ['bushy', 'MNTB']:
            psd[-1].Ro1 = 107.85
            psd[-1].Ro2 = 0.6193
            psd[-1].Rc1 = 3.678
            psd[-1].Rc2 = 0.3212
            gNAR = 0.036 * AN_Po_Ratio * nmda_ratio # 0.36*AN_Po_Ratio*nmda_ratio
            psdn[-1].gmax = psdn[-1].gmax * gNAR
            #if k == 0:
            #    print "Bushy NMDAR set to %8.2f" % psdn[-1].gNAR
        if cellname == 'stellate':
            psd[-1].Ro1 = 39.25
            psd[-1].Ro2 = 4.40
            psd[-1].Rc1 = 0.667
            psd[-1].Rc2 = 0.237
            psd[-1].PA = 0.1
            gNAR = 1 * AN_Po_Ratio * nmda_ratio
            psdn[-1].gmax = psdn[-1].gmax * gNAR

    h.pop_section()
    par = {'Ro1': ('r', psd[0].Ro1),
           'Ro2': ('r', psd[0].Ro2),
           'Rc1': ('r', psd[0].Rc1),
           'Rc2': ('r', psd[0].Rc2), }
    #parn = {'gNAR': ('4', gNAR), }
    return (psd, psdn, par, {})
