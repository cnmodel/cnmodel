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
    ampa_params : dict
        Dictionary containing kinetic parameters for AMPA mechanism. Suggested
        keys are Ro1, Ro2, Rc1, Rc2, and PA.
        
    Notes
    =====
    
    *ampa_gmax* and *nmda_gmax* should be provided as the maximum *measured*
    conductances; these will be automatically corrected for the maximum open
    probability of the receptor mechanisms.
    
    GluPSD does not include a cleft mechanism because AMPATRUSSELL implements
    its own cleft and NMDA_Kampa is slow enough that a cleft would have 
    insignificant effect.
    """
    def __init__(self, section, terminal, ampa_gmax, nmda_gmax,
                 gvar=0, eRev=0, ampa_params=None):
        PSD.__init__(self, section, terminal)
        
        ampa_params = {} if ampa_params is None else ampa_params
        
        # and then make a set of postsynaptic receptor mechanisms
        ampa_psd = []
        nmda_psd = []
        relsite = terminal.relsite
        self.section.push()
        for i in range(0, terminal.n_rzones):
            # create mechanisms
            ampa = h.AMPATRUSSELL(0.5, self.section) # raman/trussell AMPA with rectification
            nmda = h.NMDA_Kampa(0.5, self.section) # Kampa state model NMDA receptors

            # Connect terminal to psd
            relsite.setpointer(relsite._ref_XMTR[i], 'XMTR', ampa)
            relsite.setpointer(relsite._ref_XMTR[i], 'XMTR', nmda)
            
            # Set any extra ampa parameters provided by the caller
            # (Ro1, Ro2, Rc1, Rc2, PA, ...)
            for k,v in ampa_params.items():
                setattr(ampa, k, v)
            
            # add a little variability - gvar is CV of amplitudes
            v = 1.0 + gvar * np.random.standard_normal()
            
            # set gmax and eRev for each postsynaptic receptor mechanism
            ampa.gmax = ampa_gmax * v
            ampa.Erev = eRev
            nmda.gmax = nmda_gmax * v
            nmda.Erev = eRev
            nmda.vshift = 0
            
            ampa_psd.append(ampa)
            nmda_psd.append(nmda)
        
        h.pop_section()
        
        self.ampa_psd = ampa_psd
        self.nmda_psd = nmda_psd
        self.all_psd = nmda_psd + ampa_psd

    @property
    def n_psd(self):
        """The number of postsynaptic densities represented by this object.
        """
        return len(self.ampa_psd)

    def record(self, *args):
        """Create a new set of vectors to record parameters for each release
        site. Allowed parameters are 'i', 'g', and 'Open'.
        """
        self.vectors = {'ampa': [], 'nmda': []}
        for receptor in self.vectors:
            for mech in getattr(self, receptor+'_psd'):
                vec = {}
                for var in args:
                    vec[var] = h.Vector()
                    vec[var].record(getattr(mech, '_ref_'+var))
                self.vectors[receptor].append(vec)
        
    def get_vector(self, receptor, var, i=0):
        """Return an array from a previously recorded vector. 
        
        *receptor* may be 'ampa' or 'nmda'
        *var* may be 'i', 'g', or 'Open'
        *i* is the integer index of the psd (if this is a multi-site synapse)
        """
        v = self.vectors[receptor][i][var]
        return np.array(v)
