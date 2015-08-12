import numpy as np
from neuron import h

from .psd import PSD


class Exp2PSD(PSD):
    """
    Simple double-exponential PSD.
    """
    def __init__(self, section, terminal):
        PSD.__init__(self, section, terminal)
        self.syn = h.Exp2Syn(0.5, sec=section)
        self.syn.tau1 = 0.1
        self.syn.tau2 = 0.3
        self.syn.e = 0

        terminal.connect(self.syn, weight=0.01)
 
    @property
    def n_psd(self):
        """The number of postsynaptic densities represented by this object.
        """
        return 1

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
