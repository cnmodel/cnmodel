import numpy as np
from neuron import h

from .psd import PSD


class Exp2PSD(PSD):
    """
    Simple double-exponential PSD from Neuron (fast).
    """
    def __init__(self, section, terminal, loc=0.5):
        """
        Parameters
        ----------
        section : Section
            The postsynaptic section in which to insert the receptor mechanism.
        terminal : Terminal
            The presynaptic Terminal instance
        loc : float, default=0.5
            Position on the postsynaptic section to insert the mechanism, from [0..1]. 
        
        """
        PSD.__init__(self, section, terminal)
        self.syn = h.Exp2Syn(loc, sec=section)
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
        site.
        
        
        Parameters
        ----------
        \*args : 
            Allowed parameters are 'i' (current), 'g' (conductnace), and 'Open' (open probability).

        """
        self.vectors = {'ampa': [], 'nmda': []}
        for receptor in self.vectors:
            for mech in getattr(self, receptor+'_psd'):
                vec = {}
                for var in args:
                    vec[var] = h.Vector()
                    vec[var].record(getattr(mech, '_ref_'+var))
                self.vectors[receptor].append(vec)
        
    def get_vector(self, var):
        """Return an array from a previously recorded vector. 
        
        Parameters
        ----------
        receptor : str
             May be 'ampa' or 'nmda'
        var : str
            Allowed parameters are 'i' (current), 'g' (conductance), and 'Open' (open probability).
        i : int, default=0
             The integer index of the psd (if this is a multi-site synapse)
        
        """
        v = self.vectors[receptor][i][var]
        return np.array(v)
