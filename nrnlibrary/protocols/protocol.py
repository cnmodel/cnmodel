from neuron import h
import numpy as np

class Protocol(object):
    """
    Base class providing common tools for running, analyzing, and displaying
    simulations.
    """
    def __init__(self):
        self.reset()

    def reset(self):
        self._vectors = {}
    
    def __setitem__(self, name, variable):
        """
        Record *variable* during the next run.
        """
        vec = h.Vector()
        self._vectors[name] = vec
        vec.record(variable)
        
    def __getitem__(self, name):
        """
        Return a np array for previously recorded data given *name*.
        """
        return np.array(self._vectors[name])


