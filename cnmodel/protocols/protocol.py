from neuron import h
import numpy as np
from ..util import random, custom_init

class Protocol(object):
    """
    Base class providing common tools for running, analyzing, and displaying
    simulations.
    """
    def __init__(self):
        self.reset()

    def reset(self):
        self._vectors = {}
        
    def run(self, seed=None):
        """
        Run this protocol. 
        
        Subclasses should extend this method.
        """
        if seed is not None:
            random.set_seed(seed)
        self.reset()
    
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

    def custom_init(self, vinit=-60.):
        return custom_init(vinit)
