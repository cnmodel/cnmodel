from neuron import h
import numpy as np
from ..util import random

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

    # moved to util
    # def custom_init(self, vinit=-60.):
    #     """
    #     perform a custom initialization of the current section to vinit
    #     """
    #     initdur = 1e6
    #     tdt = h.dt
    #     dtstep = 1e3
    #     h.finitialize()
    #     h.t = -initdur
    #     tmp = h.cvode.active()
    #     if tmp != 0:
    #         h.cvode.active(0)
    #     h.dt = dtstep
    #     while h.t < 0:
    #         h.fadvance()
    #     if tmp != 0:
    #         h.cvode.active(1)
    #     h.t = 0
    #     if h.cvode.active():
    #         h.cvode.re_init()
    #     else:
    #         h.fcurrent()
    #     h.frecord_init()
    #     h.dt = tdt
    #     h.fcurrent()