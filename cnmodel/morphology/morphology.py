
class Morphology(object):
    """
    Base class for morphological structure
    Provides the ability to read .hoc files for NEURON simulations.
    
    """
    def init():
        pass
        
    
    def hocReader(self, filename=None):
        self.morphology = hoc_reader(hoc=filename)
    
    def swcReader(self, filename=None):
        raise ValueError("swcReader not implmented.")
        
    
