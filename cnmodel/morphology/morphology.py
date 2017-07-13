
class Morphology(object):
    """
    Base class for morphological structure
    Provides the ability to read .hoc files for NEURON simulations.
    
    """
    def init():
        pass
        
    
    def hocReader(self, filename=None):
        """
        Access the hoc reader
        
        Parameters
        ----------
        filename : str (default: None)
            The hoc file to read and save as the morphology object
        """
        self.morphology = hoc_reader(hoc=filename)
    
    def swcReader(self, filename=None):
        """
        Access the swc reader (***NOT IMPLEMENTED***)
        
        Parameters
        ----------
        filename : str (default: None)
            The swc file to read and save as the morphology object
        """
        raise ValueError("swcReader not implmented.")
        
