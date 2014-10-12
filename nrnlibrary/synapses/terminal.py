class Terminal(object):
    """
    Base class for axon terminals. A terminal has a single postsynaptic 
    neuron, but may have multiple release zones. It defines a release mechanism
    with a NetCon input (triggering from presynaptic voltage or calcium level)
    and either NetCon or pointer output for driving a PSD.
    
    """
    def __init__(self, section):
        self._section = section
        
    @property
    def section(self):
        """ The cell section this terminal is attached to.
        """
        return self._section
    
    @property
    def cell(self):
        """ The cell this terminal is attached to.
        """
        from ..cells import Cell
        return Cell.from_section(self.section)
    