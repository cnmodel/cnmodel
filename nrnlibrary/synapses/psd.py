

class PSD(object):
    """
    Base class for postsynaptic density mechanisms, possibly including cleft. 
    May accept either NetCon or pointer inputs from a Terminal, and directly
    modifies the membrane potential and/or ion concentrations of the 
    postsynaptic cell.    
    """
    def __init__(self, section, terminal):
        self._section = section
        self._terminal = terminal
        
    @property
    def section(self):
        """ The cell section this PSD is attached to.
        """
        return self._section
    
    @property
    def cell(self):
        """ The cell this PSD is attached to.
        """
        from ..cells import Cell
        return Cell.from_section(self.section)

    @property
    def terminal(self):
        """ The presynaptic terminal connected to this PSD.
        """
        return self._terminal
