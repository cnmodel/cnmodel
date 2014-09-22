import numpy as np

from .population import Population
from .. import cells

class SGC(Population):
    def __init__(self, species='mouse'):
        # Completely fabricated cell distribution: uniform from 4kHz to 90kHz.
        size = 10000
        super(SGC, self).__init__(species, size, fields=[('cf', float)])
        self._cells['cf'] = 4000 * 2**np.linspace(0, 4.5, size)
    
    def create_cells(self, cells):
        """ Instantiate each cell in *cells*, which is a list of indexes into
        self.cells.
        """
        for i in cells:
            self._cells[i]['cell'] = cells.SGC.create(species=self.species)
        
    def connect_pop_to_cell(self, pop, index):
        # SGC does not support any inputs
        assert len(self.connections) == 0
