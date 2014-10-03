import numpy as np

from .population import Population
from .. import cells

class SGC(Population):
    def __init__(self, species='mouse'):
        # Completely fabricated cell distribution: uniform from 4kHz to 90kHz.
        # Evenly divided between SR groups
        size = 10000
        fields = [
            ('cf', float),
            ('sr', int),  # 0=low sr, 1=mid sr, 2=high sr
        ]
        super(SGC, self).__init__(species, size, fields=fields)
        self._cells['cf'] = 4000 * 2**np.linspace(0, 4.5, size)
        self._cells['sr'] = np.arange(size) % 3
    
    def create_cell(self, cell_rec):
        """ Return a single new cell to be used in this population. The 
        *cell_rec* argument is the row from self.cells that describes the cell 
        to be created.
        """
        return cells.SGC.create(species=self.species)
        
    def connect_pop_to_cell(self, pop, index):
        # SGC does not support any inputs
        assert len(self.connections) == 0
