import scipy.stats
import numpy as np

from .population import Population
from .. import cells


class DStellate(Population):
    type = 'dstellate'
    
    def __init__(self, species='mouse', **kwds):
        # Completely fabricated cell distribution: uniform from 2kHz to 64kHz.
        # Note that `cf` is the mean value used when selecting SGCs to connect;
        # it is NOT the measured CF of the cell (although it should be close).
        size = 500
        fields = [
            ('cf', float),
        ]
        super(DStellate, self).__init__(species, size, fields=fields, **kwds)
        self._cells['cf'] = 2000 * 2**np.linspace(0, 5.0, size)
    
    def create_cell(self, cell_rec):
        """ Return a single new cell to be used in this population. The 
        *cell_rec* argument is the row from self.cells that describes the cell 
        to be created.
        """
        return cells.DStellate.create(species=self.species, **self._cell_args)
