import scipy.stats
import numpy as np

from .population import Population
from .. import cells


class TStellate(Population):
    type = 'tstellate'
    
    def __init__(self, species='mouse', **kwds):
        # Completely fabricated cell distribution: uniform from 4kHz to 90kHz.
        # Note that `cf` is the mean value used when selecting SGCs to connect;
        # it is NOT the measured CF of the cell (although it should be close).
        size = 5000
        fields = [
            ('cf', float),
        ]
        super(TStellate, self).__init__(species, size, fields=fields, **kwds)
        self._cells['cf'] = 4000 * 2**np.linspace(0, 4.5, size)
    
    def create_cell(self, cell_rec):
        """ Return a single new cell to be used in this population. The 
        *cell_rec* argument is the row from self.cells that describes the cell 
        to be created.
        """
        return cells.TStellate.create(species=self.species, **self._cell_args)
