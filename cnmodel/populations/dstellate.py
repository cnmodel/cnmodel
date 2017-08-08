import scipy.stats
import numpy as np

from .population import Population
from .. import cells


class DStellate(Population):
    type = 'dstellate'
    
    def __init__(self, species='mouse', **kwds):
        # Note that `cf` is the mean value used when selecting SGCs to connect;
        # it is NOT the measured CF of the cell (although it should be close).
        freqs = self._get_cf_array(species)
        fields = [
            ('cf', float),
        ]
        super(DStellate, self).__init__(species, len(freqs), fields=fields, **kwds)
        self._cells['cf'] = freqs
    
    def create_cell(self, cell_rec):
        """ Return a single new cell to be used in this population. The 
        *cell_rec* argument is the row from self.cells that describes the cell 
        to be created.
        """
        return cells.DStellate.create(species=self.species, **self._cell_args)
