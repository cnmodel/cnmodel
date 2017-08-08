import scipy.stats
import numpy as np

from .population import Population
from .. import cells


class Pyramidal(Population):
    type = 'pyramidal'
    
    def __init__(self, species='mouse', **kwds):  # ***** NOTE Species - no direct data for mouse (uses RAT data)
        # Note that `cf` is the mean value used when selecting SGCs to connect;
        # it is NOT the measured CF of the cell (although it should be close).
        freqs = self._get_cf_array(species)
        fields = [
            ('cf', float),
        ]
        super(Pyramidal, self).__init__(species, len(freqs), fields=fields, **kwds)
        self._cells['cf'] = freqs
    
    
    def create_cell(self, cell_rec):
        """ Return a single new cell to be used in this population. The 
        *cell_rec* argument is the row from self.cells that describes the cell 
        to be created.
        """
        return cells.Pyramidal.create(species=self.species, **self._cell_args)
