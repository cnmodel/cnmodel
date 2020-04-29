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
            ('input_sr', list),   # distribution probability of SGC SR groups
        ]
        super(DStellate, self).__init__(species, len(freqs), fields=fields, **kwds)
        self._cells['cf'] = freqs
        self._cells['input_sr'] = [np.tile([1., 1., 1.], len(freqs))]
    
    def create_cell(self, cell_rec):
        """ Return a single new cell to be used in this population. The 
        *cell_rec* argument is the row from self.cells that describes the cell 
        to be created.
        """
        return cells.DStellate.create(species=self.species, **self._cell_args)
    
    def connection_stats(self, pop, cell_rec):
        """ The population *pop* is being connected to the cell described in 
        *cell_rec*. Return the number of presynaptic cells that should be
        connected and a dictionary of distributions used to select cells 
        from *pop*. 
        """
        size, dist = Population.connection_stats(self, pop, cell_rec)
        
        from .. import populations

        if isinstance(pop, populations.SGC):
            dist['sr'] = (pop.cells['sr'] < 2).astype(float)
        return size, dist