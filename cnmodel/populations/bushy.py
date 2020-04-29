import scipy.stats
import numpy as np

from .population import Population
from .. import cells


class Bushy(Population):
    """Population of bushy cells.
    
    Cells are distributed uniformly from 2kHz to 64kHz.
    
    Note that `cf` is the mean value used when selecting SGCs to connect;
    it is NOT the measured CF of the cell (although it should be close).
    """
    type = 'bushy'
    
    def __init__(self, species='mouse', **kwds):
        freqs = self._get_cf_array(species)
        fields = [
            ('cf', float),
            ('input_sr', list),   # distribution probability of SGC SR groups
            ('sr', int),
        ]
        super(Bushy, self).__init__(species, len(freqs), fields=fields, **kwds)
        self._cells['cf'] = freqs
        self._cells['input_sr'] = [np.tile([1., 1., 1.], len(freqs))]

    def create_cell(self, cell_rec):
        """ Return a single new cell to be used in this population. The 
        *cell_rec* argument is the row from self.cells that describes the cell 
        to be created.
        """
        return cells.Bushy.create(species=self.species, **self._cell_args)
        
    def connection_stats(self, pop, cell_rec):
        """ The population *pop* is being connected to the cell described in 
        *cell_rec*. Return the number of presynaptic cells that should be
        connected and a dictionary of distributions used to select cells 
        from *pop*. 
        """
        size, dist = Population.connection_stats(self, pop, cell_rec)
        
        from .. import populations

        if isinstance(pop, populations.SGC):
            # only select SGC inputs from a single SR group
            # (this relationship is hypothesized based on reconstructions of
            # endbulbs)
            sr_vals = pop.cells['sr']
            u = np.random.choice(sr_vals)  # assign input sr for this cell
#            print('u: ', u)
            # pick from one sr group for all inputs, with prob same as distribution in nerve
            dist['sr'] = (sr_vals == u).astype(float)
            self._cells['sr'] = u

        return size, dist
