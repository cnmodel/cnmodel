import scipy.stats
import numpy as np

from .population import Population
from .. import cells


class DStellate(Population):
    def __init__(self, species='mouse'):
        # Completely fabricated cell distribution: uniform from 4kHz to 90kHz.
        # Note that `cf` is the mean value used when selecting SGCs to connect;
        # it is NOT the measured CF of the cell (although it should be close).
        size = 500
        fields = [
            ('cf', float),
        ]
        super(DStellate, self).__init__(species, size, fields=fields)
        self._cells['cf'] = 4000 * 2**np.linspace(0, 4.5, size)
    
    def create_cell(self, cell_rec):
        """ Return a single new cell to be used in this population. The 
        *cell_rec* argument is the row from self.cells that describes the cell 
        to be created.
        """
        return cells.DStellate.create(species=self.species)

    def connection_stats(self, pop, cell_rec):
        """ The population *pop* is being connected to the cell described in 
        *cell_rec*. Return the number of presynaptic cells that should be
        connected and a dictionary of distributions used to select cells 
        from *pop*. 
        """
        from .. import populations
        
        # Fabricated convergence distributions (how many presynaptic 
        # cells to connect) and input distributions (which presynaptic
        # cells to connect). 
        cf = cell_rec['cf']
        if isinstance(pop, populations.SGC):
            size = np.random.randint(10, 20)
            dist = {'cf': scipy.stats.norm(loc=cf, scale=cf / 10.)}
            
        else:
            raise TypeError("Cannot connect population %s to %s" % (pop, self))

        return size, dist

