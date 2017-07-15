import scipy.stats
import numpy as np

from .population import Population
from .. import cells


class Pyramidal(Population):
    type = 'pyramidal'
    
    def __init__(self, species='rat', **kwds):  # ***** NOTE SPecies - no diret data for mouse
        # Completely fabricated cell distribution: uniform from 2kHz to 64kHz.
        # Note that `cf` is the mean value used when selecting SGCs to connect;
        # it is NOT the measured CF of the cell (although it should be close).
        size = 3000
        fields = [
            ('cf', float),
        ]
        super(Pyramidal, self).__init__(species, size, fields=fields, **kwds)
        self._cells['cf'] = 2000 * 2**np.linspace(0, 5.0, size)
    
    def create_cell(self, cell_rec):
        """ Return a single new cell to be used in this population. The 
        *cell_rec* argument is the row from self.cells that describes the cell 
        to be created.
        """
        return cells.Pyramidal.create(species=self.species, **self._cell_args)

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
            size = np.random.randint(4, 10)
            dist = {'cf': scipy.stats.norm(loc=cf, scale=cf / 50.)}
            
        elif isinstance(pop, populations.DStellate):
            size = np.random.randint(10, 20)
            dist = {'cf': scipy.stats.norm(loc=cf, scale=cf / 5.)}

        elif isinstance(pop, populations.Tuberculoventral):
            size = np.random.randint(10, 20)
            dist = {'cf': scipy.stats.norm(loc=cf, scale=cf / 50.)}  # narrow population

        # elif isinstance(pop, populations.TStellate):  # can be added once we know TStellate projections to DCN
        #     size = np.random.randint(10, 20)
        #     dist = {'cf': scipy.stats.norm(loc=cf, scale=cf / 50.)}  # narrow population
            
        else:
            raise TypeError("Cannot connect population %s to %s" % (pop, self))

        return size, dist
