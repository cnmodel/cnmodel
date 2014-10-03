import scipy.stats
import numpy as np

from .population import Population
from .. import cells


class Bushy(Population):
    def __init__(self, species='mouse'):
        # Completely fabricated cell distribution: uniform from 4kHz to 90kHz.
        # Note that `cf` is the mean value used when selecting SGCs to connect;
        # it is NOT the measured CF of the cell (although it should be close).
        size = 5000
        fields = [
            ('cf', float),
            ('sgc_sr', int),   # preferred SR group for SGC inputs
        ]
        super(Bushy, self).__init__(species, size, fields=fields)
        self._cells['cf'] = 4000 * 2**np.linspace(0, 4.5, size)
        self._cells['sgc_sr'] = np.arange(size) % 3
    
    def create_cell(self, cell_rec):
        """ Return a single new cell to be used in this population. The 
        *cell_rec* argument is the row from self.cells that describes the cell 
        to be created.
        """
        return cells.Bushy.create(species=self.species)
        
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
            size = np.random.randint(2, 6)
            
            # only select inputs from a single SR group
            sr_vals = pop.cells['sr']
            sr_dist = (sr_vals == cell_rec['sgc_sr']).astype(float)
            
            dist = {'cf': scipy.stats.norm(loc=cf, scale=cf / 100.),
                    'sr': sr_dist}
            
        elif isinstance(pop, populations.DStellate):
            size = np.random.randint(4, 10)
            dist = {'cf': scipy.stats.norm(loc=cf, scale=cf / 10.)}
            
        elif isinstance(pop, populations.TStellate):
            has_input = np.random.randint(10) == 0
            size = 0 if not has_input else np.random.randint(1, 4)
            dist = {'cf': scipy.stats.norm(loc=cf, scale=cf / 50.)}
            
        else:
            raise TypeError("Cannot connect population %s to %s" % (pop, self))

        return size, dist
