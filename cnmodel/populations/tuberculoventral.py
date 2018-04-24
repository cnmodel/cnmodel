import scipy.stats
import numpy as np

from .population import Population
from .. import cells

class Tuberculoventral(Population):
    type = 'tuberculoventral'
    
    def __init__(self, species='mouse', **kwds):
        # Note that `cf` is the mean value used when selecting SGCs to connect;
        # it is NOT the measured CF of the cell (although it should be close).
        freqs = self._get_cf_array(species)
        fields = [
            ('cf', float),
        ]
        super(Tuberculoventral, self).__init__(species, len(freqs), fields=fields, **kwds)
        self._cells['cf'] = freqs
    
    def create_cell(self, cell_rec):
        """ Return a single new cell to be used in this population. The 
        *cell_rec* argument is the row from self.cells that describes the cell 
        to be created.
        """
        return cells.Tuberculoventral.create(species=self.species, **self._cell_args)

    def connection_stats(self, pop, cell_rec):
        """ The population *pop* is being connected to the cell described in 
        *cell_rec*. Return the number of presynaptic cells that should be
        connected and a dictionary of distributions used to select cells 
        from *pop*. 
        """
        size, dist = Population.connection_stats(self, pop, cell_rec)
        
        from .. import populations

        if isinstance(pop, populations.SGC):
            # only select SGC inputs from low- and medium SR. See:
            #     Spectral Integration by Type II Interneurons in Dorsal Cochlear Nucleus
            #     George A. Spirou, Kevin A. Davis, Israel Nelken, Eric D. Young
            #     Journal of Neurophysiology Aug 1999, 82 (2) 648-663;
            dist['sr'] = (pop.cells['sr'] < 2).astype(float)

        return size, dist
