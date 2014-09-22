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
        ]
        super(Bushy, self).__init__(species, size, fields=fields)
        self._cells['cf'] = 4000 * 2**np.linspace(0, 4.5, size)
    
    def create_cells(self, cells):
        """ Instantiate each cell in *cells*, which is a list of indexes into
        self.cells.
        """
        for i in cells:
            self._cells[i]['cell'] = cells.Bushy.create(species=self.species)
        
    def connect_pop_to_cell(self, pop, cell_index):
        """ Connect cells in a presynaptic population to the cell in this 
        population at *cell_index*. Return the presynaptic indexes of cells
        that were connected.
        """
        cell = self._cells[cell_index]['cell']
        
        if isinstance(pop, populations.SGC):
            # Fabricated input distribution (this should be log-normal, or 
            # perhaps cf should be expressed in octaves?)
            cf = self._cells[cell_index]['cf']
            dist = scipy.stats.norm(loc=cf, scale=cf / 10.)
            
            # Fabricated convergence distribution
            size = np.random.randint(3) + 2
            
            # Select SGCs from distribution, create, and connect to this cell
            pre_cells = pop.select(size=size, create=False, cf=dist)
            for j in pre_cells:
                sgc = pop.get_cell(j)
                sgc.connect(cell)
                
        else:
            raise TypeError("Cannot connect population %s to %s" % (pop, self))
