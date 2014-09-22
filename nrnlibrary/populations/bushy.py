import scipy.stats
from .population import Population
from .. import cells


class Bushy(Population):
    def __init__(self, species='mouse'):
        # Completely fabricated cell distribution: uniform from 4kHz to 90kHz.
        # Note that `cf` is the mean value used when selecting SGCs to connect;
        # it is NOT the measured CF of the cell (although it should be close).
        size = 5000
        super(Bushy, self).__init__(species, size, fields=[('cf', float)])
        self._cells['cf'] = 4000 * 2**np.linspace(0, 4.5, size)
    
    def create_cells(self, cells):
        """ Instantiate each cell in *cells*, which is a list of indexes into
        self.cells.
        """
        for i in cells:
            self._cells[i]['cell'] = cells.Bushy.create(species=self.species)
        
    def resolve_inputs(self, depth=1):
        """ For each _real_ cell in the population, select a set of 
        presynaptic partners from each connected population and generate a 
        synapse from each.
        """
        for i in self.unresolved_cells():
            cf = self._cells[i]['cf']
            cell = self._cells[i]['cell']
            
            # select cells from each population to connect to this cell
            for pop in self.connections:
                if isinstance(pop, populations.SGC):
                    # Fabricated input distribution (this should be log-normal, or 
                    # perhaps cf should be expressed in octaves?)
                    dist = scipy.stats.norm(loc=cf, scale=cf / 10.)
                    
                    # Fabricated convergence distribution
                    size = np.random.randint(3) + 2
                    
                    # Select SGCs from distribution, create, and connect to this cell
                    for j in pop.select(size=size, create=False, cf=dist):
                        sgc = pop.get_cell(j)
                        sgc.connect(cell)
                        
                else:
                    raise TypeError("Cannot connect population %s to %s" % (pop, self))

        # recursively resolve inputs in connected populations
        if depth > 1:
            for pop in self.connections:
                pop.resolve_inputs(depth-1)
        