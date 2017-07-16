import logging
import numpy as np

from .population import Population
from .. import cells


class SGC(Population):
    """A population of spiral ganglion cells.
    
    The cell distribution is uniform from 4kHz to 90kHz, evenly divided between
    spontaneous rate groups.
    """
    type = 'sgc'
    
    def __init__(self, species='mouse', model='dummy', **kwds):
        # Completely fabricated cell distribution: uniform from 4kHz to 90kHz.
        # Evenly divided between SR groups
        size = 10000
        next_seed = 0
        fields = [
            ('cf', float),
            ('sr', int),  # 0=low sr, 1=mid sr, 2=high sr
        ]
        super(SGC, self).__init__(species, size, fields=fields, model=model, **kwds)
        self._cells['cf'] = 4000 * 2**np.linspace(0, 4.5, size)
        self._cells['sr'] = np.arange(size) % 3
    
    def set_seed(self, seed):
        self.next_seed = seed
    
    def create_cell(self, cell_rec):
        """ Return a single new cell to be used in this population. The 
        *cell_rec* argument is the row from self.cells that describes the cell 
        to be created.
        """
        return cells.SGC.create(species=self.species, cf=cell_rec['cf'],
                                sr=cell_rec['sr'], **self._cell_args)
        
    def connect_pop_to_cell(self, pop, index):
        # SGC does not support any inputs
        assert len(self.connections) == 0

    def set_sound_stim(self, stim):
        """Set a sound stimulus to generate spike trains for all (real) cells
        in this population.
        """
        real = self.real_cells()
        logging.info("Assigning spike trains to %d SGC cells..", len(real))
        for i, ind in enumerate(real):
            logging.info("Assigning spike train to SGC %d (%d/%d)", ind, i, len(real))
            cell = self.get_cell(ind)
            cell.set_sound_stim(stim, self.next_seed)
            self.next_seed += 1
