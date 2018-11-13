import logging
import numpy as np
import pyqtgraph.multiprocess as mp

from .population import Population
from .. import cells


class SGC(Population):
    """A population of spiral ganglion cells.
    
    The cell distribution is uniform from 2kHz to 64kHz, evenly divided between
    spontaneous rate groups.
    """
    type = 'sgc'
    
    def __init__(self, species='mouse', model='dummy', **kwds):
        # Completely fabricated cell distribution: uniform from 2kHz to 40kHz,
        # evenly divided between SR groups. We only go up to 40kHz because the
        # auditory periphery model does not support >40kHz.
        freqs = self._get_cf_array(species)
        fields = [
            ('cf', float),
            ('sr', int),  # 0=low sr, 1=mid sr, 2=high sr
        ]
        super(SGC, self).__init__(species, len(freqs), fields=fields, model=model, **kwds)
        self._cells['cf'] = freqs
        # evenly distribute SR groups
        self._cells['sr'] = np.arange(len(freqs)) % 3
    
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

    def set_sound_stim(self, stim, parallel=False):
        """Set a sound stimulus to generate spike trains for all (real) cells
        in this population.
        """
        real = self.real_cells()
        logging.info("Assigning spike trains to %d SGC cells..", len(real))
        if not parallel:
            for i, ind in enumerate(real):
                #logging.info("Assigning spike train to SGC %d (%d/%d)", ind, i, len(real))
                cell = self.get_cell(ind)
                cell.set_sound_stim(stim, self.next_seed)
                self.next_seed += 1

        else:
            seeds = range(self.next_seed, self.next_seed + len(real))
            self.next_seed = seeds[-1] + 1
            tasks = zip(seeds, real)
            trains = [None] * len(tasks)
            # generate spike trains in parallel
            with mp.Parallelize(enumerate(tasks), trains=trains, progressDialog='Generating SGC spike trains..') as tasker:
                for i, x in tasker:
                    seed, ind = x
                    cell = self.get_cell(ind)
                    train = cell.generate_spiketrain(stim, seed)
                    tasker.trains[i] = train
            # collected all trains; now assign to cells
            for i,ind in enumerate(real):
                cell = self.get_cell(ind)
                cell.set_spiketrain(trains[i])
            