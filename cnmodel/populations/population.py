import logging
import scipy.stats
import numpy as np

"""
Todo: 

* Need to work out most of the API details here, probably by starting with a
  specific use case and working backward

* Distributions of cell properties

* Mechanisms for automatically recording from neurons
    - record all Vm for all real neurons
    - record spike trains
    - record per-synapse currents
"""

class Population(object):
    """
    A Population represents a group of cell all having the same type. 

    Populations provide methods for:
    
    * Adding cells to the population with characteristic distributions.
    * Connecting the cells in one population to the cells in another.
    * Automatically adding cells to satisfy connectivity requirements when 
      connecting populations together.
    
    Populations have a concept of a "natural" underlying distribution of
    neurons, and behave as if all neurons in this distribution already exist
    in the model. However, initially all neurons are virtual, and are only 
    instantiated to become a part of the running model if the neuron provides 
    synaptic input to another non-virtual neuron, or if the user explicitly 
    requests a recording of the neuron.
    
    Subclasses represent populations for a specific cell type, and at least
    need to reimplement the `create_cell` and `connection_stats` methods.
    """
    def __init__(self, species, size, fields, **kwds):
        self._species = species
        self._post_connections = []  # populations this one connects to
        self._pre_connections = []  # populations connecting to this one
        
        # numpy record array with information about each cell in the 
        # population
        fields = [
            ('cell', object), 
            ('input_resolved', bool),
            ('connections', object),  # {pop: [cells], ...}
        ] + fields
        self._cells = np.zeros(size, dtype=fields)
        self._cell_indexes = {}
        self._cell_args = kwds

    @property
    def cells(self):
        """ The array of cells in this population. 
        
        For all populations, this array has a 'cell' field that is either 0
        (for virtual cells) or a Cell instance (for real cells). 
        
        Extra fields may be added by each Population subclass.
        """
        return self._cells.copy()
    
    @property
    def species(self):
        return self._species
    
    def unresolved_cells(self):
        """ Return indexes of all real cells whose inputs have not been 
        resolved.
        """
        real = self._cells['cell'] != 0
        unresolved = self._cells['input_resolved'] == False
        return np.argwhere(real & unresolved)[:,0]

    def real_cells(self):
        """ Return indexes of all real cells in this population.
        
        Initially, all cells in the population are virtual--they are accounted
        for, but not actually instantiated as part of the NEURON simulation.
        Virtual cells can be made real by calling `get_cell()`. This method
        returns the indexes of all cells for which `get_cell()` has already
        been invoked.
        """
        return np.argwhere(self._cells['cell'] != 0)[:,0]

    def connect(self, *pops):
        """ Connect this population to any number of other populations. 
        
        A connection is unidirectional; calling ``pop1.connect(pop2)`` can only
        result in projections from pop1 to pop2.
        
        Note that the connection is purely symbolic at first; no cells are 
        actually connected by synapses at this time.
        """
        self._post_connections.extend(pops)
        for pop in pops:
            pop._pre_connections.append(self)

    @property
    def pre_connections(self):
        """ The list of populations connected to this one.
        """
        return self._pre_connections[:]

    def cell_connections(self, index):
        """ Return a dictionary containing, for each population, a list of 
        cells connected to the cell in this population at *index*.
        """
        return self._cells[index]['connections']

    def resolve_inputs(self, depth=1):
        """ For each _real_ cell in the population, select a set of 
        presynaptic partners from each connected population and generate a 
        synapse from each.
        
        Although it is allowed to call ``resolve_inputs`` multiple times for
        a single population, each individual cell will only resolve its inputs
        once. Therefore, it is recommended to create and connect all 
        populations before making any calls to ``resolve_inputs``.
        """
        for i in self.unresolved_cells():
            # loop over all cells whose presynaptic inputs have not been resolved
            cell = self._cells[i]['cell']
            logging.info("Resolving inputs for %s %d", self, i)
            self._cells[i]['connections'] = {}
            
            # select cells from each population to connect to this cell
            for pop in self._pre_connections:
                pre_cells = self.connect_pop_to_cell(pop, i)
                logging.info("  connected %d cells from %s", len(pre_cells), pop)
                assert pre_cells is not None
                self._cells[i]['connections'][pop] = pre_cells
            self._cells[i]['input_resolved'] = True

        # recursively resolve inputs in connected populations
        if depth > 1:
            for pop in self.pre_connections:
                pop.resolve_inputs(depth-1)

    def connect_pop_to_cell(self, pop, cell_index):
        """ Connect cells in a presynaptic population to the cell in this 
        population at *cell_index*, and return the presynaptic indexes of cells
        that were connected.
        
        This method is responsible for choosing pairs of cells to be connected
        by synapses, and may be overridden in subclasses.
        
        The default implementation calls `self.connection_stats()` to determine
        the number and selection criteria of presynaptic cells.
        """
        cell_rec = self._cells[cell_index]
        cell = cell_rec['cell']
        
        size, dist = self.connection_stats(pop, cell_rec) 

        # Select SGCs from distribution, create, and connect to this cell
        # todo: select sgcs with similar spont. rate?
        pre_cells = pop.select(size=size, create=False, **dist)
        for j in pre_cells:
            pre_cell = pop.get_cell(j)
            # use default settings for connecting these. 
            pre_cell.connect(cell)
        return pre_cells

    def connection_stats(self, pop, cell_rec):
        """ The population *pop* is being connected to the cell described in 
        *cell_rec*.
        
        This method is responsible for deciding the distributions of presynaptic
        cell properties for any given postsynaptic cell (for example, a cell 
        with cf=10kHz might receive SGC input from 10 cells selected from a 
        normal distribution centered at 10kHz). 
        
        Most subclasses of Population should reimplement this method with the
        appropriate connectivity statistics for that cell type.
        
        Return a tuple (size, dist) with the following values:
        
        * size: integer giving the number of cells that should be selected from
          the presynaptic population and connected to the postsynaptic cell.
        * dist: dictionary of {property_name: distribution} pairs that describe
          how cells should be selected from the presynaptic population. See
          keyword arguments to `select()` for more information on the content
          of this dictionary.
        """
        raise NotImplementedError()
    
    def select(self, size, create=False, **kwds):
        """ Return a list of indexes for cells matching the selection criteria.
        
        The *size* argument specifies the number of cells to return.
        
        If *create* is True, then any selected cells that are virtual will be
        instantiated.
        
        Each keyword argument must be the name of a field in self.cells. Values
        may be:
        
        * A distribution (see scipy.stats), in which case the distribution 
          influences the selection of cells
        * An array giving the probability to assign to each cell in the
          population
        * A number, in which case the cell(s) with the closest match 
          are returned. If this is used, it overrides all other criteria except
          where they evaluate to 0.
        
        If multiple distributions are provided, then the product of the survival
        functions of all distributions determines the probability of selecting 
        each cell.
        """
        if len(kwds) == 0:
            raise TypeError("Must specify at least one selection criteria")
        
        full_dist = np.ones(len(self._cells))
        nearest = None
        nearest_field = None
        for field, dist in kwds.items():
            if np.isscalar(dist):
                if nearest is not None:
                    raise Exception("May not specify multiple single-valued selection criteria.")
                nearest = dist
                nearest_field = field
            elif isinstance(dist, scipy.stats.distributions.rv_frozen):
                vals = self._cells[field]
                full_dist *= dist.pdf(vals)
            elif isinstance(dist, np.ndarray):
                full_dist *= dist
            else:
                raise TypeError("Distributed criteria must be array or rv_frozen.")
                
        # Select cells nearest to the requested value, but only pick from 
        # cells with nonzero probability. 
        if nearest is not None:
            cells = []
            mask = full_dist == 0
            err = np.abs(self._cells[nearest_field] - nearest)
            for i in range(size):
                err[mask] = np.inf
                cell = np.argmin(err)
                mask[cell] = True
                cells.append(cell)
            
        # Select cells randomly from the specified combined probability 
        # distribution
        else:
            cells = []
            full_dist /= full_dist.sum()
            vals = np.random.uniform(size=size)
            vals.sort()
            cumulative = np.cumsum(full_dist)
            for val in vals:
                cell = np.argwhere(cumulative >= val)[0,0]
                cells.append(cell)
            
        if create:
            self.create_cells(cells)
        
        return cells

    def get_cell(self, i, create=True):
        """ Return the cell at index i. If the cell is virtual, then it will 
        be instantiated first unless *create* is False.
        """
        if create and self._cells[i]['cell'] == 0:
            self.create_cells([i])
        return self._cells[i]['cell']
    
    def get_cell_index(self, cell):
        """Return the index of *cell*.
        """
        return self._cell_indexes[cell]
        
    def create_cells(self, cell_inds):
        """ Instantiate each cell in *cell_inds*, which is a list of indexes into
        self.cells.
        """
        for i in cell_inds:
            if self._cells[i]['cell'] != 0:
                continue
            cell = self.create_cell(self._cells[i])
            self._cells[i]['cell'] = cell
            self._cell_indexes[cell] = i
            
    def create_cell(self, cell_rec):
        """ Return a single new cell to be used in this population. The 
        *cell_rec* argument is the row from self.cells that describes the cell 
        to be created.
        
        Subclasses must reimplement this method.
        """
        raise NotImplementedError()

    def __str__(self):
        return "<Population %s (%d/%d real)>" % (type(self).__name__, (self._cells['cell'] != 0).sum(), len(self._cells))

    def __getstate__(self):
        """Return a picklable copy of self.__dict__. 
        
        Note that we remove references to the actual cells in order to allow pickling.
        """
        state = self.__dict__.copy()
        state['_cells'] = state['_cells'].copy()
        mask = state['_cells']['cell'] != 0
        state['_cells'][mask] = [str(cell) for cell in state['_cells'][mask]]
        
        return state
        
        
