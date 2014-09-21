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
    
    """
    def __init__(self, species, size, fields):
        self._species = species
        self._connections = []
        
        # numpy record array with information about each cell in the 
        # population
        fields = [('cell', object), ('resolved', bool)] + fields
        self._cells = np.empty(size, dtype=fields)
    
    def connect(self, *pops):
        """ Connect this population to any number of other populations. 
        
        A connection is unidirectional; calling ``pop1.connect(pop2)`` can only
        result in projections from pop1 to pop2.
        
        Note that the connection is purely symbolic at first; no cells are 
        actually connected by synapses at this time.
        """
        self._connections.extend(pops)

    def resolve_inputs(self):
        """ For each _real_ cell in the population, select a set of 
        presynaptic partners from each connected population and generate a 
        synapse from each.
        
        Although it is allowed to call ``resolve_inputs`` multiple times for
        a single population, each individual cell will only resolve its inputs
        once. Therefore, it is recommended to create and connect all 
        populations before making any calls to ``resolve_inputs``.
        
        """
    
    def select(self, **kwds):
        """ Return a list of neurons matching the selection criteria. Neurons
        may be either real or virtual.
        
        Individual populations must reimplement this method because the 
        availability of selection criteria may vary between cell types.
        """
        raise NotImplementedError()
        
