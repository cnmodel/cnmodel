class Terminal(object):
    """
    Base class for axon terminals. A terminal has a single postsynaptic 
    neuron, but may have multiple release zones. It defines a release mechanism
    with a NetCon input (triggering from presynaptic voltage or calcium level)
    and either NetCon or pointer output for driving a PSD.
    
    """
