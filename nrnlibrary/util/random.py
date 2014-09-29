import numpy as np
import hashlib, struct

_current_seed = 0

def set_seed(seed):
    """
    Set the random seed to be used globally. If a string is supplied, it 
    will be converted to int using hash().
    
    This immediately seeds the numpy RNG. Any other RNGs must be seeded using
    current_seed()
    """
    if isinstance(seed, str):
        seed = struct.unpack('=I', hashlib.md5(seed).digest()[:4])[0]
    np.random.seed(seed)
    assert seed < 2**64  # neuron RNG fails if seed is too large
    global _current_seed
    _current_seed = seed
    return seed
    
def current_seed():
    """
    Return the currently-set global random seed. 
    """
    return _current_seed
