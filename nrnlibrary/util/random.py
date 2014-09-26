import numpy as np

_current_seed = 0

def set_seed(seed):
    """
    Set the random seed to be used globally. If a string is supplied, it 
    will be converted to int using hash().
    
    This immediately seeds the numpy RNG. Any other RNGs must be seeded using
    current_seed()
    """
    if isinstance(seed, str):
        seed = abs(hash(seed))
    np.random.seed(seed)
    global _current_seed
    _current_seed = seed
    return seed
    
def current_seed():
    """
    Return the currently-set global random seed. 
    """
    print "SEED:", _current_seed
    return _current_seed
