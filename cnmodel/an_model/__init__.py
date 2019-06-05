from .wrapper import model_ihc, model_synapse, seed_rng
try:
    from .wrapper import get_matlab
    MATLAB_FOUND = True
except:
    MATLAB_FOUND = False

from .cache import get_spiketrain
