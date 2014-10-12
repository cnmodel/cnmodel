from .stochastic_terminal import StochasticTerminal
from .psd import PSD



class Synapse(object):
    def __init__(self, terminal, psd):
        self.terminal = terminal
        self.psd = psd
        
