from .stochastic_terminal import StochasticTerminal
from .psd import PSD



class Synapse(object):
    def __init__(self, pre_cell, pre_opts, post_cell, post_opts, psd_type='multisite'):
        pre_opts['term_type'] = psd_type
#        post_opts['psd_type'] = psd_type # psd_type is not a valid keyword in make_psd...
        post_opts = {}
        
        terminal = pre_cell.make_terminal(post_cell, **pre_opts)
        psd = post_cell.make_psd(terminal, **post_opts)
        self.terminal = terminal
        self.psd = psd
        
