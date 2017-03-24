from neuron import h

from .terminal import Terminal


class SimpleTerminal(Terminal):
    """
    """
    def __init__(self, pre_sec, target_cell, spike_source=None, loc=0.5):
        Terminal.__init__(self, pre_sec)
        if spike_source is None:
            spike_source = pre_sec(loc)._ref_v
        self.spike_source = spike_source

    def connect(self, post, weight):
        thresh = -20
        delay = 0.5
        self.netcon = h.NetCon(self.spike_source, post, thresh, delay, weight)
        self.netcon.weight[0] = weight
        