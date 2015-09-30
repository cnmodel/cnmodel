from neuron import h

from .terminal import Terminal


class SimpleTerminal(Terminal):
    """
    """
    def __init__(self, pre_sec, target_cell, spike_source=None):
        Terminal.__init__(self, pre_sec)
        if spike_source is None:
            spike_source = pre_sec(0.5)._ref_v
        self.spike_source = spike_source
        self.pre_sec = pre_sec

    def connect(self, post, weight):
        thresh = -20
        delay = 0.5
        self.netcon = h.NetCon(self.spike_source, post, thresh, delay, weight, sec=self.pre_sec)
        self.netcon.weight[0] = weight
        