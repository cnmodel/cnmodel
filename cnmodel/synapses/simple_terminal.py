from neuron import h

from .terminal import Terminal


class SimpleTerminal(Terminal):
    """
    Simple terminal using netcon.
    """
    def __init__(self, pre_sec, target_cell, spike_source=None, loc=0.5):
        """
        Parameters
        ----------
        pre_sec : :obj: `NEURON Section`
            The presynaptic section that is monitored for spikes. The voltage
            in this section is monitored to trigger the postsynaptic conductance
            as the spike source
        spike_source : :obj: `NEURON Section`
            Overrides the pre_sec as the spike source for this terminal.
        terminal : :obj: `Synapse Terminal`
            The presynaptic Terminal instance
        loc : float, default=0.5
            Position on the postsynaptic section to insert the mechanism, from [0..1].    
        """
        
        Terminal.__init__(self, pre_sec)
        if spike_source is None:
            spike_source = pre_sec(loc)._ref_v
        self.spike_source = spike_source
        self.pre_sec = pre_sec

    def connect(self, post, weight):
        """
        Connect this terminal to a postsynaptic cell section
        The synaptic delay is 0.5 msec, and the presynaptic
        action potential threshold is -20 mV.
        
        Parameters
        ----------
        post : :obj: `NEURON Section`
        
        weight : float
            Strength of the connection
        
        """
        thresh = -20
        delay = 0.5
        self.netcon = h.NetCon(self.spike_source, post, thresh, delay, weight, sec=self.pre_sec)
        self.netcon.weight[0] = weight
        