from .stochastic_terminal import StochasticTerminal
from .psd import PSD

from .. import cells

# utility class to create parameter lists... 
# create like: p = Params(abc=2.0, defg = 3.0, lunch='sandwich')
# reference like p.abc, p.defg, etc.
class Params(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

class Synapse(object):
    def connect(self, pre_sec, post_sec, debug=False):

        self.pre_sec = pre_sec
        self.post_sec = post_sec
        self.pre_cell = cells.cell_from_section(pre_sec)
        self.post_cell = cells.cell_from_section(post_sec)

        #
        # set parameters according to the target cell type
        #
        if isinstance(self.post_cell, cells.TStellate):
            nANTerminals_ReleaseZones = 1
        elif isinstance(self.post_cell, cells.DStellateEager):
            nANTerminals_ReleaseZones = 1
        elif isinstance(self.post_cell, cells.DStellate):
            nANTerminals_ReleaseZones = 1
        elif isinstance(self.post_cell, cells.Bushy):
            nANTerminals_ReleaseZones = 100
        else:
            raise ValueError("Unknown cell type '%s'" % self.post_cell)

        self.zones_per_terminal = nANTerminals_ReleaseZones
        #
        # make a multisite synapse (consists of an axon section and a coh point process)
        # synapse is like a calyx or end bulb of Held
        #
        #calyx_zones = 1
        #debug = True
        #synType = 'epsc'
        #ipscdynamics = 1

        ANTerminals_Delay = 0.0 # latency between AP and mode of release distribution, in milliseconds.    
        ANTerminals_Latency = 0.5 # latency 
        thresh = -30 # mV - AP detection on the presynaptic side.

        # set parameter control for the stochastic release of vesicles...
        # this structure is passed to stochastic synapses, and replaces several variables 
        # that were previously defined in the call to that function.
        vPars = Params(LN_Flag=1, LN_t0=10.0, LN_A0=0.05, LN_tau=35, LN_std=0.05,
                    Lat_Flag=1, Lat_t0=10.0, Lat_A0=0.140, Lat_tau=21.5,
                    latency=ANTerminals_Latency)

        #printParams(vPars)
        # LC: disabled this -- let the user decide which sections to connect.
        #if hasattr(post_cell, 'len') and len(post_cell) > 1: # handle possible multi-segment cells
            #TC = post_cell[0].soma # soma
            #SOMA = post_cell[0].soma
        #else:
            #TC = post_cell.soma
            #SOMA = post_cell.soma
        #if isinstance(post_cell, cells.DStellateEager): 
            ## well, these always have a dendrite... 
            #TC = post_cell.dendrite
            #SOMA = post_cell.soma


        # create a presynaptic ending with nzones release sites
        #NOTE: stochastic_pars must define parameters used by multisite, including:
            #.delay is the netcon delay between the presynaptic AP and the start of release events
            #.Latency is the latency to the mean release event... this could be confusing.
        term = StochasticTerminal(parent_section=pre_sec,
                                    target_cell=self.post_cell,
                                    nzones=nANTerminals_ReleaseZones,
                                    stochastic_pars=vPars,
                                    calcium_pars=None,
                                    identifier=1, debug=debug)
        
        # ****************
        # delay CV is about 0.3 (Sakaba and Takahashi, 2001), assuming delay is about 0.75
        # This only adjusts the delay to each "terminal" from a single an - so is path length,
        # but is not the individual site release latency. That is handled in COH4
        #        if cellname == 'bushy':
        #            delcv = 0.3*(delay/0.75) # scale by specified delay
        #            newdelay = delay+delcv*np.random.standard_normal()
        #            print "delay: %f   newdelay: %f   delcv: %f" % (delay, newdelay, delcv)
        #            netcons[-1].delay = newdelay # assign a delay to EACH zone that is different...
        #*****************
        psd = PSD(pre_sec=pre_sec,
                    post_sec=post_sec,
                    terminal=term,
                    eRev=0, debug=debug,
                    delay=ANTerminals_Delay,
                    thresh=thresh, gvar=0.3,
                    nmda_ratio=0.0, identifier=1,
                    ) # set gVar to 0 for testing
        
        self.terminal = term
        self.psd = psd
        



