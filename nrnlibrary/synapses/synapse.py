from .stochastic_terminal import StochasticTerminal
from .psd import PSD



class Synapse(object):
    def __init__(self, terminal, psd):
        self.terminal = terminal
        self.psd = psd
        
    #def connect(self, pre_sec, post_sec, debug=False):

        #self.pre_sec = pre_sec
        #self.post_sec = post_sec
        #self.pre_cell = cells.cell_from_section(pre_sec)
        #self.post_cell = cells.cell_from_section(post_sec)

        
        ## ****************
        ## delay CV is about 0.3 (Sakaba and Takahashi, 2001), assuming delay is about 0.75
        ## This only adjusts the delay to each "terminal" from a single an - so is path length,
        ## but is not the individual site release latency. That is handled in COH4
        ##        if cellname == 'bushy':
        ##            delcv = 0.3*(delay/0.75) # scale by specified delay
        ##            newdelay = delay+delcv*np.random.standard_normal()
        ##            print "delay: %f   newdelay: %f   delcv: %f" % (delay, newdelay, delcv)
        ##            netcons[-1].delay = newdelay # assign a delay to EACH zone that is different...
        ##*****************
        #psd = PSD(pre_sec=pre_sec,
                  #post_sec=post_sec,
                  #terminal=term,
                  #eRev=0, debug=debug,
                  #gvar=0.3,
                  #nmda_ratio=0.0, identifier=1,
                  #) # set gVar to 0 for testing
        
        #self.terminal = term
        #self.psd = psd
        



