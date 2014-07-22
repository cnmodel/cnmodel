#!/usr/bin/python
#
# Synapse definitions for models.
#
# This file includes a number of different synapse definitions and default conductances
# for point models. Most are models from the lab for neurons of the cochlear nucleus.
# the synaptic receptor models are gleaned from the literature and sometimes fitted to the 
# cochlear nucleus data.
#
# Paul B. Manis, Ph.D. 2009 (August - November 2009)
#
from neuron import h
import numpy as np 

from .. import util
from .. import cells


runQuiet = True
veryFirst = 1
# utility class to create parameter lists... 
# create like: p = Params(abc=2.0, defg = 3.0, lunch='sandwich')
# reference like p.abc, p.defg, etc.
class Params(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


#def template_SGAxon(debug=False, message=None):
    #axon = cells.HH(message='')
    #return (axon,)

# define the "Calyx" template.
# Each axon ends in (and thus includes) a calyx; therefore we make a list of axons
# based on the template.
# incoming argument: which axon (index) to address SGs
# Technically, this isn't the axon, but the terminal....

def template_Calyx_Billup(debug=False, nzones=1, message=None):
    calyx = cells.HH(message='  >> creating calyx_Billup')
    coh = []
    calyx.push()
    for k in range(0, nzones):
        coh.append(h.COH2(0.5, calyx))
    h.pop_section()
    return (calyx, coh)




def template_iGluR_PSD(sec, nReceptors=1, debug=False, cellname=None, message=None, nmda_ratio=1):
    """
    Create an ionotropic Glutamate receptor "PSD"
    Each PSD has receptors for each active zone, which must be matched (connected) to presynaptic
    terminals. Each PSD recetpor consists of an AMPATRUSSELL and an NMDA_KAMPA receptor
    Inputs:
        sec: The template requires a segment to insert the receptors into
        nReceptors: The number of receptor sites to insert
        debug: flag for debugging (prints extra information)
        cellname: Bushy/MNTB/stellate: determines ampa receptor kinetics
        message: Not used.
        nmda_ratio: The relative conductance of the open NMDA receptors to the open AMPA receptors.
    Outputs:
        (psd, psdn, par, parn)
        psd is the list of PSDs that were created (AMPA)
        psdn is the list of NMDA PSDs (same number as psd, just the NMDARs)
        par: dictionary of AMPAR kinetics as inserted
        parn: NMDA ratio
    Side Effecdts: None
    """
    psd = []
    psdn = []
    sec.push()
    AN_Po_Ratio = 23.2917 # ratio of open probabilities for AMPA and NMDAR's at peak currents
    for k in range(0, nReceptors):
        psd.append(h.AMPATRUSSELL(0.5, sec)) # raman/trussell AMPA with rectification
        psdn.append(h.NMDA_Kampa(0.5, sec)) # Kampa state model NMDA receptors

        if cellname in ['bushy', 'MNTB']:
            psd[-1].Ro1 = 107.85
            psd[-1].Ro2 = 0.6193
            psd[-1].Rc1 = 3.678
            psd[-1].Rc2 = 0.3212
            psdn[-1].gNAR = 0.036 * AN_Po_Ratio * nmda_ratio # 0.36*AN_Po_Ratio*nmda_ratio
            #if k == 0:
            #    print "Bushy NMDAR set to %8.2f" % psdn[-1].gNAR
        if cellname == 'stellate':
            psd[-1].Ro1 = 39.25
            psd[-1].Ro2 = 4.40
            psd[-1].Rc1 = 0.667
            psd[-1].Rc2 = 0.237
            psd[-1].PA = 0.1
            psdn[-1].gNAR = 1 * AN_Po_Ratio * nmda_ratio

    h.pop_section()
    par = {'Ro1': ('r', psd[0].Ro1),
           'Ro2': ('r', psd[0].Ro2),
           'Rc1': ('r', psd[0].Rc1),
           'Rc2': ('r', psd[0].Rc2), }
    parn = {'gNAR': ('4', psdn[0].gNAR), }
    return (psd, psdn, par, parn)


# the following templates are a bit more complicated.
# The parameter names as defined in the model are returned
# but also those that are involved in the forward binding reactions are
# listed separately - this is to allow us to run curve fits that adjust
# only a subset of the parameters a  time - e.g., the rising phase, then
# the falling phase with the rising phase fixed.
# the dictionary selection is made by selectpars in glycine_fit.py.
#

def template_Gly_PSD_exp(sec, debug=False, nReceptors=2, cellname=None, message=None):
    psd = []
    sec.push()
    for k in range(0, nReceptors):
        psd.append(h.GLY2(0.5, sec))
    h.pop_section()
    par = ['alpha', 'beta']
    p = []
    for n in par:
        p.append(eval('psd[0].' + n)) # get default values from the mod file
    return (psd, par, p)


def template_Gly_PSD_State_Glya5(sec, debug=False, nReceptors=2, psdtype=None, message=None):
    psd = []
    sec.push()
    for k in range(0, nReceptors):
        psd.append(h.GLYa5(0.5, sec))
    h.pop_section()
    par = {'kf1': ('r', psd[0].kf1), # retreive values in the MOD file
           'kf2': ('r', psd[0].kf2),
           'kb1': ('r', psd[0].kb1),
           'kb2': ('r', psd[0].kb2),
           'a1': ('f', psd[0].a1),
           'b1': ('f', psd[0].b1),
           'a2': ('f', psd[0].a2),
           'b2': ('f', psd[0].b2)}
    return (psd, par)


def template_Gly_PSD_State_Gly6S(sec, debug=False, nReceptors=2, psdtype=None, message=None):
    psd = []
    sec.push()
    for k in range(0, nReceptors):
        psd.append(h.Gly6S(0.5, sec)) # simple using Trussell model 6 states with desens
        if not runQuiet:
            print 'Gly6S psdtype: ', psdtype
        if psdtype == 'glyslow': # fit on 8 March 2010, error =  0.164, max open: 0.155
            psd[-1].Rd = 1.177999
            psd[-1].Rr = 0.000005
            psd[-1].Rb = 0.009403
            psd[-1].Ru2 = 0.000086
            psd[-1].Ro1 = 0.187858
            psd[-1].Ro2 = 1.064426
            psd[-1].Ru1 = 0.028696
            psd[-1].Rc1 = 0.103625
            psd[-1].Rc2 = 1.730578
    h.pop_section()
    par = {'Rb': ('n', psd[0].Rb), # retrive values in the MOD file
           'Ru1': ('r', psd[0].Ru1),
           'Ru2': ('r', psd[0].Ru2),
           'Rd': ('f', psd[0].Rd),
           'Rr': ('f', psd[0].Rr),
           'Ro1': ('f', psd[0].Ro1),
           'Ro2': ('r', psd[0].Ro2),
           'Rc1': ('f', psd[0].Rc1),
           'Rc2': ('r', psd[0].Rc2)}
    return (psd, par)


def template_Gly_PSD_State_PL(sec, debug=False, nReceptors=2, cellname=None,
                              psdtype=None, message=None):
    psd = []
    sec.push()
    for k in range(0, nReceptors):
        psd.append(h.GLYaPL(0.5, sec)) # simple dextesche glycine receptors
        if not runQuiet:
            print 'PL psdtype: ', psdtype
        if psdtype == 'glyslow':
            psd[-1].a1 = 0.000451
            psd[-1].a2 = 0.220
            psd[-1].b1 = 13.27
            psd[-1].b2 = 6.845
            psd[-1].kon = 0.00555
            psd[-1].koff = 2.256
            psd[-1].r = 1.060
            psd[-1].d = 55.03
        if psdtype == 'glyfast': # fit from 3/5/2010. error =  0.174 maxopen = 0.0385
            psd[-1].a1 = 1.000476
            psd[-1].a2 = 0.137903
            psd[-1].b1 = 1.700306
            psd[-1].koff = 13.143132
            psd[-1].kon = 0.038634
            psd[-1].r = 0.842504
            psd[-1].b2 = 8.051435
            psd[-1].d = 12.821820
    h.pop_section()
    par = {'kon': ('r', psd[0].kon), # retrive values in the MOD file
           'koff': ('r', psd[0].koff),
           'a1': ('r', psd[0].a1),
           'b1': ('r', psd[0].b1),
           'a2': ('f', psd[0].a2),
           'b2': ('f', psd[0].b2),
           'r': ('f', psd[0].r),
           'd': ('f', psd[0].d)}
    return (psd, par)


def template_Gly_PSD_State_GC(sec, debug=False, nReceptors=2, psdtype=None, message=None):
    psd = []
    sec.push()
    for k in range(0, nReceptors):
        psd.append(h.GLYaGC(0.5, sec)) # simple dextesche glycine receptors
        if psdtype == 'glyslow':
            psd[-1].k1 = 12.81    # (/uM /ms)	: binding
            psd[-1].km1 = 0.0087#	(/ms)	: unbinding
            psd[-1].a1 = 0.0195#	(/ms)	: opening
            psd[-1].b1 = 1.138    #(/ms)	: closing
            psd[-1].r1 = 6.13    #(/ms)	: desense 1
            psd[-1].d1 = 0.000462 # (/ms)	: return from d1
            psd[-1].r2 = 0.731# (/ms)     : return from deep state
            psd[-1].d2 = 1.65 #(/ms)     : going to deep state
            psd[-1].r3 = 3.83 #(/ms)     : return from deep state
            psd[-1].d3 = 1.806 #(/ms)     : going to deep state
            psd[-1].rd = 1.04 #(/ms)
            psd[-1].dd = 1.004 # (/ms)
    h.pop_section()
    par = {'k1': ('r', psd[0].k1), # retrive values in the MOD file
           'km1': ('r', psd[0].km1),
           'a1': ('r', psd[0].a1),
           'b1': ('r', psd[0].b1),
           'r1': ('f', psd[0].r1),
           'd1': ('f', psd[0].d1),
           'r2': ('f', psd[0].r2),
           'd2': ('f', psd[0].d2),
           'r3': ('f', psd[0].r3),
           'd3': ('f', psd[0].d3),
           'rd': ('f', psd[0].rd),
           'dd': ('f', psd[0].dd)}
    return (psd, par)





class Synapse(object):
    def connect(self, pre_sec, post_sec, debug=False):
        
        self.AN_Po_Ratio = 23.2917 # ratio of open probabilities for AMPA and NMDAR's at peak currents
        self.AMPA_Max_Po = 0.44727
        self.NMDARatio = 0.0

        self.pre_sec = pre_sec
        self.post_sec = post_sec
        self.pre_cell = cells.cell_from_section(pre_sec)
        self.post_cell = cells.cell_from_section(post_sec)

        #
        # set parameters according to the target cell type
        #
        if isinstance(self.post_cell, cells.TStellate):
            nANTerminals_ReleaseZones = 1
            AN_gMax = 4600.0 / self.AMPA_Max_Po # correct because Po in model is not 1.0

        elif isinstance(self.post_cell, cells.DStellateEager):
            nANTerminals_ReleaseZones = 1
            AN_gMax = 4600.0 / self.AMPA_Max_Po # correct because Po in model is not 1.0

        elif isinstance(self.post_cell, cells.DStellate):
            nANTerminals_ReleaseZones = 1
            AN_gMax = 4600.0 / self.AMPA_Max_Po # correct because Po in model is not 1.0

        elif isinstance(self.post_cell, cells.Bushy):
            nANTerminals_ReleaseZones = 100
            # normally an_gmax would be about 22 nS (22000) total
            # However, if we do this as a "per release site" conductance, we have ~100 pA / 60 mV = 1.7 nS
            #AN_gMax = 22000.0/(AMPA_Max_Po) # correct because Po in model is not 1.0
            AN_gMax = 1700.0 / (self.AMPA_Max_Po)

        else:
            raise ValueError("Unknown cell type '%s'" % self.post_cell)

        self.zones_per_terminal = nANTerminals_ReleaseZones
        #
        # make a multisite synapse (consists of an axon section and a coh point process)
        # synapse is like a calyx or end bulb of Held
        #
        #calyx_zones = 1
        #debug = True
        self.psdType = 'ampa'
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
                                    identifier=1, debug=debug,
                                    psdtype=self.psdType)
        
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
        psd = PSD(parent_section=pre_sec,
                    target_section=post_sec,
                    terminal=term,
                    eRev=0, debug=debug,
                    delay=ANTerminals_Delay,
                    thresh=thresh, psdtype=self.psdType, gmax=AN_gMax, gvar=0.3,
                    nmda_ratio=0.0, identifier=1,
                    ) # set gVar to 0 for testing
        
        self.terminal = term
        self.psd = psd
        
        # adjust NMDA receptors to match postsynaptic cell
        self.adjust_nmda()


    def adjust_nmda(self):
        k = 0
        kNMDA = -1
        kAMPA = -1
        for p in self.psd.psd:
            if p.hname().find('NMDA', 0, 6) >= 0:
                if isinstance(self.post_cell, cells.TStellate):
                    p.gNAR = 1.28 * self.AN_Po_Ratio * self.NMDARatio # for T-stellate cells, this yields correct AMPA, NMDA ratio of 1.13 at +40 mV
                elif isinstance(self.post_cell, cells.Bushy):
                    p.gNAR = 0.36 * self.AN_Po_Ratio * self.NMDARatio # for bushy cells, this yields correct AMPA, NMDA ratio of 0.429 at +40 mV
                    #if p is psd[0]:
                    #    print "NMDAR's for bushy cells set to %8.3f" % p.gNAR
                elif isinstance(self.post_cell, cells.DStellate):
                    p.gNAR = 1.28 * self.AN_Po_Ratio * self.NMDARatio # same as for T-stellate (no other data)
                p.vshift = 0
                if kNMDA == -1:
                    kNMDA = k # save the first instance where we have an NMDA receptor
            else:
                if kAMPA == -1: # not NMDA, so get AMPA 
                    kAMPA = k
            k = k + 1

        self.kNMDA = kNMDA
        self.kAMPA = kAMPA



class Terminal(object):
    """
    Base class for axon terminals. A terminal has a single postsynaptic 
    neuron, but may have multiple release zones. It defines a release mechanism
    with a NetCon input (triggering from presynaptic voltage or calcium level)
    and either NetCon or pointer output for driving a PSD.
    
    """
    
    
class StochasticTerminal(object):
    """
    Axon terminal with multi-site sctochastic release mechanism.    
    """
    def __init__(self, parent_section, target_cell, nzones=1, celltype='bushy', message=None,
                type='lognormal', identifier=0, stochastic_pars=None, calcium_pars=None,
                debug=False, psdtype=None, select=None):
        """
        This routine creates a (potentially) multisite synapse with:
            A MultiSiteSynapse release mechanism that includes stochastic release, with a lognormal
                release latency distribution.
            A "cleft" mechanism (models diffusion of transmitter). Note that the cleft is inserted as part of the
                presynaptic section, but is not connected to the postsynaptic side yet.
        Inputs:
            parent_section: the section where the synaptic mechanisms should be inserted.
            nzones: the number of activate zones to insert into the section.
            celltype: bushy or (anything else), sets the duration and amplitude of the transmitter transient
                generated by these synapses
            message: a message to type out when instantiating (mostly for verification of code flow)
            type: 'lognormal' sets the release event latency distribution to use a lognormal function. Currently,
                no other function is supported.
            identifier: an identifier to associate with these release sites so we can find them later.
            stochastic_pars: A dictionary of parameters (Param class) used to specifiy the stochastic behavior of this site,
                including release latency, stdev, and lognormal distribution paramaters
            calcium_pars: A dictionary of parameters (Param class) to determine the calcium channels in this section.
                If None, then no calcium channels are inserted; otherwise, a P-type calcium conductance and a dynamic
                mechanism are inserted, and their conductance is set.
        Outputs: a list with 3 variables:
        terminal, relsite, cleft
            terminal: this is the pointer to the terminal section that was inserted (same as parent_section if it was
                specified)
            relsite: a list of the nzones release sites that were created
            cleft: a list of the nzones cleft mechanisms that were created.
        """
        if stochastic_pars is None:
            raise TypeError
        global veryFirst
        mu = u'\u03bc'
        sigma = u'\u03c3'
        message='  >> creating terminal with %d release zones using lognormal release latencies (coh4)' % nzones
        if debug:
            print message
        terminal = parent_section
        terminal.push()
        if calcium_pars is not None:
            terminal.insert('cap') # insert calcium channel density
            terminal().cap.pcabar = calcium_pars.Ca_gbar
            terminal.insert('cad')
            
        # Create point process to simulate multiple independent release zones.
        relsite = h.MultiSiteSynapse(0.5, sec=terminal)
        relsite.nZones = nzones
        relsite.rseed = 2 # int(np.random.random_integers(1,1024))
        relsite.latency = stochastic_pars.latency
        relsite.latstd = stochastic_pars.LN_std
        if debug is True:
            relsite.debug = 1
        relsite.Identifier = identifier
        # if type == 'gamma':
        #     gd = gamma.rvs(2, size=10000)/2.0 # get a sample of 10000 events with a gamma dist of 2, set to mean of 1.0
        #     if relsite.latstd > 0.0:
        #         gds = relsite.latency+std*(gd-1.0)/gd.std() # scale standard deviation
        #     else:
        #         gds = relsite.latency*np.ones((10000,1))
        # if type == 'lognormal':
        #     if std > 0.0:
        #         gds = lognormal(mean=0, sigma=relsite.latstd, size=10000)
        #     else:
        #         gds = np.zeros((10000, 1))
        # use the variable latency mode of COH4. And, it is lognormal no matter what.
        # the parameters are defined in COH4.mod as follows
        #    Time course of latency shift in release during repetitive stimulation
        #   Lat_Flag = 0 (1) : 0 means fixed latency, 1 means lognormal distribution
        #   Lat_t0 = 0.0 (ms) : minimum time since simulation start before changes in latency are calculated
        #   Lat_A0 = 0.0 (ms) : size of latency shift from t0 to infinity
        #   Lat_tau = 100.0 (ms) : rate of change of latency shift (from fit of a+b(1-exp(-t/tau)))
        #   : Statistical control of log-normal release shape over time during repetive stimulation
        #   LN_Flag = 0 (1) : 0 means fixed values for all time
        #   LN_t0 = 0.0 (ms) : : minimum time since simulation start before changes in distribution are calculated
        #   LN_A0 = 0.0 (ms) : size of change in sigma from t0 to infinity
        #   LN_tau = 100.0 (ms) : rate of change of sigma over time (from fit of a+b*(1-exp(-t/tau)))

        relsite.LN_Flag = stochastic_pars.LN_Flag # enable use of lognormal release latency
        relsite.LN_t0 = stochastic_pars.LN_t0
        relsite.LN_A0 = stochastic_pars.LN_A0
        relsite.LN_tau = stochastic_pars.LN_tau
        relsite.Lat_Flag = stochastic_pars.Lat_Flag
        relsite.Lat_t0 = stochastic_pars.Lat_t0
        relsite.Lat_A0 = stochastic_pars.Lat_A0
        relsite.Lat_tau = stochastic_pars.Lat_tau
        if veryFirst == 1 and debug is True:
            veryFirst = 0
            #mpl.figure(2)
        if celltype in ['bushy', 'MNTB']:
            relsite.TDur = 0.10
            relsite.TAmp = 0.770
        else: # stellate
            relsite.TDur = 0.25
            relsite.TAmp = 1.56625
        h.pop_section()
        self.relsite = relsite
        self.n_rzones = nzones

        if psdtype is not None:
            self.setPsdType(psdtype, target_cell, select)

    def setPsdType(self, psdtype, target_cell, select=None):
        if psdtype == 'ampa':
            self.setDF(target_cell, 'epsc') # set the parameters for release
        elif psdtype.startswith('gly'):
            self.setDF(target_cell, 'ipsc', select) # set the parameters for release
        

    ################################################################################
    # The following routines set the synapse dynamics, based on measurements and fit
    # to the Dittman-Regehr model.
    ################################################################################

    def setDF(self, target_cell, synapsetype, select=None):
        """ set the parameters for the calyx release model ...
            These paramteres were obtained from an optimized fit of the Dittman-Regehr
            model to stimulus and recovery data for the synapses at 100, 200 and 300 Hz,
            for times out to about 0.5 - 1.0 second. Data from Ruili Xie and Yong Wang.
            Fitting by Paul Manis
        """
        if isinstance(target_cell, cells.Bushy):
            if synapsetype == 'epsc':
                self.bushy_epsc()
            if synapsetype == 'ipsc':
                if select is None:
                    self.bushy_ipsc_average()
                else:
                    self.bushy_ipsc_single(select=select)
        elif isinstance(target_cell, cells.TStellate):
            if synapsetype == 'epsc':
                self.stellate_epsc()
            if synapsetype == 'ipsc':
                self.stellate_ipsc()

    def bushy_epsc(self):
        """ data is average of 3 cells studied with recovery curves and individually fit """
        self.relsite.F = 0.29366
        self.relsite.k0 = 0.52313 / 1000.0
        self.relsite.kmax = 19.33805 / 1000.0
        self.relsite.taud = 15.16
        self.relsite.kd = 0.11283
        self.relsite.taus = 17912.2
        self.relsite.ks = 11.531
        self.relsite.kf = 17.78
        self.relsite.tauf = 9.75
        self.relsite.dD = 0.57771
        self.relsite.dF = 0.60364
        self.relsite.glu = 2.12827

    def stellate_epsc(self):
        """ data is average of 3 cells studied with recovery curves and individually fit """
        self.relsite.F = 0.43435
        self.relsite.k0 = 0.06717 / 1000.0
        self.relsite.kmax = 52.82713 / 1000.0
        self.relsite.taud = 3.98
        self.relsite.kd = 0.08209
        self.relsite.taus = 16917.120
        self.relsite.ks = 14.24460
        self.relsite.kf = 18.16292
        self.relsite.tauf = 11.38
        self.relsite.dD = 2.46535
        self.relsite.dF = 1.44543
        self.relsite.glu = 5.86564


    def stellate_ipsc(self):
        """ data is average of 3 cells studied with recovery curves and individually fit, 100 Hz """
        self.relsite.F = 0.23047
        self.relsite.k0 = 1.23636 #/ 1000.0
        self.relsite.kmax = 45.34474 #/ 1000.0
        self.relsite.taud = 98.09
        self.relsite.kd = 0.01183
        self.relsite.taus = 17614.50
        self.relsite.ks = 17.88618
        self.relsite.kf = 19.11424
        self.relsite.tauf = 32.28
        self.relsite.dD = 2.52072
        self.relsite.dF = 2.33317
        self.relsite.glu = 3.06948

    def bushy_ipsc_average(self):
        """average of 16 Bushy cells. Done differently than other averages.
        The individual fits were compiled, and an average computed for just the 100 Hz data
        across the individual fits. This average was then fit to the function of Dittman and Regeher
        (also in Xu-Friedman's papers). 
        The individual cells show a great deal of variability, from straight depression, to 
        depression/facilitaiton mixed, to facilation alone. This set of parameters generates
        a weak facilitaiton followed by depression back to baseline.
        """
        print "USING average kinetics for Bushy IPSCs"

        # average of 16cells for 100 Hz (to model); no recovery.
        self.relsite.F = 0.18521
        self.relsite.k0 = 2.29700
        self.relsite.kmax = 27.6667
        self.relsite.taud = 0.12366
        self.relsite.kd = 0.12272
        self.relsite.taus = 9.59624
        self.relsite.ks = 8.854469
        self.relsite.kf = 5.70771
        self.relsite.tauf = 0.37752
        self.relsite.dD = 4.00335
        self.relsite.dF = 0.72605
        self.relsite.glu = 5.61985

        # estimates for 400 Hz
        # self.relsite.F = 0.09
        # self.relsite.k0 = 1.2;
        # self.relsite.kmax = 30.;
        # self.relsite.taud = 0.01
        # self.relsite.kd = 0.75
        # self.relsite.taus = 0.015
        # self.relsite.ks = 1000.0
        # self.relsite.kf = 5.0
        # self.relsite.tauf = 0.3
        # self.relsite.dD = 1.0
        # self.relsite.dF = 0.025
        # self.relsite.glu = 4

        # self.relsite.F = 0.085426
        # self.relsite.k0 = 1.199372
        # self.relsite.kmax = 24.204277
        # self.relsite.taud = 0.300000
        # self.relsite.kd = 1.965292
        # self.relsite.taus = 2.596443
        # self.relsite.ks = 0.056385
        # self.relsite.kf = 0.721157
        # self.relsite.tauf = 0.034560
        # self.relsite.dD = 0.733980
        # self.relsite.dF = 0.025101
        # self.relsite.glu = 3.877192

        # average of 8 cells, all at 100 Hz with no recovery
        # self.relsite.F =    0.2450
        # self.relsite.k0 =   1.6206/1000.0
        # self.relsite.kmax = 26.0607/1000.0
        # self.relsite.taud = 0.0798
        # self.relsite.kd =   0.9679
        # self.relsite.taus = 9.3612
        # self.relsite.ks =   14.3474
        # self.relsite.kf =    4.2168
        # self.relsite.tauf =  0.1250
        # self.relsite.dD =    4.2715
        # self.relsite.dF =    0.6322
        # self.relsite.glu =   8.6160

        # average of 5 cells, mostly 100 Hz, but some 50, 200 and 400
        # self.relsite.F =    0.15573
        # self.relsite.k0 =   2.32272/1000.
        # self.relsite.kmax = 28.98878/1000.
        # self.relsite.taud = 0.16284
        # self.relsite.kd =   2.52092
        # self.relsite.taus = 17.97092
        # self.relsite.ks =   19.63906
        # self.relsite.kf =   7.44154
        # self.relsite.tauf = 0.10193
        # self.relsite.dD =   2.36659
        # self.relsite.dF =   0.38516
        # self.relsite.glu =  8.82600

        #original average - probably skewed.
        # self.relsite.F = 0.23382
        # self.relsite.k0 = 0.67554/1000.0
        # self.relsite.kmax = 52.93832/1000.0
        # self.relsite.taud = 8.195
        # self.relsite.kd = 0.28734
        # self.relsite.taus = 17.500
        # self.relsite.ks = 4.57098
        # self.relsite.kf = 16.21564
        # self.relsite.tauf = 123.36
        # self.relsite.dD = 2.21580
        # self.relsite.dF = 1.17146
        # self.relsite.glu = 1.90428


    def bushy_ipsc_single(self, select=None):
        """ data is from 31aug08b (single cell, clean dataset)"""
        print "Using bushy ipsc"

        if select is None or select > 4 or select <= 0:
            bushy_ipsc_average()
            return

        if select is 1: # 30aug08f
            print "using 30aug08f ipsc"
            self.relsite.F = 0.221818
            self.relsite.k0 = 0.003636364
            self.relsite.kmax = 0.077562107
            self.relsite.taud = 0.300000
            self.relsite.kd = 1.112554
            self.relsite.taus = 3.500000
            self.relsite.ks = 0.600000
            self.relsite.kf = 3.730452
            self.relsite.tauf = 0.592129
            self.relsite.dD = 0.755537
            self.relsite.dF = 2.931578
            self.relsite.glu = 1.000000

        if select is 2: #30aug08h
            print "using 30aug08H ipsc"
            self.relsite.F = 0.239404
            self.relsite.k0 = 3.636364 / 1000.
            self.relsite.kmax = 16.725479 / 1000.
            self.relsite.taud = 0.137832
            self.relsite.kd = 0.000900
            self.relsite.taus = 3.500000
            self.relsite.ks = 0.600000
            self.relsite.kf = 4.311995
            self.relsite.tauf = 0.014630
            self.relsite.dD = 3.326148
            self.relsite.dF = 0.725512
            self.relsite.glu = 1.000000

        if select is 3:
            print "using IPSC#3 "
            self.relsite.F = 0.29594
            self.relsite.k0 = 0.44388 / 1000.0
            self.relsite.kmax = 15.11385 / 1000.0
            self.relsite.taud = 0.00260
            self.relsite.kd = 0.00090
            self.relsite.taus = 11.40577
            self.relsite.ks = 27.98783
            self.relsite.kf = 30.00000
            self.relsite.tauf = 0.29853
            self.relsite.dD = 3.70000
            self.relsite.dF = 2.71163
            self.relsite.glu = 4.97494



class PSD(object):
    """
    Base class for postsynaptic density mechanisms, possibly including cleft. 
    May accept either NetCon or pointer inputs from a Terminal, and directly
    modifies the membrane potential and/or ion concentrations of the 
    postsynaptic cell.    
    
    """

    def __init__(self, parent_section, target_section, terminal,
                    psdtype='ampa', message=None, debug=False,
                    thresh=0.0, gmax=1.0, gvar=0, eRev=0,
                    delay=0,
                    nmda_ratio=1.0, identifier=0):
        """ This routine generates the synaptic connections from one presynaptic
            input onto a postsynaptic cell.
            Each connection is a stochastic presynaptic synapse ("presynaptic") with
            depression and facilitation.
            Each synapse can have  multiple releasesites "releasesites" on the
            target cell, as set by "NRZones".
            Each release site releases transmitter using "cleftXmtr"
            Each release site's transmitter is in turn attached to a PSD at each ending ("psd")
            Each psd can have a different conductance centered about the mean of
            gmax, according to a gaussian distribution set by gvar.
        """
        parent_cell = cells.cell_from_section(parent_section)
        target_cell = cells.cell_from_section(target_section)
        
        #if cellname not in ['bushy', 'MNTB', 'stellate']:
            #raise TypeError
            #exit()
        #if debug:
            #print "\nTarget cell  = %s, psdtype = %s" % (cellname, psdtype)
        #printParams(stochastic_pars)
        glyslowPoMax = 0.162297  # thse were measured from the kinetic models in Synapses.py, as peak open P for the glycine receptors
        glyfastPoMax = 0.038475  # also later verified, same numbers...
        if psdtype == 'glyfast':
            gmax /= glyfastPoMax  # normalized to maximum open probability for this receptor
        if psdtype == 'glyslow':
            gmax /= glyslowPoMax  # normalized to max open prob for the slow receptor.
        
        self.netcons = [] # build list of connections from individual release sites to the mother calyx
        #        print "Stochastic syn: j = %d of n_fibers = %d n_rzones = %d\n" % (j, n_fibers, n_rzones)
        relzone = terminal.relsite
        n_rzones = terminal.n_rzones
        
        #
        # Create cleft mechanisms
        # 
        clefts = []
        if not runQuiet:
            print "adding %d clefts to terminal " % n_rzones, parent_section
        for k in range(0, n_rzones):
            cl = h.cleftXmtr(0.5, sec=target_section)
            clefts.append(cl) # cleft
        
        # and then make a set of postsynaptic zones on the postsynaptic side
        #        print 'PSDTYPE: ', psdtype
        if psdtype == 'ampa':
            (psd, psdn, par, parn) = template_iGluR_PSD(target_section, nReceptors=n_rzones,
                                                        nmda_ratio=nmda_ratio)
        elif psdtype == 'glyslow':
            (psd, par) = template_Gly_PSD_State_Gly6S(target_section, nReceptors=n_rzones,
                                                        psdtype=psdtype)
        elif psdtype == 'glyfast':
            (psd, par) = template_Gly_PSD_State_PL(target_section, nReceptors=n_rzones,
                                                    psdtype=psdtype)
        elif psdtype == 'glyGC':
            (psd, par) = template_Gly_PSD_State_GC(target_section, nReceptors=n_rzones,
                                                    psdtype=psdtype)
        elif psdtype == 'glya5':
            (psd, par) = template_Gly_PSD_State_Glya5(target_section, nReceptors=n_rzones,
                                                        psdtype=psdtype)
        elif psdtype == 'glyexp':
            (psd, par) = template_Gly_PSD_exp(target_section, nReceptors=n_rzones,
                                                psdtype=psdtype)
        else:
            print "**PSDTYPE IS NOT RECOGNIZED: [%s]\n" % (psdtype)
            exit()
            # connect the mechanisms on the presynaptic side together
        if debug:
            print 'parent_section: ', parent_section
        
        parent_section.push()

        self.netcons.append(h.NetCon(parent_section(0.5)._ref_v, relzone, thresh, delay, 1.0))
        self.netcons[-1].weight[0] = 1
        self.netcons[-1].threshold = -30.0

        for k in range(0, n_rzones): # 2. connect each release site to the mother axon
            if psdtype == 'ampa': # direct connection, no "cleft"
                relzone.setpointer(relzone._ref_XMTR[k], 'XMTR', psd[k])
                relzone.setpointer(relzone._ref_XMTR[k], 'XMTR', psdn[k]) # include NMDAR's as well at same release site
            elif psdtype == 'glyslow' or psdtype == 'glyfast' or psdtype == 'glyexp' or psdtype == 'glya5' or psdtype == 'glyGC':
                relzone.setpointer(relzone._ref_XMTR[k], 'pre',
                                    clefts[k]) # connect the cleft to the release of transmitter
                clefts[k].preThresh = 0.1
                if isinstance(target_cell, cells.TStellate):
                    clefts[k].KV = 531.0 # set cleft transmitter kinetic parameters
                    clefts[k].KU = 4.17
                    clefts[k].XMax = 0.731
                elif isinstance(target_cell, cells.Bushy):
                    clefts[k].KV = 1e9 # really fast for bushy cells.
                    clefts[k].KU = 4.46
                    clefts[k].XMax = 0.733
                else:
                    print('Error in cell name, dying:   '), cellname
                    exit()
                clefts[k].setpointer(clefts[k]._ref_CXmtr, 'XMTR', psd[k]) #connect transmitter release to the PSD
            else:
                print "PSDTYPE IS NOT RECOGNIZED: [%s]\n" % (psdtype)
                exit()
            v = 1.0 + gvar * np.random.standard_normal()
            psd[k].gmax = gmax * v # add a little variability - gvar is CV of amplitudes
            psd[k].Erev = eRev # set the reversal potential
            if psdtype == 'ampa': # also adjust the nmda receptors at the same synapse
                psdn[k].gmax = gmax * v
                psdn[k].Erev = eRev
        h.pop_section()
        par = list(par)
        if psdtype == 'ampa': # include nmda receptors in psd
            psd.extend(psdn)
            par.extend(parn)
        if not runQuiet:
            if message is not None:
                print message
                
        self.psd = psd
        self.clefts = clefts
        self.par = par




#
# Currently unused functions:
#

#def print_model_parameters(psd, par, message=None):
    #if message is not None:
        #print "%s" % message
    #for i in range(0, len(par)):
        #for j in range(0, len(psd[i])):
            #print '%s: ' % (psd[i].hname()),
            #for k in par[i]:
                #v = eval('psd[i][j].' + k)
                #print '%8s  = %12f' % (k, v)


#def ParameterInit(delay=0.6, latency=0.0, std=0.0):
    #""" convenience function to set up a default parameter group for stochastic synapses, in which
    #there is no time variance in release or release latency
    #Note the following: these parameters are passed to coh4.mod. 
    #In that routine, if LN_Flag is 0, then there is no time evolution of the distribution,
    #BUT the sigma value is passed (called "std") and is used to jitter the release events.
    #If you don't want this, set LN_std to 0.0
    #"""
    #vPars = Params(LN_Flag=0, LN_t0=10.0, LN_A0=0., LN_tau=10.0, LN_std=std,
                   #Lat_Flag=0, Lat_t0=10.0, Lat_A0=0., Lat_tau=10.,
                   #delay=delay, latency=latency)
    #return (vPars)


#def printParams(p):
    #"""
    #print the parameter block created in Parameter Init"""
    #u = dir(p)
    #print "-------- Synapse Parameter Block----------"
    #for k in u:
        #if k[0] != '_':
            #print "%15s = " % (k), eval('p.%s' % k)
    #print "-------- ---------------------- ----------"


