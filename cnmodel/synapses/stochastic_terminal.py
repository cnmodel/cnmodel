from __future__ import print_function
from neuron import h

from .terminal import Terminal
from ..util import random_seed

# utility class to create parameter lists... 
# create like: p = Params(abc=2.0, defg = 3.0, lunch='sandwich')
# reference like p.abc, p.defg, etc.
class Params(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


class StochasticTerminal(Terminal):
    """
    Axon terminal with multi-site sctochastic release mechanism.
    """
    def __init__(self, pre_sec, target_cell, nzones=1, multisite=True, 
                 message=None, type='lognormal', identifier=0,
                 stochastic_pars=None, calcium_pars=None, delay=0., debug=False,
                 select=None, spike_source=None, spike_section=None, dep_flag=1):
        """
        This routine creates a (potentially) multisite synapse using a NEURON mod file with:
            - A MultiSiteSynapse release mechanism that includes stochastic release, with a lognormal
              release latency distribution. The Facilitation and Depression of release are governed
              by parameters obtaine from fitting the Dittman-Kreitzer-Regher (DKR) model (J Neurosci. 2000 Feb 15;20(4):1374-85.)
              to experimental data at various
              frequencies.
            - A "cleft" mechanism (models diffusion of transmitter). Note that the cleft is inserted as part of the
              presynaptic section, but is not connected to the postsynaptic side yet.
                 
        Each release site is stochastically independent, but all sites within a terminal are drive by the
        same presynaptic action potentials.
        
        Turning off the depression and facilitation in the kinetic portion of the model substantially decreases
        the time the terminal mechanism takes to run. 
                 
        Parameters
        ----------
        pre_sec : :obj:`section`
             The NEURON section where the synaptic mechanisms should be inserted.
        target_cell : :obj:`Cell`
             The target cell object that the synapse will innervate. 
        nzones : int
              The number of activate zones to insert into the section.
        multisite : bool, default: True
            A flag that determines whether the terminal actually creates multiple 
            release zones (True) or just creates a single release zone that
            varies its amplitude based on the depression/facilitation state.
        message : str
              A message to   when instantiating (mostly for verification of code flow).
        type: str (default: 'lognormal')
            'lognormal' sets the release event latency distribution to use a lognormal function. Currently,
            no other function is supported.
        identifier : int (default: 0)
             An identifier to associate with these release sites so we can find them later.
        stochastic_pars : dict (default: None)
             A dictionary of parameters (Param class) used to specifiy the stochastic behavior of this site,
             including release latency, stdev, and lognormal distribution paramaters
        calcium_pars : dict (default: None)
             A dictionary of parameters (Param class) to determine the calcium channels in this section.
             If None, then no calcium channels are inserted; otherwise, a P-type calcium conductance and a dynamic
             mechanism are inserted, and their conductance is set.
        delay : float (default: 0)
             Delay time in msec from action potential until transmitter release for this terminal.
        debug : bool (default: False)
             Flag to print stuff out while debugging. 
        spike_source : :obj:`section` (default: None)
             The input spike source to use in net con - default is to use pre_sec when set to None.
        dep_flag : int (default: 1)
             Set to 1 for depression mechanism (slows computation), 0 to turn off the depression calculations

        Returns
        -------
        list
            the list contains the terminal, the relsites, and the list of cleft mechanisms:
                 
            - terminal: this is the pointer to the terminal section that was inserted (same as pre_sec if it was
              specified)
            - relsite: a list of the nzones release sites that were created
            - cleft: a list of the nzones cleft mechanisms that were created.
        
        """
        Terminal.__init__(self, pre_sec)
        
        
        # set parameter control for the stochastic release of vesicles...
        # this structure is passed to stochastic synapses, and replaces several variables 
        # that were previously defined in the call to that function.
            
        thresh = -30 # mV - AP detection on the presynaptic side.
        
        ANTerminals_Latency = 0.5 # latency 
        vPars = Params(LN_Flag=1, LN_t0=10.0, LN_A0=0.05, LN_tau=35, LN_std=0.05,
                    Lat_Flag=1, Lat_t0=10.0, Lat_A0=0.140, Lat_tau=21.5,
                    latency=ANTerminals_Latency)
        #NOTE: stochastic_pars must define parameters used by multisite, including:
            #.delay is the netcon delay between the presynaptic AP and the start of release events
            #.Latency is the latency to the mean release event... this could be confusing.
        
        if stochastic_pars is None:
            stochastic_pars = vPars
            
        message='  >> creating StochasticTerminal with %d release zones using lognormal release latencies' % nzones
        if debug:
            print(message)
        terminal = pre_sec
        #terminal.push()
        if calcium_pars is not None:
            terminal.insert('cap') # insert calcium channel density
            terminal().cap.pcabar = calcium_pars.Ca_gbar
            terminal.insert('cad')
            
        # Create point process to simulate multiple independent release zones.
        relsite = h.MultiSiteSynapse(0.5, sec=terminal)
        relsite.nZones = nzones
        if multisite:
            relsite.multisite = 1
            relsite.rseed = random_seed.current_seed()  # use global random seed
            relsite.latency = stochastic_pars.latency
            relsite.latstd = stochastic_pars.LN_std
            self.n_rzones = nzones
        else:
            relsite.multisite = 0
            self.release_rng = h.Random(random_seed.current_seed())
            self.release_rng.uniform(0, 1)
            relsite.setUniformRNG(self.release_rng)
            self.n_rzones = 1
        
        relsite.Dep_Flag = dep_flag  # control synaptic dynamics
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
            #mpl.figure(2)
            
        h.pop_section()
        self.relsite = relsite
        sourcesec = pre_sec
        if spike_section is not None:
            sourcesec = spike_section
        if spike_source is None:
            spike_source = sourcesec(0.5)._ref_v
        if spike_source == 'cai':
            spike_source = sourcesec(0.5)._ref_cai
        self.netcon = h.NetCon(spike_source, relsite, thresh, delay, 1.0, sec=sourcesec)
        self.netcon.weight[0] = 1
        self.netcon.threshold = -30.0

        self.setPsdType(target_cell, select)

    def setPsdType(self, target_cell, select=None):
        """
        Assign a postsynpatic density type - selection of receptors to be
        associated with the presynaptic terminal.
        
        Parameters
        ----------
        target_cell : :obj:`NEURON section`
            Define the target cell so that the correct PSD is inserted
        
        select : str
            Not used in current implementation.
        
        """
        # TODO: must resurrect this for inhibitory synapses.
        # (and move constants out to their respective synapse locations)
        #elif psdtype.startswith('gly'):
            #self.setDF(target_cell, 'ipsc', select) # set the parameters for release
        self.setDF(target_cell, 'epsc') # set the parameters for release
        

    ################################################################################
    # The following routines set the synapse dynamics, based on measurements and fit
    # to the DKR model.
    ################################################################################

    def setDF(self, target_cell, synapsetype, select=None):
        """
        Set the facilitation/depression parameters for the multisite release model.
        The parameters used here were obtained from an optimized fit of the DKR
        model to stimulus and recovery data for the auditory nerve synapses onto bushy
        and T-stellate cells at 100, 200 and 300 Hz (simultaneous fitting), for times out to about
        0.5 - 1.0 second. Data were collected by Dr. Ruili Xie and Dr. Yong Wang.
        Fitting by Paul Manis.

        Parameters
        ----------
        target_cell : :obj:`NEURON section`
            Define the target cell so that the correct release kinetics are used
        
        synapsetype : str
            String defining the type of synapse: 'ipsc' or 'epsc'
        
        select : str (default: None)
            Select kinetcs from a particular example cell.
        
        """
        from .. import cells
        if isinstance(target_cell, cells.Bushy):
            if synapsetype == 'ipsc':
                if select is None:
                    self.bushy_ipsc_average()
                else:
                    self.bushy_ipsc_single(select=select)
        elif isinstance(target_cell, cells.TStellate):
            if synapsetype == 'ipsc':
                self.stellate_ipsc()

    def set_params(self, **params):
        """Set arbitrary parameters on the release mechanism.
        
        Parameters
        ----------
        \**params : dict
            dictionary of parameters to set on the release mechanism.
        
        """
        for k, v in params.items():
            setattr(self.relsite, k, v)

    def stellate_ipsc(self):
        """ Facilitation/Depression parameters for DKR model for IPSCs onto stellate cells.
        IPSCs were derived from stimulation of the dorsal cochlear nucleus, and so represent
        the glycinergic synapses from tuberculoventral neurons.
        Data is average of 3 cells studied with recovery curves and individually fit, at 100 Hz.
        """
        self.relsite.F = 0.23047
        self.relsite.k0 = 1.23636
        self.relsite.kmax = 45.34474
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
        """ Facilitation/Depression parameters for DKR model for IPSCs onto bushy cells.
        IPSCs were derived from stimulation of the dorsal cochlear nucleus, and so represent
        the glycinergic synapses from tuberculoventral neurons.
        Data for the kinetcs are the average of 16 bushy cells.
        The individual fits were compiled, and an average computed for just the 100 Hz data
        across the individual fits. This average was then fit to the DKR model
        (also see Matthew Xu-Friedman's papers). 
        The individual cells show a great deal of variability, from straight depression, to 
        mixed depression/facilitaiton, to facilation alone. This set of parameters generates
        a weak facilitation followed by depression back to baseline.
        
        There are other data sets in the source that are commented out.
        """
        #print( "USING average kinetics for Bushy IPSCs")

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


    def bushy_ipsc_single(self, select=None):
        """ 
        Facilitation/Depression parameters for DKR model for IPSCs onto 3 individual bushy cells.
        IPSCs were derived from stimulation of the dorsal cochlear nucleus, and so represent
        the glycinergic synapses from tuberculoventral neurons.
        
        Parameters
        ----------
        select : int, default=None
            select which cell's measured kinetics will be used:
        
            0: Use average 
            1: use data from 30aug08f
            2: select = 2, use data from 30aug08h
            3: select = 3, use data from 31aug08b (single cell, clean dataset)
        
        """
       # print ("Using bushy ipsc")

        if select is None or select > 4 or select <= 0:
            self.bushy_ipsc_average()
            return

        if select is 1: # 30aug08f
#            print ("using 30aug08f ipsc")
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
#            print ("using 30aug08H ipsc")
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
#            print ("using IPSC#3 ")
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
