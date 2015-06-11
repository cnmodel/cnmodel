TITLE Multisite synapse

COMMENT
-----------------------------------------------------------------------------
Multi-site synapse with independent release sites. Each site operates independently
and releases a vesicle upon presynaptic depolarization with a probability 
determined by the history of activity, using the Dittman and Regehr (1998, 2000)
model.

Revised from coh2.mod, coh3.mod, and coh4.mod.
The Dittman and Regeher (1998, 2000) release model with
facilitation closely fits  Auditory Nerve data from mouse over
a wide range of frequencies.
The model DOES NOT include the postsynaptic receptors or desensitization, since
these should be treated separately (couple XMTR to an AMPA receptor model,
such as the Trussell-Raman model)

Range variables:
nZones: is the number of active zones simulated in this calyx model. Each zone
        can be connected to a separate PSD.
nVesicles: the size of the release pool in the presynaptic terminal
F (0.4): The base release probability
k0 (1/1.75): /s, baseline recovery rate from depletion (slow rate)
kmax (1/0.025): /s, maximal recovery rate from depletion (fast rate)
td (0.05) : time constant for fast calcium-dependent recovery, sec
kd (0.7) : affinity of fast recovery process for calcium sensor
kf (0.5) : affinity of facilitation process
tf (0.01) : rate of facilitation process (slow) seconds
dD (0.02): calcium that drives recovery (ca influx per AP)
dF (0.02): calcium that drives facilitiation

Added latency and variable delay (latstd, latency standard deviation in msec)
around the mean spike time. 4/5/2011 pbm.

Version 4 uses a log-normal distribution to determine release latencies. 
The calcuation is built-in instead of being passed through an array. 
The lognormal distribution describes the individual vesicle release time
course at this synapse as measured by Isaacson and Walmsley, 1996. Note that
they used a gamma distribution in some plots, but the lognormal distribution 
seems to fit their published data at least as well. 
The parameters of the distribution, as well as the release latency,
are controlled by an exponential function whose parameters are initialized at
run time. 
10/19/2011 Paul B. Manis, UNC Chapel Hill

ENDCOMMENT

DEFINE MAX_ZONES 100  : maximum number of zones in this model
DEFINE EVENT_N 10000   : number of entries in the Event Distribution (e.g., as sampled)

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
THREADSAFE
    POINT_PROCESS MultiSiteSynapse
    RANGE F, k0, kmax, taud, kd, tauf, kf, taus, ks
    RANGE nZones, multisite, nVesicles, rseed, latency, latstd, debug
    RANGE dD, dF, XMTR, glu, CaDi, CaFi
    RANGE Fn, Dn
    RANGE TTotal
    RANGE nRequests, nReleases
    RANGE Identifier : just a number so we can report which instance is active
    RANGE TDur, TAmp
    : Distributions for stochastic release and testing (Sept, Oct, 2011):
    RANGE EventLatencies, EventTime : returns the first EVENT_N latencies and absolute times at which they were used
    RANGE ev_index : count in the EventLatencies (in case we are "short")
    : parameters for latency shift during repetitive stimulation (Oct 19, 2011)
    RANGE Dep_Flag : Depression flag (0 to remove depression; 1 to allow DKR control of facilitation and depression)
    RANGE Lat_Flag, Lat_t0, Lat_A0, Lat_tau : Lat_Flag  = 0 means fixed latency (set by "latency" above)
                                            : otherwise, latency = latency for t < Lat_t0
                                            :            latency = latency + Lat_A0*(1-exp(-(t-Lat_t0)/Lat_tau))
    : parameters for lognorm distribution shift during repetitive stimulation (Oct 19, 2011)
    RANGE LN_Flag, LN_t0, LN_A0, LN_tau : LN_Flag  = 0 means fixed sigma as well
                                            : otherwise, sigma = latstd for t < LN_t0
                                            :            sigma = latstd + LN_A0*(1-exp(-(t-LN_t0)/LN_tau))

    : externally assigned pointers to RNG functions
    POINTER binomial_rng  : for deciding the number of active synapses when multisite==0
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (mM) = (milli/liter)
    (uM) = (micro/liter)
}

PARAMETER {
    dt      (ms)
    n0 = 1      (1)        : initial size of RRVP at each release site
    TAmp = 1.0 (mM)          : amplitude of transmitter pulse
    TDur = 0.5    (ms)    : duration of transmitter pulse
    dD = 0.02 (1)       : calcium influx driving recovery per AP
    dF = 0.02 (1)       : calcium influx driving facilitation per AP
    F  = 0.5 (1)        : basal facilitaiton
    k0 = 0.0005714(/ms)       : slow recovery from depletion (1.0/1.75)
    kmax = 0.040 (/ms)      : fast recovery from depletion (1/0.025)
    taud = 50.0 (ms)        : time constant for fast calcium dependent recovery
    kd = 0.7 (1)       : affinity of fast recovery process for calcium sensor
    tauf = 10.0 (ms)        : rate of slow facilitiation process
    kf = 0.5 (1)       : affinity of slow facilitation process
    taus = 1 (ms)
    ks = 0.5 (1)      : not used but defined anyway
    glu = 1 (mM)
    rseed (1)        : random number generator seed (for SCOP module)
    latency = 0.0 (ms)
    latstd = 0.0 (ms)
    : Time course of latency shift in release during repetitive stimulation
    Lat_Flag = 0 (1) : 0 means fixed latency, 1 means lognormal distribution
    Lat_t0 = 0.0 (ms) : minimum time since simulation start before changes in latency are calculated
    Lat_A0 = 0.0 (ms) : size of latency shift from t0 to infinity
    Lat_tau = 100.0 (ms) : rate of change of latency shift (from fit of a+b(1-exp(-t/tau)))
    : Statistical control of log-normal release shape over time during repetive stimulation
    LN_Flag = 0 (1) : 0 means fixed values for all time
    LN_t0 = 0.0 (ms) : : minimum time since simulation start before changes in distribution are calculated
    LN_A0 = 0.0 (ms) : size of change in sigma from t0 to infinity
    LN_tau = 100.0 (ms) : rate of change of sigma over time (from fit of a+b*(1-exp(-t/tau)))
    
    : control flags     - if debug is 1, show all, if 2, just "some"
    debug = 0
    Identifier = 0
    Dep_Flag = 1 (1) : 1 means use depression calculations; 0 means always set release probability to F
}

ASSIGNED {
    : Externally set assignments
    nZones (1)    : number of zones in the model
    multisite (1)  : whether zones are modeled individualy (1) or as a single, variable-amplitude zone (0)
    nRequests (1) 
    nReleases (1)
    EventLatencies[EVENT_N] (0)
    EventTime[EVENT_N] (0)

    : Internal calculated variables
    Fn (1)
    Dn (1)
    CaDn (1)
    CaFn (1)
    CaDi (1)
    CaFi (1)
    eta (1)
    tSpike (ms)            : time of last spike
    tstep(ms)
    tRelease[MAX_ZONES] (ms)    : time of last release
    ZoneLatency[MAX_ZONES] (ms) : latency to release in this zone
    TTotal(0)
    tz(1)
    tspike (ms)
    latzone (ms)
    vesicleLatency (ms)
    sigma (ms)
    gindex (0)
    ev_index (0)
    scrand (0)    
    
    binomial_rng
}

: Function prototypes needed to assign RNG function pointers
VERBATIM
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM

: Return a pick from binomial distribution
FUNCTION rand_binomial() {
VERBATIM
    _lrand_binomial = nrn_random_pick(_p_binomial_rng);
ENDVERBATIM
}

: Function to allow RNG to be externaly set
PROCEDURE setBinomialRNG() {
VERBATIM
 {
    void** pv = (void**)(&_p_binomial_rng);
    *pv = nrn_random_arg(1);
 }
ENDVERBATIM
}


STATE {
    nVesicles[MAX_ZONES]    (1) : vesicles in RRVP
    XMTR[MAX_ZONES]       (mM)       : pulse of neurotransmitter
    TCount[MAX_ZONES] (1)  : timesteps remaining before zone becomes inactive
}

INITIAL {

    TTotal = 0
    nRequests = 0
    nReleases = 0
    set_seed(rseed)
:    VERBATIM
:        fprintf(stdout, "MultiSiteSynapse: Calyx #%d Initialized with Random Seed: %d\n", (int)Identifier, (int)rseed);
:    ENDVERBATIM
    nVesicles[0] = n0
    XMTR[0] = 0
    tRelease[0] = 0
    tSpike = -1000.0
    latzone = 0.0
    sigma = 0.0
    vesicleLatency = 0.0
    gindex = 0
    ev_index = 0
    scrand = 0.0
    CaDi = 1.0
    CaFi = 0.0
    CaDn = 1.0
    CaFn = 0.0
    Fn = F
    Dn = 1.0
    FROM i = 0 TO (nZones-1) {
        nVesicles[i] = n0
        XMTR[i] = 0
        tRelease[i] = 0
    }
    update(t-tSpike)
}

BREAKPOINT {
    SOLVE release
}

PROCEDURE release() {
    : Control the release process and transmitter concentration. 
    : The syanpse can release one vesicle per AP per zone, with a probabiliyt 0<p<1.
    : The probability, p, is defined by the time evolution of a Dittman-Regher model
    : of release, whose parameters are set during initialization.
    : The vesicle can be released over a variable time interval defined by a lognormal
    : distribution, plus a fixed latency.
    : Once released, the transmitter packet has a defined smooth time course in the "cleft"
    : represented by the product of rising and falling exponentials.
    tz = 0
    : update glutamate in cleft
    FROM i = 0 TO (nZones-1) { : for each zone in the synapse
        if (TCount[i] > 0 && t >= tRelease[i]) {
            TCount[i] = TCount[i] - 1
            if (TCount[i] < 0) { : count down and done
                XMTR[i] = 0.0
                TTotal = 0
                TCount[i] = 0
            }     : end of transmitter release pulse
            else {
                tz = t-(tRelease[i]+ZoneLatency[i]) : time since onset of release
                if (tz < 0) {
                    XMTR[i] = 0
                }
                else {
                    : calculate glutamate waveform
                    XMTR[i] = TAmp * (1.0-exp(-tz/(TDur/3.0))) * exp(-(tz-(TDur/3.0))/TDur) 
                } : done updating transmitter concentration
            }
        }
    }
}


PROCEDURE update(tstep (ms)) {
    : update the facilitation and depletion variables
    : from the Dittman-Regehr model. 
    : Updates are done with each new presynaptic AP event.
    if(tstep > 0.0) {
        CaDi = CaDi + dD
        CaFi = CaFi + dF
        CaDn = CaDi * exp (-tstep/taud)
        CaFn = CaFi * exp (-tstep/tauf)
        eta = (kd/CaDi + 1.0)/(kd/CaDi + exp(-tstep/taud))
        eta = eta^(-(kmax-k0)*taud)
        Dn = 1.0-(1.0-(1.0-Fn)*Dn)*exp(-k0*tstep)*eta
        Fn = F + (1.0-F)/(1.0+kf/CaFn)
        CaDi = CaDn
        CaFi = CaFn
    }
    if (Dep_Flag == 0) { : no depression
        Dn = 1.0  : set to 1
        Fn = F   : set to initial value to set release probability constant
    }
:    VERBATIM
:        if (debug >= 2 ){
:            fprintf(stdout, "update start t = %f ts=%f: F=%7.2f CaDi = %g CaFi = %g\n", \
:            t, tstep, F, CaDi, CaFi);
:            fprintf(stdout, "    vars:  taud=%g: tauf=%g kd = %g kmax= %g\n", taud, tauf, kd, kmax);
:            fprintf(stdout, "    CaDi = %g CaFi = %g\n", CaDi, CaFi);
:            fprintf(stdout, "    CaDn = %g CaFn = %g\n", CaDn, CaFn);
:            fprintf(stdout, "    eta: %g\n", eta);
:            fprintf(stdout, "    Fn=%7.2f Dn: %7.2f  CaDi = %g CaFi = %g,\n", \
:            Fn, Dn, CaDi, CaFi);
:        }
:    ENDVERBATIM
}


NET_RECEIVE(weight) {
    : Connect to here when a spike occurs...
    update(t - tSpike) : see if we need to do an update on the event
    tSpike = t      : save the time of spike
    TTotal = 0 : reset total transmitter from this calyx for each release
:    VERBATIM
:      if (debug == 1) {
:        fprintf(stderr, "  ---> Spike at t = %9.3f\n", tSpike);
:      }
:    ENDVERBATIM
    nRequests = nRequests + 1 : count the number of inputs that we received
    
    latzone = 0.0
    FROM i = 0 TO (nZones-1) { : for each zone in the synapse
        
        : now handle vesicle release...
        if(tRelease[i] < tSpike) {
            scrand = scop_random()
            : look to make release if we have not already (single vesicle per zone per spike)
            : check for release and release probability - assume infinite supply of vesicles
            if (scrand  < Fn*Dn) { 
                nReleases = nReleases + 1 : count number of releases since inception
                tRelease[i] = tSpike     : time of release
                TTotal = TTotal + 1     : count total releases this trial.

:               Compute the median latency for this vesicle.
                if (Lat_Flag == 0 || t < Lat_t0) {
                    vesicleLatency = latency : use a fixed value
                }
                else {
                    vesicleLatency = latency + Lat_A0*(1-exp(-(t-Lat_t0)/Lat_tau)) : latency rises during train
                }
:               Now compute distribution around that latency
:               if LN_Flag is 1, we adjust the sigma values over time, otherwise we just use a fixed value.
:               The following math applies:
:               lognormal dist = exp(X), where X = gaussian(mu, sigma).
:               The median of the lognormal dist. is e^mu; Note that if mu is 0, then the median is 1.0
:               The mode of the lognormal dist is e^(u-sigma^2). 
:               Note that normrand is a SCoP function, see http://cns.iaf.cnrs-gif.fr/files/scopman.html
                if (LN_Flag == 0 || t < LN_t0) { : use fixed sigma in lognormal distribution for all time.
                    if (latstd > 0.0) {
                        latzone = normrand(0.0, latstd) : set a latency for the zone with one draw from the distribution
                        latzone = exp(latzone) - 1.0 + vesicleLatency : the latency should not be too short.... 
                    }
                    else {
                        latzone = vesicleLatency : fixed value
                    }
                }
                else {
                    sigma = latstd + LN_A0*(1-exp(-(t-LN_t0)/LN_tau)) : time-dependent std shift
                    latzone = normrand(0.0, sigma)
                    latzone = exp(latzone)-1.0 + vesicleLatency
                }
                if (latzone < 0.0) { : this is to be safe... must have causality.
                    latzone = 0.0
                }
                if (ev_index < EVENT_N) { : save event distribution list for verification
                    EventLatencies[ev_index] = latzone
                    EventTime[ev_index] = t
                    ev_index = ev_index + 1
                }
                TCount[i] = (latzone / dt) + (5.0*TDur / dt) :
                ZoneLatency[i] = latzone : save the current zone latency... 
            }
        }
    }
}
