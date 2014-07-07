TITLE Calyx of Held Version 3

COMMENT
-----------------------------------------------------------------------------
Revised from coh2.mod.
This version uses the Dittman and Regeher (1998, 2000) release model with
facilitation, which can be closely fit to Auditory Nerve data from mouse over
a range of frequencies.
The model DOES NOT include the postsynaptic receptors or desensitization, since
these should be treated separately (couple XMTR to an AMPA receptor model,
such as the Trussell model)

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

added latency and variable delay (latstd, latency standard deviation in msec)
around the mean spike time. 4/5/2011 pbm.

ENDCOMMENT

DEFINE MAX_ZONES 1000  : maximum number of zones in this model
DEFINE GAMMA_N 10000   : number of entries in the gamma distribution
DEFINE EVENT_N 1000		: number of entries in the Event Distribution (e.g., as sampled)

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS COH3
	RANGE F, k0, kmax, taud, kd, tauf, kf, taus, ks
	RANGE nZones, nVesicles, rseed, latency, latstd, debug
	RANGE dD, dF, XMTR, glu
	RANGE Fn, Dn
	RANGE TTotal
	RANGE nRequests, nReleases
	RANGE Identifier : just a number so we can report which instance is active
	RANGE TDur, TAmp
	RANGE LatencyDist : input the release latency distribution here
	RANGE EventDist : returns the first 1000 events times used
	RANGE ScopDist : returns the first 1000 random numbers generated and checked.
	RANGE ev_index : count in the EventDist (in case we are "short")
	RANGE sc_index : count in the EventDist (in case we are "short")
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)
}

PARAMETER {
	dt	  (ms)
	n0 = 1	  (1)		: initial size of RRVP at each release site
	TAmp = 1.0 (mM)		  : amplitude of transmitter pulse
	TDur = 0.5	(ms)	: duration of transmitter pulse
	dD = 0.02 (1)	   : calcium influx driving recovery per AP
	dF = 0.02 (1)	   : calcium influx driving facilitation per AP
	F  = 0.5 (1)		: basal facilitaiton
	k0 = 0.0005714(/ms)	   : slow recovery from depletion (1.0/1.75)
	kmax = 0.040 (/ms)	  : fast recovery from depletion (1/0.025)
	taud = 50.0 (ms)		: time constant for fast calcium dependent recovery
	kd = 0.7 (1)	   : affinity of fast recovery process for calcium sensor
	tauf = 10.0 (ms)		: rate of slow facilitiation process
	kf = 0.5 (1)	   : affinity of slow facilitation process
	taus = 1 (ms)
	ks = 0.5 (1)	  : not used but defined anyway
	glu = 1 (mM)
    rseed (1)		: random number generator seed (for SCOP module)
	latency = 0.0 (ms)
	latstd = 0.0 (ms)
	: control flags	 - if debug is 1, show all, if 2, just "some"
  debug = 0
  Identifier = 0
}

ASSIGNED {
	: Externally set assignments
	nZones (1)			: number of zones in the model
	nRequests (1) 
	nReleases (1)
	LatencyDist[GAMMA_N] (0)
	EventDist[EVENT_N] (0)
	ScopDist[EVENT_N] (0)

	: Internal calculated variables
	Fn (1)
	Dn (1)
	CaDn (1)
	CaFn (1)
	CaDi (1)
	CaFi (1)
	eta (1)
	tSpike (ms)			: time of last spike
	tstep(ms)
	tRelease[MAX_ZONES] (ms)	: time of last release
	relThisSpike[MAX_ZONES] (0)	: a flag that lets us limit  
								:release to one vesicle per site per spike.
	ZoneLatency[MAX_ZONES] (ms) : latency to release in this zone
	nTotal (1)			: total (or mean) RRVP size
	TTotal(0)
	hasReleased(0)
	tz(1)
	tspike(ms)
	latzone(ms)
	iZone (1)
	gindex (0)
	ev_index (0)
	sc_index(0)
	gammasize (GAMMA_N)
	scrand (0)
}

STATE {
	nVesicles[MAX_ZONES]	(1) : vesicles in RRVP
	XMTR[MAX_ZONES]	   (mM)	   : pulse of neurotransmitter
	TCount[MAX_ZONES] (1)
	ZoneActive[MAX_ZONES] (1)  : indicate when the zone is currently "active"
}

INITIAL {

	VERBATIM
		if(debug == 1) {
			fprintf(stdout, "Initial for COH3, instance %d\n", (int)Identifier);
			fprintf(stdout, "Latency: %f   Var: %f\n", latency, latstd);
		}
	//	fprintf(stdout, "Debug flag: %d for instance %d\n", (int)debug, (int)Identifier);
	ENDVERBATIM
	nTotal = 0
	TTotal = 0
	nRequests = 0
	nReleases = 0
	set_seed(rseed)
	nVesicles[0] = n0
	XMTR[0] = 0
	tRelease[0] = 0
	hasReleased = 0
	relThisSpike[0] = 0
	ZoneActive[0] = 0
	tSpike = -1000.0
	latzone = 0.0
	gindex = 0
	ev_index = 0
	sc_index = 0
	scrand = 0.0
	CaDi = 1.0
	CaFi = 0.0
	CaDn = 1.0
	CaFn = 0.0
	Fn = F
	Dn = 1.0
	iZone = 0
	FROM i = 0 TO (nZones-1) {
		VERBATIM
			if (debug == 1) {
		//		fprintf(stderr, " Initializing release for zone= %.0f of %.0f\n", iZone, nZones);
			}
		ENDVERBATIM
		nVesicles[i] = n0
		XMTR[i] = 0
		tRelease[i] = 0
		relThisSpike[i] = 0
		ZoneActive[i] = 0
		iZone = iZone + 1
	}
	update(t-tSpike)
}

BREAKPOINT {
	SOLVE release
}

PROCEDURE release() {

	nTotal = 0
	tz = 0
	latzone = 0.0
	gindex = 1
	VERBATIM
		if (debug == 1) {
// 			  fprintf(stderr, "** Entering release @ t=%f, nZones = %d\n", t, (int)nZones);
		}
	ENDVERBATIM

	iZone = 0
	FROM i = 0 TO (nZones-1) { : for each zone in the synapse
		VERBATIM
			if (debug == 1) {
/* 
				  fprintf(stderr, "      Checking Zones, zone= %.0f of %d, ", iZone, (int)nZones);
				  fprintf(stderr, " relthisspike = %.0f, trel: %.3f, tspike: %.3f\n", \
				  relThisSpike[(int)iZone], tRelease[(int)iZone], tSpike);
 */
			}
		ENDVERBATIM
		scrand = scop_random()
		if (sc_index < EVENT_N) {
					ScopDist[sc_index] = scrand
					sc_index = sc_index + 1
				}
		VERBATIM
			if (debug == 1) {
// 				fprintf(stderr, "        scop_random= %8.3f, fn*dn = %8.3f\n", scrand, Fn*Dn);
			}
		ENDVERBATIM

		: now handle vesicle release...
		if(relThisSpike[iZone] == 0 && (tRelease[iZone] < tSpike)) {
		: look to make release if we have not already (single vesicle per zone per spike)
			VERBATIM
			if (debug == 1) {
				fprintf(stderr, "Zone %.0f: Checking Probability\n", iZone);
			}
			ENDVERBATIM	
: check for release and release probability - assume infinite supply of vesicles
			if (scrand  < Fn*Dn) { 
				VERBATIM
				if (debug == 1) {
					fprintf(stderr, "   *** Zone %d RELEASES at %7.2f.  scrnd: %8.3f   fndn: %8.3f\n", \
					(int)iZone, tSpike+latzone, scrand, Fn*Dn);
				}
				ENDVERBATIM
				nReleases = nReleases + 1 : count number of releases since inception
				tRelease[iZone] = tSpike	 : time of release
				TTotal = TTotal + 1		: count total releases this trial.
				hasReleased = 0 : keep track of whether release has happened yet,variable latency
				relThisSpike[iZone] = 1	 : keep track of release on this spike
				ZoneActive[iZone] = 1 : set the zone as "active"
				VERBATIM
					gindex = (int)(10000.0*scop_random()); // select a value randomly
//					fprintf(stderr, "       gindex: %d", (int)gindex);
				ENDVERBATIM
				if (gindex <= 0) {
					gindex = 1
				}
				if (gindex > GAMMA_N) {
					gindex = GAMMA_N
				}
				latzone = LatencyDist[gindex]
:				latzone = normrand(latency, latstd) : set a latency for the zone 
:				: (normrand is a SCoP function, see http://cns.iaf.cnrs-gif.fr/files/scopman.html)
:				: latzone = latency + latzone*latstd
:				if (latzone < 0.0) {
:					latzone = 0.0
:				}
				VERBATIM
					if (debug == 1) {
						fprintf(stderr, "     gindex: %d zone: %d  latency: %g \n", \
						(int)gindex, (int)iZone, latzone);
					}
				ENDVERBATIM
				if (ev_index < EVENT_N) {
					EventDist[ev_index] = latzone
					ev_index = ev_index + 1
				}
				TCount[iZone] = (latzone / dt) + (5.0*TDur / dt) :
				ZoneLatency[iZone] = latzone : save the current zone latency... 
				XMTR[iZone] = XMTR[iZone] + TAmp : bump transmitter
			}
			else {
				: flag that we did not release, but looked anyway (failure of release).
				: Prevents additional releases
				relThisSpike[iZone] = 2 
				VERBATIM
					if (debug == 1) {
						fprintf(stderr, "  ### Zone %d: Looked but did not release - scrnd: %8.3f   fndn: %8.3f\n",\
						 (int) iZone, scrand, Fn*Dn);
					}
				ENDVERBATIM
			}
		}
		: regardless of whether or not a release was triggered, update glutamate in cleft

		if (ZoneActive[iZone] > 0) {
			if (t >= tRelease[iZone]) { : do not start until we pass the zone latency...
				TCount[iZone] = TCount[iZone] - 1
				if (TCount[iZone] < 0) { : count down and done
					XMTR[iZone] = 0.0
					ZoneActive[iZone] = 0
					TTotal = 0
					TCount[iZone] = 0
					VERBATIM
						if (debug == 1) {
						//	fprintf(stderr, " Reset XMTR\n");
							}
					ENDVERBATIM
				}	 : end of transmitter release pulse
				else {
					tz = t-(tRelease[iZone]+ZoneLatency[iZone]) : time since onset of release
					if (tz < 0) {
						XMTR[iZone] = 0
					}
					else {
						VERBATIM
						if (hasReleased == 0 && debug == 1) {
						//	fprintf(stderr,  "First latency is %8.4f\n", t);
						}
						ENDVERBATIM
						hasReleased = 1
						: calculate glutamate waveform
						XMTR[iZone] = TAmp * (1.0-exp(-tz/(TDur/3.0))) * exp(-(tz-(TDur/3.0))/TDur) 
					}
					VERBATIM
						if (debug == 1) {
						//	fprintf(stderr, " Updated XMTR\n");
						}
					ENDVERBATIM
				}
			}
		}
	iZone = iZone + 1
	}						   : end of Zone loop - all vesicle states have been updated
}


PROCEDURE update(tstep (ms)) {

: update the facilitation and depletion variables
	VERBATIM
		if (debug == 2 ){
		fprintf(stdout, "update start t = %f ts=%f: F=%7.2f CaDi = %g CaFi = %g\n", \
		t, tstep, F, CaDi, CaFi);
		fprintf(stdout, "	vars:  taud=%g: tauf=%g kd = %g kmax= %g\n", taud, tauf, kd, kmax);
		}
	ENDVERBATIM
	if(tstep < 0.0) {
		return(0)
	}
	CaDi = CaDi + dD
	CaFi = CaFi + dF
	VERBATIM
		if (debug >= 2 ){
			fprintf(stdout, "	 CaDi = %g CaFi = %g\n", CaDi, CaFi);
		}
	ENDVERBATIM
	CaDn = CaDi * exp (-tstep/taud)
	CaFn = CaFi * exp (-tstep/tauf)
	VERBATIM
		if (debug >= 2 ){
			fprintf(stdout, "	 CaDn = %g CaFn = %g\n", CaDn, CaFn);
		}
	ENDVERBATIM
	eta = (kd/CaDi + 1.0)/(kd/CaDi + exp(-tstep/taud))
	eta = eta^(-(kmax-k0)*taud)
	VERBATIM
		if (debug >= 2 ){
			fprintf(stdout, "eta: %g\n", eta);
		}
	ENDVERBATIM
	Dn = 1.0-(1.0-(1.0-Fn)*Dn)*exp(-k0*tstep)*eta
	Fn = F + (1.0-F)/(1.0+kf/CaFn)
	CaDi = CaDn
	CaFi = CaFn
	VERBATIM
		if (debug >= 2 ){
			fprintf(stdout, "	  Fn=%7.2f Dn: %7.2f  CaDi = %g CaFi = %g,\n", \
			Fn, Dn, CaDi, CaFi);
		}
	ENDVERBATIM
}


NET_RECEIVE(weight) {
: Connect to here when a spike occurs...
	update(t - tSpike) : see if we need to do an update on the event
	tSpike = t	  : save the time of spike
	TTotal = 0 : reset total transmitter from this calyx for each release
	FROM i = 0 TO nZones-1 {
		relThisSpike[i] = 0	 : reset release indicator for each zone
	}
	VERBATIM
	  if (debug == 1) {
	    fprintf(stderr, "  ---> Spike at t = %9.3f\n", tSpike);
	  }
	ENDVERBATIM
	nRequests = nRequests + 1 : count the number of inputs that we received
}
