TITLE Calyx of Held

COMMENT
-----------------------------------------------------------------------------

 Model of vesicle mobilization and release at multiple
 release sites in the calyx of Held.

 - basic enhanced replenishment model

 - each release site produces a pulse of transmitter, T, when
   a vesicle is released
    - T used as the pointer to the AMPA receptor

B. Graham, Dept. of Computing Science & Maths, University of Stirling
(Contact: b.graham@cs.stir.ac.uk)
(previously IANC, Division of Informatics, University of Edinburgh)

CNS 2000 Version (19/11/02)
-----------------------------------------------------------------------------


P. Manis: modifications: (Sept and October 2004, January 2005).
RSIZE is static allocation variable; nzones is actual number of zones used
pV0 is different for each instance of the calyx.
Modified to limit vesicular release to a single event
at each site for each AP (prior code didn't do this).

ENDCOMMENT

DEFINE SSIZE 10001
DEFINE RSIZE 1000  : maximum number of zones in this model

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    POINT_PROCESS COH
    RANGE CGLU, ntot, Ttot, spike, PR0, nzones, Camp, rseed
    GLOBAL KC, KO, R
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (mM) = (milli/liter)
}

PARAMETER {
    dt    (ms)
    : rseed = 0   : random number seed
    n0 = 1    (1)   : initial size of RRVP at each release site
    ke = 80.0 (/mM /ms)   : 80 rate of enhanced replenishment
    kd = 0.002  (/ms)  : 0.002 background depletion of rrvp
    Cres = 0 (mM)   : amplitude of local residual [Ca]
    : a shift from 22 to 33 deg c at a q10 of 3 yields a rate differential of 3.35 CDur and Cndur are adjusted from the original
	: values by this much.
	Cdur = 0.25 (ms)   : duration of local [Ca] transient (was 1)
    Cnamp = 0.01 (mM)   : amplitude of distant [Ca] transient
    Cnres = 0 (mM)  : amplitude of distant residual [Ca]
    Cndur = 0.5 (ms)  : duration of  distant [Ca] transient (was 2)
    Tamp = 1 (mM)   : amplitude of transmitter pulse
    Tdur = 0.1  (ms)  : 1 duration of transmitter pulse (Modified, 9/5/04 PBM to 0.5; and 3/15/05 to 0.1)
}

ASSIGNED {
    spike[SSIZE]    (ms)    : list of spike times
    rseed  : random number seed to start ending with (usually, run number)
	PR0		 : iniitial release prob for this ending
    nzones  (1) : number of zones in the model
    index           : index to spike times
    tspike      (ms)    : time of last spike
    trel[RSIZE] (ms)    : time of last release
    relthisspike[RSIZE] (0)   : block multiple release per spike coming in.
	ntot        (1) : total (or mean) RRVP size
    Ttot        (mM)    : number of releases
    Camp (mM) : amplitude of local [Ca] transient following AP
    						: default used was 0.1; 0.05 gives less depression consistent
    						: with VCN data

    km          (/ms)   : background replenishment rate
    R   (1)      : probability of vesicle release
    KO[2]       (/mM /ms) : gate opening rates
    KC[2]       (/ms)     : gate closing rates
    inf[2] tau[2] fac[2]
	up (1) : pv (just for viewing/testing)
	ui (1) : index holder...
}

STATE {
    n[RSIZE]    (1) : vesicles in RRVP
	Rmax (1) : max value of R (at peak Camp)
    RO[2]   (1) : gates for release
    CGLU[RSIZE]    (mM)    : pulse of neurotransmitter
    C       (mM)    : [Ca] release transient
    Cn      (mM)    : [Ca] mobilization transient
    pv[RSIZE]   (1) : scaling for release probability
    Ccount (1)  : count of time steps for Ca transient (BPG 10-1-02)
    Cncount (1)  : count of time steps for mob Ca transient (BPG 13-1-02)
    Tcnt[RSIZE] (1)  : count of time steps for T transient (BPG 13-1-02)
}

INITIAL {
    index = 0
       n[0] = n0
      CGLU[0] = 0
      trel[0] = 0
      relthisspike[0] = 0
	  pv[0] = PR0
    FROM i = 0 TO nzones-1 {
      n[i] = n0
      CGLU[i] = 0
      trel[i] = 0
      relthisspike[i] = 0
	  pv[i] = PR0
    }
    R = 0
    C = Cres
    Cn = Cnres
    ntot = 0
    Ttot = 0
    set_seed(rseed)
:  VERBATIM
:   fprintf(stdout, "COH: seed set to %f, PR=0 = %f\n", rseed, pv[0]);
:	ENDVERBATIM
	km = n0 * kd        : background replenishment rate
    KO[0] = 150 : fast gate opening
    KO[1] = 1 : slow gate opening
    KC[0] = 30 : fast gate closing
    KC[1] = 0.1 : slow gate closing

	C = Camp
	prel() : now compute R (fraction of gates open at peak calcium)
	Rmax = 0.000275953 : magic number - should have been R, but testing reveals calculated value is too low for P to work correctly.
	: Note: this value only works for Cres = 0.05, which is the VCN endbulb model (PBM 4/2/05).
	: The value is obtained by uncommenting the VERBATIM code below in one_rel(i), and reading R This is the kinetic
	: parameter at the moment of a spike (after C is incremented.).
	: Values were obtained for a single, isolated input spike.
	:
	C = Cres : return to resting calcium
	}

BREAKPOINT {
    SOLVE release
}

PROCEDURE release() {

: detect presynaptic spike and reset calcium/times  and release.
: Logic changed so that t must be within dt of spike time for this to happen
: old logic was just > spike[index].

:	if (index < SSIZE && t>=spike[index] && C < Camp) {    : presynaptic spike
	if ((index < SSIZE) && ((t-spike[index])>=0) && ((t-spike[index]) < dt) && (C < Camp)) {    : presynaptic spike
			C = Camp    : calcium transient
			index=index+1 : index of next spike
			tspike = t    : time of spike
			Ccount = Cdur / dt  : number of time steps for Ca transient
			Ttot = 0
			FROM i = 0 TO nzones-1 {
				relthisspike[i] = 0  : reset release indicator for each zone
				}
	}

	prel()    : probability of release

	ntot = 0
	FROM i = 0 TO nzones-1 {
			if (unirand() < dt*km) {
				n[i] = n[i]+1
			}  : background replenishment
			if (n[i] > 0 && unirand() < dt*kd) {
				n[i] = n[i]-1
			}  :depletion
			if (Cn > 0) {
				if (unirand() < dt*ke*Cn) {
					n[i] = n[i]+1
				} : extra replenishment
			}

	 	  onerel(i)   : single vesicle release

			if (CGLU[i] > 0) {
				Tcnt[i] = Tcnt[i] - 1
				if (Tcnt[i] < 0) {CGLU[i] = 0}    : end of pulse
			}

			ntot = ntot + n[i] : total of available vesicles
	} 								: end of Zone loop - all vesicle states have been updated

	if (C > Cres) {
			Ccount = Ccount - 1
			if (Ccount < 0) {
					C = Cres  : end of release transient
					Cn = Cnamp  : beginning of mobilization transient
					Cncount = Cndur / dt  : number of time steps for Ca transient
			}
	}

	if (Cn > Cnres) {
			Cncount = Cncount - 1
			if (Cncount < 0) {Cn = Cnres}
	}

	VERBATIM
	return 0;
	ENDVERBATIM
}


PROCEDURE onerel(i) {   : one release per site per AP
		if(relthisspike[i] == 0) {  : look to make release if we have not already
			if (trel[i] < tspike) {
				: probability is scaled so that pv yields true release probability for unaffected spikes.
					ui = R*pv[i]*n[i]/Rmax
:					up = n[i]
:					VERBATIM
:					fprintf(stdout, ">>>>>>>>>>>>>>>>>>>>>>>>>n = %f   R = %g ui = %f C = %f\n", up, R, ui, C);
:					ENDVERBATIM
					if (n[i] > 0 && unirand() < ui*dt*40) { : scaling (31-1-01)
							n[i] = n[i] - 1     : release
							CGLU[i] = Tamp     : pulse of transmitter
							Tcnt[i] = Tdur / dt     : duration of pulse of transmitter
							trel[i] = t     : time of release
							Ttot = Ttot + 1     : count total releases (BPG 10-1-00)
							relthisspike[i] = 1  : keep track of release on this spike
					}
					else {
						relthisspike[i] = 2 : flag that we did not release, but looked anyway (failure of release)
					}
			}
		}
}


PROCEDURE prel() {  : Probability of release
  rates(C)
  R = 1
  FROM i=0 TO 1 {
    RO[i] = RO[i] + fac[i]*(inf[i] - RO[i])
    R = R * RO[i]
  }
}


PROCEDURE rates(C) {LOCAL a, b  :Computes gate rates at concentration C.
        TABLE inf, fac, tau DEPEND dt FROM 0 TO 0.2 WITH 200
        FROM j=0 TO 1 {
            a = KO[j] * C
            b = KC[j]
            tau[j] = 1/(a + b)
            inf[j] = a/(a + b)
            fac[j] = (1 - exp(-dt/tau[j]))
        }
}


FUNCTION unirand() {    : uniform random numbers between 0 and 1
        return(scop_random())
}
