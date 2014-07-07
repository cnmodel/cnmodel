TITLE Calyx of Held Version 2

COMMENT
-----------------------------------------------------------------------------

 Model of vesicle mobilization and release at multiple
 release sites in the calyx of Held.

 - basic enhanced replenishment model

 - each release site produces a pulse of transmitter, T, when
   a vesicle is released
	- T can then be used as the pointer to the AMPA receptor

B. Graham, Dept. of Computing Science & Maths, University of Stirling
(Contact: b.graham@cs.stir.ac.uk)
(previously IANC, Division of Informatics, University of Edinburgh)

CNS 2000 Version (19/11/02)
-----------------------------------------------------------------------------

P. Manis: modifications: (Sept and October 2004, January 2005).
MAX_ZONES is static allocation variable for the maximum # of release zones
nzones is actual number of zones used.
pV0 is different for each instance of the model
Modified to limit vesicular release to a single event
at each site for each AP (prior code didn't do this).

P. Manis, 11/2009:
Modified to use a NET RECEIVE block to detect input events.
This is a multi-site synapse driven by a single input.
Each zone in "nzones" is statistically independent.

ENDCOMMENT

DEFINE MAX_ZONES 1000  : maximum number of zones in this model

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS COH2
	RANGE ntot, Ttot, PR0, nzones, Camp, rseed, Nves
	GLOBAL KC, KO, R
	RANGE CGLU
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
	Cres = 0.05(mM)   : amplitude of local residual [Ca]. This value required to be 0.05 mm.
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
	rseed  : random number seed to start ending with (usually, run number)
	PR0		 : iniitial release prob for this ending
	nzones  (1) : number of zones in the model
	tspike      (ms)    : time of last spike
	trel[MAX_ZONES] (ms)    : time of last release
	relthisspike[MAX_ZONES] (0)   : a flag that lets us limit  release to one vesicle per site per spike.
	ntot        (1) : total (or mean) RRVP size
	Ttot        (mM)    : Total transmitter...
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
	NVes[MAX_ZONES]    (1) : vesicles in RRVP
	RScale (1) : Scaling factor for R (at peak Camp) so that release follows PV0
	RO[2]   (1) : gates for release
	CGLU[MAX_ZONES]    (mM)    : pulse of neurotransmitter
	C       (mM)    : [Ca] release transient
	Cn      (mM)    : [Ca] mobilization transient
	pv[MAX_ZONES]   (1) : scaling for release probability
	Ccount (1)  : count of time steps for Ca transient (BPG 10-1-02)
	Cncount (1)  : count of time steps for mob Ca transient (BPG 13-1-02)
	Tcnt[MAX_ZONES] (1)  : count of time steps for T transient (BPG 13-1-02)
}

INITIAL {
	NVes[0] = n0
	CGLU[0] = 0
	trel[0] = 0
	relthisspike[0] = 0
	tspike = 0.0
	pv[0] = PR0
	FROM i = 0 TO nzones-1 {
		NVes[i] = n0
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
: Initialize RO
	rates(C)
	FROM i = 0 TO 1 {
		RO[i] = inf[i]
	}
	prel() : now compute R (fraction of gates open at peak calcium)
	RScale = sqrt(R)
    C = Cres : return to resting calcium
}

BREAKPOINT {
	SOLVE release
}

PROCEDURE release() {
	prel()    : compute the probability of release at this moment
	ntot = 0
	FROM i = 0 TO nzones-1 { : for each zone in the synapse
		if (unirand() < dt*km) { : calculate background replenishment
			NVes[i] = NVes[i]+1
		}
		if (NVes[i] > 0 && unirand() < dt*kd) { : calculate vesicle depletion
			NVes[i] = NVes[i]-1
		}
		if (Cn > 0.0) { : amplitude of distant calcium controls extra replenishment
			if (unirand() < dt*ke*Cn) {
				NVes[i] = NVes[i]+1
			}
		}
		onerel(i)   : now figure out the single vesicle release
		if (CGLU[i] > 0) { : keep track of glutamate in cleft
			Tcnt[i] = Tcnt[i] - 1
			if (Tcnt[i] < 0) {
				CGLU[i] = 0
                Ttot = 0
            }    : end of transmitter release pulse
		}
		ntot = ntot + NVes[i] : total of available vesicles
	} 						   : end of Zone loop - all vesicle states have been updated
	if (C > Cres) { : compute decay of local calcium
		Ccount = Ccount - 1
		if (Ccount < 0) { : stepwise...
			C = Cres  : end of release transient
			Cn = Cnamp  : beginning of distant calcium transient
			Cncount = Cndur / dt  : number of time steps for distant Ca transient
		}
	}
	if (Cn > Cnres) { : compute decay of distant calcium
		Cncount = Cncount - 1
		if (Cncount < 0) { : stepwise decrease
			Cn = Cnres
		}
	}
	VERBATIM
	return 0;
	ENDVERBATIM
}

PROCEDURE onerel(izone) {   : one release per zone per AP. izone is the zone index
	if(relthisspike[izone] == 0) {  : look to make release if we have not already (single vesicle per zone per spike)
		if (trel[izone] < tspike) {
			ui = R*pv[izone]/RScale : probability is scaled so that pv yields true release probability for unaffected spikes.
			if (NVes[izone] > 0 && unirand() < ui) { : scaling ui*dt*40 (31-1-01)
				NVes[izone] = NVes[izone] - 1     : release decreases vesicle count at this zone
				CGLU[izone] = Tamp     : pulse the transmitter
				Tcnt[izone] = Tdur / dt     : duration of pulse of transmitter
				trel[izone] = t     : time of release
				Ttot = Ttot + 1     : count total releases (BPG 10-1-00)
				relthisspike[izone] = 1  : keep track of release on this spike
			}
			else {
				relthisspike[izone] = 2 : flag that we did not release, but looked anyway (failure of release). Prevents additional releases
			}
		}
	}
}

PROCEDURE prel() {  : calculate the probability of release
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

NET_RECEIVE(weight) {
: Connect to here when a spike occurs...
	if (C < Camp) {    : presynaptic spike
		C = Camp    : calcium transient
		tspike = t    : save the time of spike
		Ccount = Cdur / dt  : number of time steps for Ca transient
		Ttot = 0 : reset total transmitter from this calyx for each release
		FROM i = 0 TO nzones-1 {
			relthisspike[i] = 0  : reset release indicator for each zone
		}
	}
}
