TITLE simple GABAa receptors

COMMENT
-----------------------------------------------------------------------------

	Simple model for GABAa receptors
	================================

  - FIRST-ORDER KINETICS, FIT TO WHOLE-CELL RECORDINGS

    Whole-cell recorded GABA-A postsynaptic currents (Otis et al, J. Physiol. 
    463: 391-407, 1993) were used to estimate the parameters of the present
    model; the fit was performed using a simplex algorithm (see Destexhe et
    al., J. Neurophysiol. 72: 803-818, 1994).

  - SHORT PULSES OF TRANSMITTER (0.3 ms, 0.5 mM)

    The simplified model was obtained from a detailed synaptic model that 
    included the release of transmitter in adjacent terminals, its lateral 
    diffusion and uptake, and its binding on postsynaptic receptors (Destexhe
    and Sejnowski, 1995).  Short pulses of transmitter with first-order
    kinetics were found to be the best fast alternative to represent the more
    detailed models.

  - ANALYTIC EXPRESSION

    The first-order model can be solved analytically, leading to a very fast
    mechanism for simulating synapses, since no differential equation must be
    solved (see references below).



References

   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  An efficient method for
   computing synaptic conductances based on a kinetic model of receptor binding
   Neural Computation 6: 10-14, 1994.  

   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J. Synthesis of models for
   excitable membranes, synaptic transmission and neuromodulation using a 
   common kinetic formalism, Journal of Computational Neuroscience 1: 
   195-230, 1994.

-----------------------------------------------------------------------------
ENDCOMMENT



INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS GABAa_S
	POINTER pre
	RANGE C, R, R0, R1, g, gmax,  Erev, lastrelease
	NONSPECIFIC_CURRENT i
	GLOBAL Cmax, Cdur, Alpha, Beta, Prethresh, Deadtime, Rinf, Rtau
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	Cmax	= 0.5	(mM)		: max transmitter concentration
	Cdur	= 0.3	(ms)		: transmitter duration (rising phase)
	Alpha	= 10.5	(/ms mM)	: forward (binding) rate
	Beta	= 0.166	(/ms)		: backward (unbinding) rate
	Erev	= -80	(mV)		: reversal potential
	Prethresh = 0 			: voltage level nec for release
	Deadtime = 1	(ms)		: mimimum time between release events
	gmax		(umho)		: maximum conductance
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	g 		(umho)		: conductance
	C		(mM)		: transmitter concentration
	R				: fraction of open channels
	R0				: open channels at start of release
	R1				: open channels at end of release
	Rinf				: steady state channels open
	Rtau		(ms)		: time constant of channel binding
	pre 				: pointer to presynaptic variable
	lastrelease	(ms)		: time of last spike
}

INITIAL {
	R = 0
	C = 0
	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
	Rtau = 1 / ((Alpha * Cmax) + Beta)
	lastrelease = -1000
}

BREAKPOINT {
	SOLVE release
	g = gmax * R
	i = g*(v - Erev)
}

PROCEDURE release() { LOCAL q
	:will crash if user hasn't set pre with the connect statement 

	q = ((t - lastrelease) - Cdur)		: time since last release ended

						: ready for another release?
	if (q > Deadtime) {
		if (pre > Prethresh) {		: spike occured?
			C = Cmax			: start new release
			R0 = R
			lastrelease = t
		}
						
	} else if (q < 0) {			: still releasing?
	
		: do nothing
	
	} else if (C == Cmax) {			: in dead time after release
		R1 = R
		C = 0.
	}



	if (C > 0) {				: transmitter being released?

	   R = Rinf + (R0 - Rinf) * exptable (- (t - lastrelease) / Rtau)
				
	} else {				: no release occuring

  	   R = R1 * exptable (- Beta * (t - (lastrelease + Cdur)))
	}

	VERBATIM
	return 0;
	ENDVERBATIM
}

FUNCTION exptable(x) { 
	TABLE  FROM -10 TO 10 WITH 2000

	if ((x > -10) && (x < 10)) {
		exptable = exp(x)
	} else {
		exptable = 0.
	}
}
