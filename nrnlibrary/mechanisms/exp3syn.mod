COMMENT
Three state kinetic scheme synapse described by rise time tau_rise,
and decay times constant tau1 and tau2, with tau_rise having a power integer exponent, n.
The normalized peak condunductance is 1. The ratio relative magnitude of tau1 and tau2
is given by F. A is a factor that makes the amplitude 1. Erev is the synapse reversal potential.

g = A*(1-exp(-t/taur))^n * ((1-F)*exp(-t/tau1) + F*exp(-t/tau2)))

ENDCOMMENT

NEURON {
	POINT_PROCESS Exp3Syn
	RANGE taur, tau1, tau2, N, F, Erev,  A
	NONSPECIFIC_CURRENT i

}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	taur = 0.1 (ms) <1e-9,1e9>
	tau1 = 3 (ms) <1e-9,1e9>
	tau2 = 20 (ms) <1e-9,1e9>
	Erev = 0	(mV)
	F = 0.8
    N = 2
    A = 1
}

ASSIGNED {
	v (mV)
	i (nA)
}

::::::::::PROGRAM NOT YET MODIFIED PAST THIS POINT ::::::::

STATE {
	g (uS)
}

INITIAL {
    g = 0
}

BREAKPOINT {
	SOLVE state METHOD cnexp

	g = A*pow((1-exp(-t/taur)),N) * ((1-F)*exp(-t/tau1) + F*exp(-t/tau2))
	i = g*(v - Erev)
}

DERIVATIVE state {
	g' = -g/tau1
}

NET_RECEIVE(weight (uS)) {
	A = A + weight
}
