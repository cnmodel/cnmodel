: Ih current
: Created 8/6/02 - nwg

NEURON {
	SUFFIX hpkj
	NONSPECIFIC_CURRENT i
	RANGE ghbar, eh
	GLOBAL ninf, ntau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	v	 	(mV)
	
	ghbar = .0001	(S/cm2)

	eh = -30	(mV)
}

ASSIGNED {
	i (mA/cm2)
	ninf
	ntau (ms)
}

STATE {
	n
}

INITIAL {
	rates(v)
	n = ninf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	i = ghbar*n*(v - eh)
}

DERIVATIVE states {
	rates(v)
	n' = (ninf - n)/ntau
}

PROCEDURE rates(v (mV)) {
	ninf = 1/(1+exp((v+90.1(mV))/9.9(mV)))
	ntau = (1000) * (0.19 (s) + 0.72 (s)*exp(-((v-(-81.5(mV)))/11.9(mV))^2))
}