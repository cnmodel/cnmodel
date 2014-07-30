: HH Low TEA-sensitive Purkinje potassium current
: Created 8/7/02 - nwg

NEURON {
	THREADSAFE
    SUFFIX kpkj2
	USEION k READ ek WRITE ik
	RANGE gbar, ik, gk
	GLOBAL ninf, ntau
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER {
	v		(mV)
	gbar = .002	(mho/cm2)
	
	nivh = -24	(mV)
	nik = 20.4 (mV)
	
	ek (mV)
}

ASSIGNED {
	gk (mho/cm2)
    ik (mA/cm2)
	ninf (1)
	ntau		(ms)
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
	gk = gbar * n^4
    ik = gk * (v - ek)
}

DERIVATIVE states {
	rates(v)
	n' = (ninf - n) / ntau
}

PROCEDURE rates(Vm (mV)) {
	LOCAL v
	v = Vm + 11	: Account for Junction Potential
	ninf = 1/(1+exp(-(v-nivh)/nik))
	ntau = 1000 * ntau_func(v)
}

FUNCTION ntau_func(v (mV)) (ms) {
	if (v < -20) {
		ntau_func = 0.000688 (ms) + 1 (ms)/(exp((v+64.2 (mV))/6.5 (mV))+exp((v-141.5 (mV))/-34.8 (mV)))
	} else {
		ntau_func = 0.00016 (ms) + 0.0008 (ms) *exp(-0.0267 * v /(1 (mV)))
	}
}
