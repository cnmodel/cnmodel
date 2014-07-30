: HH Slow TEA-insensitive Purkinje potassium current
: Created 8/7/02 - nwg

NEURON {
	THREADSAFE
    SUFFIX kpkjslow
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
	gbar = 0.004	(mho/cm2)
	
	nivh = -16.5	(mV)
	nik = 18.4 (mV)
	
	ek (mV)
}

ASSIGNED {
    gk (mho/cm2)
	ik (mA/cm2)
	ninf
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

FUNCTION ntau_func(v (mV)) (ms){
	ntau_func = 0.000796 (ms) + 1 (ms)/(exp((v+73.2 (mV))/11.7 (mV))+exp((v-306.7 (mV))/-74.2(mV)))
}