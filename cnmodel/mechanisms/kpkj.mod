: HH TEA-sensitive Purkinje potassium current
: Created 8/5/02 - nwg

NEURON {
	THREADSAFE
    SUFFIX kpkj
	USEION k READ ek WRITE ik
	RANGE gbar, ik, gk
	GLOBAL minf, hinf, mtau, htau
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER {
	v (mV)

	gbar = .004	(mho/cm2)

	mivh = -24	(mV)
	mik = 15.4	(mV)
	mty0 = .00012851 (s)	
	mtvh1 = 100.7	(mV)
	mtk1 = 12.9	(mV)
	mtvh2 = -56.0	(mV)
	mtk2 = -23.1	(mV)
	
	hiy0 = .31	
	hiA = .78
	hivh = -5.802	(mV)
	hik = 11.2	(mV)

	ek (mV)
}

ASSIGNED {
	gk    (mho/cm2)
    ik		(mA/cm2)
	minf
	mtau		(ms)
	hinf
	htau		(ms)
}

STATE {
	m
	h
}

INITIAL {
	rates(v)
	m = minf
	h = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gk = gbar * m^3 * h
    ik =  gk * (v - ek)
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m) / mtau
	h' = (hinf - h) / htau
}

PROCEDURE rates( Vm (mV)) {
	LOCAL v
	v = Vm + 11 (mV)	: Account for Junction Potential
	minf = 1/(1+exp(-(v-mivh)/mik))
	mtau =  mtau_func(v)
	hinf = hiy0 + hiA/(1+exp((v-hivh)/hik))
	htau = 1000 * htau_func(v)
}

FUNCTION mtau_func (v (mV)) (ms) {
	if (v < -35 (mV)) {
		mtau_func = (1000)*(3.4225e-5+.00498*exp(-v/-28.29 (mV)))*3 (s)
	} else {
		mtau_func = (1000)*(mty0 + 1(s)/(exp((v+mtvh1)/mtk1)+exp((v+mtvh2)/mtk2)))
	}
}

FUNCTION htau_func(Vm (mV)) (ms) {
	if ( Vm > 0) {
		htau_func = (1000)*(0.0012 (s) + 0.0023(s)*exp(-0.141 *Vm / 1 (mV)))
	} else {
		htau_func = (1000)*(1.2202e-05(s) + .012(s) * exp(-((Vm-(-56.3 (mV)))/49.6 (mV))^2))
	}
}
	