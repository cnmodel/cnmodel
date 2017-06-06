TITLE passive  (leak) membrane channel

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

NEURON {
	THREADSAFE
	SUFFIX leak
	NONSPECIFIC_CURRENT i
	RANGE gbar, erev, q10g
}

PARAMETER {
	v (mV)
	gbar = 0.001	(mho/cm2)
	erev = -65	(mV)
    q10g = 2.0
}

ASSIGNED { 
    i	(mA/cm2)
    qg ()
}

INITIAL {
    qg = q10g^((celsius-22)/10 (degC))
}
BREAKPOINT {
	i = qg*gbar*(v - erev)
}


