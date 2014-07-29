TITLE passive  (leak) membrane channel

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

NEURON {
	THREADSAFE
	SUFFIX leak
	NONSPECIFIC_CURRENT i
	RANGE gbar, erev
}

PARAMETER {
	v (mV)
	gbar = 0.001	(mho/cm2)
	erev = -65	(mV)
}

ASSIGNED { i	(mA/cm2)}

BREAKPOINT {
	i = gbar*(v - erev)
}


