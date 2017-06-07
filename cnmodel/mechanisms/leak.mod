TITLE passive  (leak) membrane channel

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

NEURON {
	THREADSAFE
	SUFFIX leak
	NONSPECIFIC_CURRENT i
	RANGE gbar, erev, i
}

PARAMETER {
	v (mV)
	gbar = 0.001	(mho/cm2)
	erev = -65	(mV)
}

ASSIGNED { 
    i	(mA/cm2)
}

INITIAL {
    
}
BREAKPOINT {
	i = gbar*(v - erev)
}


