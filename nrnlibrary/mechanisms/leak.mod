TITLE passive  (leak) membrane channel

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

NEURON {
	SUFFIX leak
	NONSPECIFIC_CURRENT i
	RANGE g, erev
}

PARAMETER {
	v (mV)
	g = .001	(mho/cm2)
	erev = -63.6	(mV)
}

ASSIGNED { i	(mA/cm2)}

INITIAL {
	i = g*(v-erev)
}

BREAKPOINT {
	i = g*(v - erev)
}


