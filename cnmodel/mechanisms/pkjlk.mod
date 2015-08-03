TITLE Purkinje Leak Current

: A passive purkinje cell leak current
NEURON {
	SUFFIX lkpkj
	NONSPECIFIC_CURRENT i
	RANGE i, e, gbar
}
PARAMETER {
	gbar = 5e-5	(siemens/cm2)  < 0, 1e9 >
	e = -60.5	(millivolt)
}
ASSIGNED {
	i  (milliamp/cm2)
	v  (millivolt)
}
BREAKPOINT { i = gbar*(v - e) }
