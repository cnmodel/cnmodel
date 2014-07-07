TITLE Purkinje Leak Current

: A passive purkinje cell leak current
NEURON {
	SUFFIX lkpkj
	NONSPECIFIC_CURRENT i
	RANGE i, e, g
}
PARAMETER {
	g = 5e-5	(siemens/cm2)  < 0, 1e9 >
	e = -60.5	(millivolt)
}
ASSIGNED {
	i  (milliamp/cm2)
	v  (millivolt)
}
BREAKPOINT { i = g*(v - e) }
