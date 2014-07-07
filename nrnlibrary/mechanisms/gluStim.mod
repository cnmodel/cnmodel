COMMENT
gluStim

This is a point process pulse of glutamate to be attached to a receptor.

This is meant to be used with Vector Play, although it can also generate
just pulses of glutamate.
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

DEFINE NSTEP 5

NEURON {
    POINT_PROCESS gluStim
    RANGE  dur, delay, gluMax
    RANGE CGLU
}
UNITS {
    (nA) = (nanoamp)
}
PARAMETER {
    dur (ms)   <0,1e9>
	delay (ms) <0,1e9>
    gluMax (mM)
}
ASSIGNED {

}
STATE {
	CGLU   (mM)    : pulse of neurotransmitter
}

INITIAL {
    CGLU = 0 (mM)
}

BREAKPOINT {
COMMENT
if(t < delay || t > (delay+dur)) {
        CGLU = 0
    }
	if(t >= delay && t <= (delay+dur)) {
		CGLU = gluMax
	}
ENDCOMMENT
}
