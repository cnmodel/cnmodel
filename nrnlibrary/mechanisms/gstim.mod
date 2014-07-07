COMMENT
gStim

This is a point conductance, Positive values of the amplitude depolarize the cell
and in the presence of the extracellular mechanism there will be a change
in vext since i is not a transmembrane current but a current injected
directly to the inside of the cell.

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

DEFINE NSTEP 5

NEURON {
    POINT_PROCESS gStim
    RANGE  dur, delay, g, Eg
    ELECTRODE_CURRENT i
}
UNITS {
    (nA) = (nanoamp)
	(mV) = (millivolt)
	(nS) = (nanomho)
}
PARAMETER {
    dur (ms)   <0,1e9>
	delay (ms) <0,1e9>
	g (nS)
	Eg (mV)
	v (mV)
}
ASSIGNED {
i (nA)
}

INITIAL {
    i = 0
}

BREAKPOINT {
    if(t < delay) {
        i = 0
    }
	if(t >= delay && t <= (delay+dur)) {
		i = -(0.001) * g * (v - Eg)
	}
}
