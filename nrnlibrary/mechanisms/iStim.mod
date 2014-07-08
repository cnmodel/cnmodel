COMMENT
iStim

This is a point current injection (like an electrode).
Positive values of the amplitude depolarize the cell
and in the presence of the extracellular mechanism there will be a change
in vext since i is not a transmembrane current but a current injected
directly to the inside of the cell.

This is meant to be used with Vector Play...
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

DEFINE NSTEP 5

NEURON {
    THREADSAFE
    POINT_PROCESS iStim
    RANGE  dur, delay, iMax
    ELECTRODE_CURRENT i
}
UNITS {
    (nA) = (nanoamp)
}
PARAMETER {
    dur (ms)   <0,1e9>
	delay (ms) <0,1e9>
    iMax (nA)
}
ASSIGNED {
i (nA)
}

INITIAL {
    i = 0
}

BREAKPOINT {
COMMENT
if(t < delay || t > (delay+dur)) {
        i = 0
    }
	if(t >= delay && t <= (delay+dur)) {
		i = iMax
	}
ENDCOMMENT
}
