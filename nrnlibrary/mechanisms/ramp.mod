COMMENT
ramp stimulator
Since this is an electrode current, positive values of i depolarize the cell
and in the presence of the extracellular mechanism there will be a change
in vext since i is not a transmembrane current but a current injected
directly to the inside of the cell.

Ramp by P.Manis 5/11/99
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

DEFINE NSTEP 4

NEURON {
	POINT_PROCESS IRamp
	RANGE  del, dur, hold, amp, hold2
	ELECTRODE_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
}
PARAMETER {
	del (ms)	<0,1e9>
	dur (ms)	<0,1e9>
	hold (nA)
	amp (nA)
	hold2 (nA)
}
ASSIGNED {
i (nA)
}

INITIAL {
	i = 0
}

BREAKPOINT {
	if(t < del) {
		i = hold
	}
	if(t >= del && t < del+dur) {
		i = hold + (amp/dur)*(t-del)
	}
	if(t >= del+dur) {
		i = hold2
	}
}

