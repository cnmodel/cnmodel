COMMENT

Iclamp2.mod: Electrode current injection, revised for our own use
I specifically needed more step changes in a single waveform than 
were provided by the default code. This routine gives us 5 levels:
holding, 3 steps, and a final (holding?) level.

This version generates a pulse of duration at amplitude I for
onset times in the vector onset[STEP] and durations in dur[STEP]
Each pulse is independent and pulses _could_ overlap.

Since this is an electrode current, positive values of I depolarize the cell
and in the presence of the extracellular mechanism there will be a change
in vext since I is not a transmembrane current but a current injected
directly to the inside of the cell.

(modified and borrowed extensively from other Neuron code, 
2001-2002. Paul B. Manis)

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

DEFINE NSTEP 5     : maximum number of steps supported in this mechanism 

NEURON {
    POINT_PROCESS IClamp2
    RANGE  onset, dur, amp, i
    ELECTRODE_CURRENT i
}
UNITS {
    (nA) = (nanoamp)
}
PARAMETER {
}

ASSIGNED {
	i (nA)
	index
	onset[NSTEP] (ms)
	dur[NSTEP] (ms)
	amp[NSTEP] (nA)
	pu
}

INITIAL {
	i = 0 (nA)
	index = 0
	pu = 0
}
:
: at each onset time, a pulse of duration dur is generated. 
:

BREAKPOINT {
	i = 0
	FROM index = 0 TO NSTEP-1 {
    		if(t > onset[index]) {
			if(t < (dur[index]+onset[index])) {
		        	i = i + amp[index] : allows overlap
			} : end of if condition
		} :end of second if condition
	} : end of "FROM" (for) loop
} : end of breakpoint
