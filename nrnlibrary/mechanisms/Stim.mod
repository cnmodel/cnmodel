COMMENT
Since this is an electrode current, positive values of i depolarize the cell
and in the presence of the extracellular mechanism there will be a change
in vext since i is not a transmembrane current but a current injected
directly to the inside of the cell.
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

DEFINE NSTEP 5

NEURON {
    POINT_PROCESS IClamp2
    RANGE  dur1, dur2, dur3, dur4, dur5, amp1, amp2, amp3, amp4,
    amp5
    ELECTRODE_CURRENT i
}
UNITS {
    (nA) = (nanoamp)
}
PARAMETER {
    dur1 (ms)   <0,1e9>
    dur2 (ms)   <0,1e9>
    dur3 (ms)   <0,1e9>
    dur4 (ms)   <0,1e9>
    dur5 (ms)   <0,1e9>
    amp1 (nA)
    amp2 (nA)
    amp3 (nA)
    amp4 (nA)
    amp5 (nA)
}
ASSIGNED {
i (nA)
}

INITIAL {
    i = 0
}

BREAKPOINT {
    if(t < dur1) {
        i = amp1
    }
    if(t < (dur1+dur2) && t >= dur1) {
        i = amp2
    }
    if(t < (dur1+dur2+dur3) && t >= (dur1+dur2)) {
        i = amp3
    }
    if(t < (dur1+dur2+dur3+dur4) && t >= (dur1+dur2+dur3)) {
        i = amp4
    }
    if(t < (dur1+dur2+dur3+dur4+dur5) && t >= (dur1+dur2+dur3+dur4)) {
            i = amp5
        }
    if(t > (dur1+dur2+dur3+dur4+dur5)) {
        i = 0
    }
}
