COMMENT
cleftXmtr

This is simple state model that generates "cleft" transmitter, through
the following scheme:

A netreceive block receives the driving event. This forces XV (the vesicle
state) to be set to XMax to mimic the release of a vesicle.
Then:
XV --> XC --> XU
where XV is the vesicle transmitter, XC is the cleft transmitter and
XU is transmitter that has been taken up. The forward rates are finite, and the
reverse rates are 0 (XU is an absorbing state)

The forward rate kv1 mimics simple first-order diffusion across the cleft
The forward rate ku1 mimics simple first-order uptake from the cleft

The concentration XC is available to the program as Xmtr.
XMax is the max cleft concentration of transmitter.

Because vesicle release events at a single presynaptic terminal can be nearly
simultaneous, it is important that this mechanism does not have a refractory
period. We also assume that the uptake mechanism is not saturable.

Paul B. Manis, Ph.D.
UNC Chapel Hill
3 Jan 2010

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

DEFINE NSTEP 5

NEURON {
    POINT_PROCESS cleftXmtr
    POINTER pre
    RANGE  KV, KU, XMax
    RANGE CXmtr, preThresh
}
UNITS {
    (nA) = (nanoamp)
}

PARAMETER { : Parameters are chosen from best fit to stellate cell data in VCN
    KV = 531 (/ms)   <0,1e9> : release rate from vesicle
	KU = 4.17 (/ms) <0,1e3>  : uptake rate
    XMax = 0.731 (mM)
    preThresh = 0
}

ASSIGNED {
    pre
    CXmtr (mM)
    preLast (1)
    tLast
}

STATE {
	XV : Vesicle transmitter (just released)
    XC : Cleft transmitter (e.g., at receptor)
    XU : Uptake state (dead state... )
}

INITIAL {
    XV = 0
    XC = 0 (mM)
    XU = 0
    CXmtr = 0.0
    preLast = 0.0
    tLast = 0.0
}

BREAKPOINT {
    SOLVE kstates METHOD sparse
    CXmtr = XC*XMax
}

KINETIC kstates {
    ~ XV <-> XC	(KV, 0.0)
	~ XC <-> XU	(KU, 0.0)
    : note that this mechanism has no CONSERVATION : XU can accumulate as much
    : as needed.
}

NET_RECEIVE(conc (mM)) { : detect and cause a release event
	XV = XV + 1
}
