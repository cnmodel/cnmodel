TITLE Model of glycine receptors

COMMENT
-----------------------------------------------------------------------------

    Kinetic model of glycine receptors
    ===============================

    6-state gating model with
    2 open states provide dual exponential response.

                O1
                |
    C0 -- C1 -- C2 -- O2
                |
                D1

-----------------------------------------------------------------------------

    This mod file does not include mechanisms for the release and time course
    of transmitter; it is to be used in conjunction with a sepearate mechanism
    to describe the release of transmitter and that provides the concentration
    of transmitter in the synaptic cleft (to be connected to pointer C here).

    Default parameters are set for a miniature EPSC.


-----------------------------------------------------------------------------

Paul B. Manis, Ph.D. 28 Dec 2009

-----------------------------------------------------------------------------
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    POINT_PROCESS Gly6S
    POINTER XMTR
    RANGE C0, C1, C2, D1, O1, O2, Open
    RANGE Erev
    RANGE Rb, Ru1, Ru2, Rd, Rr, Ro1, Rc1, Ro2, Rc2
    RANGE g, rb, gmax
    RANGE CellType : 0 is bushy, 1 is stellate
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (umho) = (micromho)
    (mM) = (milli/liter)
    (uM) = (micro/liter)
}

PARAMETER {

    Erev     = -70    (mV) : reversal potential
    gmax     = 500  (pS)  : maximal conductance
    CellType = 0 (1)   : Cell type definition.
: Rates
: Bushy cell IPSCs: (rates can be changed externally)
: Set per Fits 8 March 2010 (see Synapses.py as well)
    Rd  = 1.177999   (/ms)  : desensitization
    Rr  = 0.000005 (/ms)    : resensitization
    Rb  = 0.009403   (/mM /ms): binding
                : diffusion limited
    Ru2 = 0.000086  (/ms)    : unbinding (2nd site)
    Ro1 = 0.187858   (/ms)  : opening (fast)
    Ro2 = 1.064426    (/ms)   : opening (slow)
    Ru1 = 0.028696  (/ms)    : unbinding (1st site)
    Rc1 = 0.103625    (/ms)  : closing
    Rc2 = 1.730578    (/ms) : closing

}

ASSIGNED {
    v       (mV)    : postsynaptic voltage
    i       (nA)    : current = g*(v - Erev)
    g       (pS)    : conductance
    XMTR    (mM)    : pointer to glutamate concentration
    rbind   (/ms)   : binding
    Open    (1)
}

STATE {
    : Channel states (all fractions)
    C0      : unbound
    C1      : single gly bound
    C2      : double gly bound
    D1       : double gly bound, desensitized
    O1      : double gly bound, open state 1
    O2      : double gly bound, open state 2
}

INITIAL {

    XMTR = 0
	C0 = 1
    C1 = 0
    C2 = 0
    D1  = 0
    O1 = 0
    O2 = 0
}

BREAKPOINT {
    SOLVE kstates METHOD sparse
:VERBATIM
:   if (XMTR > 0.1) {
:      fprintf(stderr, "t = %f   XMTR = %f\n", t, XMTR);
:   }
:   ENDVERBATIM
    Open = (O1 + O2)
    g = gmax * Open
	i = (1e-6)*g*(v-Erev)
}


KINETIC kstates {

    rbind = Rb * (1e3) * XMTR

    ~ C0 <-> C1 (rbind,Ru1)
    ~ C1 <-> C2 (rbind*2.0,Ru2)
    ~ C2 <-> O1 (Ro1,Rc1)
    ~ C2 <-> D1  (Rd,Rr)
    ~ C2 <-> O2 (Ro2,Rc2)
    CONSERVE C0+C1+C2+D1+O1+O2 = 1
}
