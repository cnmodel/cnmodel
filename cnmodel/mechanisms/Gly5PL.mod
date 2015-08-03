TITLE detailed model of Glycine receptors

COMMENT
-----------------------------------------------------------------------------

    Kinetic model of Glycine-A receptors: Pascal Legendre (Mauthner Cell)
    ====================================



    C0--C1--C2--O1
            |
            C3--O2

-----------------------------------------------------------------------------


    This mod file does not include mechanisms for the release and time course
    of transmitter; it is to be used in conjunction with a sepearate mechanism
    to describe the release of transmitter and that provides the concentration
    of transmitter in the synaptic cleft (to be connected to pointer C here).


-----------------------------------------------------------------------------

    Modified Paul Manis, UNC Chapel Hill, 2009
    Name, pointer name, kinetics are range variables, and kinetic values
    are estimated from VCN glycine receptors.

    Note: This model does not have a desensitization state.

-----------------------------------------------------------------------------

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS GLYaPL
	POINTER XMTR
	RANGE C0, C1, C2, C3, O1, O2, Open
	RANGE g, gmax, f1, f2
	RANGE Erev
    RANGE  kon, koff, a1, b1, a2, b2, r, d
    RANGE CellType : 0 for bushy, 1 for stellate
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

	Erev	= -70  (mV)	: reversal potential
	gmax	= 500  (pS)	: maximal conductance
    CellType = 1 (1) : define cell type parameters

: Rates

: Stellate cell fit (1/1/10; excellent fit)

:	kon	= 0.0236	(/uM /ms)	: binding
:	koff = 2.4	(/ms)	: unbinding
:	a1	= 1.707 (/ms)	: opening
:	b1	= 8.95	(/ms)	: closing
:	a2	= 0.325	(/ms)	: opening
:	b2	= 5.871  (/ms)	: closing
:    r = 2.019 (/ms)     : return from deep state
:    d = 28.87 (/ms)     : going to deep state

:if psdtype == 'glyfast':  fit from 3/5/2010. error =  0.174 maxopen = 0.0385
: See synapses.py
    a1 =        1.000476 (/ms) : opening
    a2 =        0.137903 (/ms) : opening
    b1 =        1.700306 (/ms) : closing
    koff =     13.143132 (/ms) : unbinding
    kon =       0.038634 (/ms) : binding
    r =         0.842504 (/ms) : return from deep state
    b2 =        8.051435 (/ms) : closing
    d =        12.821820 (/ms) : going to deep state

}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	g 		(pS)		: conductance
	XMTR 	(mM)		: pointer to glycine concentration

	f1		(/ms)    : binding
	f2		(/ms)    : binding
    koff2   (/ms)
    Open    (1)
}

STATE {
	: Channel states (all fractions)
	C0		: unbound
	C1		: single bound
	C2		: double bound
    C3      : bound but closed state to O2
	O1		: open
	O2		: open
}

INITIAL {

    XMTR = 0.0
	C0 = 1
	C1 = 0
	C2 = 0
    C3 = 0
	O1 = 0
	O2 = 0
}

BREAKPOINT {
	SOLVE kstates METHOD sparse
:VERBATIM
:   if (CGly > 0.0) {
:      fprintf(stderr, "t = %f   Xmtr = %f\n", t, XMTR);
:   }
:   ENDVERBATIM
    Open = (O1 + O2)
	g = gmax * Open
	i = (1e-6) * g * (v - Erev)
}

KINETIC kstates {

	f1 = 2.0 * kon * (1e3) * XMTR
	f2 = kon * (1e3) * XMTR
    koff2 = 2.0 * koff

	~ C0 <-> C1	(f1,koff)
	~ C1 <-> C2	(f2,koff2)
	~ C2 <-> O1	(a1,b1)
    ~ C2 <-> C3 (d, r)
	~ C3 <-> O2	(a2,b2)

	CONSERVE C0+C1+C2+C3+O1+O2 = 1
}
