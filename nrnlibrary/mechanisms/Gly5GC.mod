TITLE detailed model of Glycine receptors

COMMENT
-----------------------------------------------------------------------------

    Kinetic model of Glycine-A receptors
    ====================================

    C -- C1 -- C2 -- O1
          |     |
         D1 -- D2 -- D3

-----------------------------------------------------------------------------

    This Model is based on:
    Gentet LJ, Clements JD Binding site stoichiometry and the effects of
    phosphorylation on human alpha1 homomeric glycine receptors J Physiol (Lond)
    2002 vol. 544 (Pt 1) pp. 97-106, Figure 7.

    Written by Paul Manis, UNC Chapel Hill, 2009
    Kinetic values are estimated from VCN glycine receptors.

    This model has desensitization states.

-----------------------------------------------------------------------------

    This mod file does not include mechanisms for the release and time course
    of transmitter; it is to be used in conjunction with a sepearate mechanism
    to describe the release of transmitter and that provides the concentration
    of transmitter in the synaptic cleft (to be connected to pointer C here).

-----------------------------------------------------------------------------

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS GLYaGC
	POINTER XMTR
	RANGE C0, C1, C2, D1, D2, D3, O1,  Open
	RANGE g, gmax, f1, f2
	RANGE Erev
    RANGE  k1, km1, a1, b1, d1, r1, d2, r2, d3, r3, rd, dd
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

: Rates
: bushy cell

	k1	= 12.81	(/uM /ms)	: binding
	km1 = 0.0087	(/ms)	: unbinding
	a1	= 0.0194	(/ms)	: opening
	b1	= 1.138	(/ms)	: closing
	r1	= 5.19	(/ms)	: desense 1
	d1	= 0.000462  (/ms)	: return from d1
    r2 = 0.731 (/ms)     : return from deep state
    d2 = 1.641 (/ms)     : going to deep state
    r3 = 3.817 (/ms)     : return from deep state
    d3 = 1.806 (/ms)     : going to deep state
    rd = 1.0 (/ms)
    dd = 1.0 (/ms)
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	g 		(pS)		: conductance
	XMTR 	(mM)		: pointer to glycine concentration

	f1		(/ms)    : binding
	f2		(/ms)    : binding
    Open    (1)
}

STATE {
	: Channel states (all fractions)
	C0		: unbound
	C1		: single bound
	C2		: double bound
    D1      : desense, bound
	O1		: open
	D2		: Desense
    D3      : Desense
}

INITIAL {
	XMTR = 0
    C0 = 1
	C1 = 0
	C2 = 0
	O1 = 0
	D1 = 0
    D2 = 0
    D3 = 0
}

BREAKPOINT {
	SOLVE kstates METHOD sparse
    Open = (O1)
	g = gmax * Open
	i = (1e-6) * g * (v - Erev)
}

KINETIC kstates {

	f1 = k1 * (1e3) * XMTR
	f2 = k1 * (1e3) * XMTR

	~ C0 <-> C1	(f1,km1)
	~ C1 <-> C2	(f2,12.5*km1)
	~ C2 <-> O1	(a1,b1)
    ~ C1 <-> D1 (r1, d1)
	~ C2 <-> D2	(r2, d2)
    ~ D1 <-> D2 (rd, dd)
    ~ D2 <-> D3 (r3, d3)

	CONSERVE C0+C1+C2+D1+D2+D3+O1 = 1
}
