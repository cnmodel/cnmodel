TITLE detailed model of Glycine receptors

COMMENT
-----------------------------------------------------------------------------

    Kinetic model of Glycine-A receptors
    ====================================



    C -- C1 -- C2
         |     |
         O1    O2

-----------------------------------------------------------------------------


    This mod file does not include mechanisms for the release and time course
    of transmitter; it is to be used in conjunction with a sepearate mechanism
    to describe the release of transmitter and that provides the concentration
    of transmitter in the synaptic cleft (to be connected to pointer C here).

-----------------------------------------------------------------------------

    Based on models

    Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  Kinetic models of
    synaptic transmission.  In: Methods in Neuronal Modeling (2nd edition;
    edited by Koch, C. and Segev, I.), MIT press, Cambridge, 1998, pp. 1-25.

    (electronic copy available at http://cns.iaf.cnrs-gif.fr)

    Written by Alain Destexhe, Laval University, 1995

-----------------------------------------------------------------------------

    Modified Paul Manis, UNC Chapel Hill, 2009
    Changed name, pointer name, kinetics are range variables, and kinetic values
    are estimated from VCN glycine receptors.

    This model does not have a desensitization state.

-----------------------------------------------------------------------------

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS GLYa5
	POINTER XMTR
	RANGE C0, C1, C2, O1, O2, Open
	RANGE g, gmax, f1, f2
	RANGE Erev
    RANGE  kf1, kf2, kb1, kb2, a1, b1, a2, b2
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

: from fits to averaged ipsc data, stellate cells 1/1/10

:    kf1	= 0.002930   (/uM /ms)	: binding
:	kf2	= 0.005936   (/uM /ms)	: binding
:	kb1	= 2.793	(/ms)	: unbinding
:	kb2	= 1.445	(/ms)	: unbinding
:	a1	= 1e-6	(/ms)	: opening
:	b1	= 129.0	(/ms)	: closing
:	a2	= 5.10	(/ms)	: opening
:	b2	= 2.79  (/ms)	: closing

: from fits to averaged ipsc data, bushy cells 1/1/10

	kf1	= 0.0278   (/uM /ms)	: binding
	kf2	= 1e-6   (/uM /ms)	: binding
	kb1	= 0.000054	(/ms)	: unbinding
	kb2	= 0.000855	(/ms)	: unbinding
	a1	= 1e-6	(/ms)	: opening
	b1	= 129.0	(/ms)	: closing
	a2	= 5.10	(/ms)	: opening
	b2	= 2.79  (/ms)	: closing

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
	O1		: open
	O2		: open
}

INITIAL {
	C0 = 1
	C1 = 0
	C2 = 0
	O1 = 0
	O2 = 0
    XMTR = 0.0
}

BREAKPOINT {
	SOLVE kstates METHOD sparse
    Open = (O1 + O2)
	g = gmax * Open
	i = (1e-6) * g * (v - Erev)
}

KINETIC kstates {

	f1 = kf1 * (1e3) * XMTR
	f2 = kf2 * (1e3) * XMTR

	~ C0 <-> C1	(f1,kb1)
	~ C1 <-> C2	(f2,kb2)
	~ C1 <-> O1	(a1,b1)
	~ C2 <-> O2	(a2,b2)

	CONSERVE C0+C1+C2+O1+O2 = 1
}
