COMMENT
-----------------------------------------------------------------------------

  Leak potassium current
  ----------------------

  This mechanism was written to be used as a potassium channel that is
  open or closed by neuromodulators.  

  WARNING: this current is NOT inserted as a standard current, but as a 
  point process (same way as a synapse or current injection).

  Procedure for insertion:

	objectvar kl
	kl = new kleak()

	access <compartment_name>
	kl.loc(0.5)

	kl.gmax = ...



  A. Destexhe , The Salk Institute, Feb 1994.

-----------------------------------------------------------------------------
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS kleak
	RANGE gmax
	GLOBAL Erev
	NONSPECIFIC_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	gmax	= 0.004	(umho)		: maximum conductance (microSiemens)
	Erev	= -100	(mV)		: reversal potential (potassium)
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
}

INITIAL {
}

BREAKPOINT {
	i = gmax * (v - Erev)
}

