COMMENT
alpha function synapse implemented as continuously integrated
kinetic scheme a la Srinivasan and Chiel (Neural Computation) so that
one can give many stimuli which summate.

Onset times are placed in the vector onset[SIZE]
Conductance located in state variable G
The amplitude of each individual alpha function is given by stim,
stim * t * exp(-t/tau).
The last onset time should be a very large number so stim stops getting
added to state A
ENDCOMMENT

DEFINE SIZE 1000

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}


NEURON {
	POINT_PROCESS KSyn100
	RANGE tau, stim, e, i,onset
	NONSPECIFIC_CURRENT i
}




UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	tau = 0.1 (ms)
	stim = 0.05 (umho)
	e=0	(mV)
	v	(mV)
	
}

ASSIGNED {
	index
	i (nA)
	bath (umho)
	k (/ms)
	onset[SIZE] (ms)
}

STATE {
	A (umho)
	G (umho)
}

INITIAL {
	k = 1/tau
	A = 0
	G = 0
	index=0
}

? current
BREAKPOINT {
	SOLVE conductance
	i = G*(v - e)
}

: at each onset time a fixed quantity of material is added to state A
: this material moves through G with the form of an alpha function

PROCEDURE conductance() { 
	LOCAL x
	while(index < SIZE && t>onset[index]) {
		index=index+1
		A = A + stim
	}
	SOLVE state METHOD sparse
	
	VERBATIM
	return 0;
	ENDVERBATIM
}

? kinetics
KINETIC state {
	~ A <-> G	(k, 0)
	~ G <-> bath	(k, 0)
}



























