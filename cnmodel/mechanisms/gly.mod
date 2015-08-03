TITLE Gly synapse

COMMENT
	MODIFIED to be a faster GLY synapse, taken from GABA synapse
	Paul B. Manis - 7 Feb 2000

        simple alpha-synapse that generates a single PSP   
        *********************************************
        reference:	McCormick, Wang & Huguenard (1993) 
			Cerebral Cortex 3(5), 387-398
        found in:       cat reticular nucleus of thalamus
        *********************************************
	Assembled for MyFirstNEURON by Arthur Houweling


ENDCOMMENT
					       
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS GlySynapse
	USEION cl READ ecl VALENCE 1
	: negative valence not accepted by nrnivmodl
	RANGE onset, gmaxIPSP, e, g, i, w
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) 	= (nanoamp)
	(mV)	= (millivolt)
	(nS) 	= (nanomho)
}

PARAMETER {
	onset= 25	(ms)
	gmaxIPSP= 0	(nS)
	w= 1 				: weight factor for gmaxIPSP
	ecl		(mV)
	v		(mV)
	celsius		(degC)
}

ASSIGNED { 
	i 		(nA)  
	g 		(nS)
	tadj
}

UNITSOFF
INITIAL {
        tadj = 3^((celsius-23.5)/10)
}    

BREAKPOINT { LOCAL tt
	tt= (t-onset)*tadj
	if ((t>onset)&&(tt<740)) {
	: the exp() function does not accept arguments smaller than -745
	  g = w*gmaxIPSP * exp(-tt/15) * (1-exp(-tt/0.5))/0.84
	}
	else {g = 0}
	: -ecl because negative valences can not be specified
	i = g * (v-(-ecl))	
}
UNITSON

