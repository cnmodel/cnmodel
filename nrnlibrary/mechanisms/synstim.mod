TITLE Synapse Stimulator

COMMENT
	Written for MyFirstNEURON by Arthur Houweling
ENDCOMMENT
					       
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS SynStim
	RANGE fAmpa,fNmda,fGabaA,fGabaB,sAmpa,sNmda,sGabaA,sGabaB,fGly,onset
}

UNITS {
	PI = (pi) (1)
}

PARAMETER {
	fAmpa
	fNmda
	fGabaA
	fGabaB
	fGly
	onset=	0 (ms)
}

ASSIGNED { 
	sAmpa   
	sNmda
	sGabaA
	sGabaB
	sGly
}

LOCAL c

INITIAL {
	c= PI*0.002
}

UNITSOFF
BREAKPOINT { LOCAL tt 
	if (t>onset) { tt= t-onset }
	else { tt= 0 }
	sAmpa= sin(c*fAmpa*tt)
	sNmda= sin(c*fNmda*tt)
	sGabaA= sin(c*fGabaA*tt)
	sGabaB= sin(c*fGabaB*tt)
	sGly = sin(c*fGly*tt)
}
UNITSON
