TITLE pyrna.mod   DCN pyramidal cell model sodium channel
 
COMMENT

Revised version of DCN Pyramidal cell model sodium channel

This model implements part of a Dorsal Cochlear Nucleus Pyramidal point cell
based on kinetic data from Kanold and Manis (1999) and Kanold's dissertation (1999)

-- 15 Jan 1999 P. Manis

This mechanism is the fast sodium channel portion of the model.


Orignal: 2/10/02. P. Manis.

Extraced from Pyr.mod, 7/24/2014.

ENDCOMMENT
 
 
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}
 

NEURON {
    THREADSAFE
    SUFFIX napyr
    USEION na READ ena WRITE ina
    RANGE gna, minf, hinf, ninf, gbar	: sodium channels and delayed rectifier
    RANGE mtau, htau, ntau    : time constants for sodium channels and delayed rectifier
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
 v (mV)
 celsius (degC)
 dt (ms)
 ek (mV) : = -81.5 (mV)
 ena (mV) : = 50.0 (mV)
 gbar =  0.02857 (mho/cm2)	<0,1e9>
 mtau0 = 0.05 (ms) <0.01,100>
 htau0 = 0.5 (ms) <0.1,100>
 ntau = 0.5 (ms) <0.1,100>
 }
 
STATE {
        m h 
        }
 
ASSIGNED {
	gna (mho/cm2)
	ina (mA/cm2)
	minf hinf mtau htau
}

LOCAL mexp, hexp

BREAKPOINT {
	SOLVE states METHOD cnexp
	gna = gbar*m*m*h
	ina = gna*(v - ena)
}

UNITSOFF 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m) / mtau
	h' = (hinf - h) / htau
}
 
LOCAL q10


PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
	                      :Call once from HOC to initialize inf at resting v.
	LOCAL  alpha, beta, sum
	TABLE minf, mtau, hinf, htau  DEPEND celsius FROM -200 TO 100 WITH 400

UNITSOFF
	q10 = 3^((celsius - 22)/10)

	: "m" sodium activation system
	minf = na_m(v)
	mtau = na_mt(v)
	
	: "h" sodium inactivation system
	hinf = na_h(v)
	htau = na_ht(v)
	
}
 
: Make these as functions so we can view them from hoc, although this 
: may slow things down a bit

FUNCTION na_m(x) { : sodium activation
	na_m = 1/(1+exp(-(x+38)/3.0))	: POK version
:	na_m = alphbet(x,35,0,5,-10) :de Schutter (doesn't work well in our version)
:	na_m = na_m/(na_m + alphbet(x,7,0,65,20))
}

FUNCTION na_mt(x) { : sodium activation with taus
	na_mt = mtau0 : flat time constants
:	na_mt = alphbet(x,35,0,5,-10)
:	na_mt = 1/(na_mt + alphbet(x,7,0,65,20))
}

FUNCTION na_h(x) { : sodium inactivation
	na_h = 1/(1+exp((x+43)/3.0))	: flat time constants (POK version)
:	na_h = alphbet(x,0.225,1,80,10)	
:	na_h = na_h/(na_h + alphbet(x,7.5,0,-3,-18))
}

FUNCTION na_ht(x) { : sodium inactivation tau
	na_ht = htau0 : POK: flat time constants
:	na_ht = alphbet(x,0.225,1,80,10) : de Schutter version (doesn't work well with other stuff)
:	na_ht = 1/(na_ht + alphbet(x,7.5,0,-3,-18))
}


FUNCTION alphbet(x,A,B,C,D) { 	: alpha/beta general functions for
				: transcrbing GENESIS models
alphbet = A/(B+exp((x+C)/D))
}

UNITSON

