TITLE pyr.mod   DCN pyramidal cell model  
 
COMMENT

Revised version of DCN Pyramidal cell model based on new hh.hoc file in NEURON

This model implements a Dorsal Cochlear Nucleus Pyramidal point cell
based on kinetic data from Kanold and Manis (1999) and Kanold's dissertation (1999)

-- 15 Jan 1999 P. Manis

Added export of start states for some variables to do perturbation tests
These start values replace the "inf" values used in the initialization procedure
Note that if the start variable is set to a value less than 0, 
then the default initialization will be done. Typically I use a value of -1 for this flagging
Note also that it is possible to set the initial values > 1 but this is meaningless in terms of
the present equations. 
-- 5 Feb 1999 P. Manis

Removed slow Ih current 30 Jan 2000. P. Manis
- also renamed variables to saner forms
Added Patrick's version of ih as ihd

Added persistent sodium current from deSchutter and Bower, J. Neurophys.
71:375, 1994.
Changed to Nap current from Quadroni and Knopfel, 1994


2/10/02. P. Manis.

ENDCOMMENT
 
 
UNITS {
 (mA) = (milliamp)
 (mV) = (millivolt)
}
 

NEURON {
 SUFFIX kdpyr
 USEION k READ ek WRITE ik
 RANGE gbar, gk		: delayed rectifier
 RANGE ntau: time constants delayed rectifier
 RANGE kd_avh
 
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
 v (mV)
 celsius (degC)
 dt (ms)
 ek (mV) : = -81.5 (mV)
 gbar = 0.006667 (mho/cm2)	<0,1e9>
 ntau = 0.5 (ms) <0.1,100>
 kd_avh = -40
 }
 
STATE {
         n
        }
 
ASSIGNED {
	gk (mho/cm2)
	ik (mA/cm2)
	ninf
}

LOCAL nexp
 
? currents
BREAKPOINT {
	SOLVE states METHOD cnexp
	gk = gbar*n*n
	ik = gk*(v - ek)
}
? currents

UNITSOFF 
 
INITIAL {
	rates(v)
	n = ninf
}


? states
DERIVATIVE states {  
	rates(v)
	n' = (ninf - n) / ntau
}
 
LOCAL q10


? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
	                      :Call once from HOC to initialize inf at resting v.
	LOCAL  alpha, beta, sum
	TABLE ninf, ntau DEPEND celsius FROM -200 TO 100 WITH 400

UNITSOFF
	q10 = 3^((celsius - 22)/10)

	
	: "n" potassium activation system
    ninf = kd_m(v)
}
 
: Make these as functions so we can view them from hoc, although this 
: may slow things down a bit
FUNCTION kd_m(x) { : potassium activation
	kd_m = 1/(1+exp(-(x-kd_avh)/3))		: flat time constants
}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

FUNCTION alphbet(x,A,B,C,D) { 	: alpha/beta general functions for
				: transcrbing GENESIS models
alphbet = A/(B+exp((x+C)/D))
}

UNITSON

