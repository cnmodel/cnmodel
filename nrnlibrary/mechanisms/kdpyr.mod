TITLE kdpyr.mod   DCN pyramidal cell model, delayed rectifier  
    
COMMENT

This is part of a model implements a Dorsal Cochlear Nucleus Pyramidal point cell
based on kinetic data from Kanold and Manis (1999) and Kanold's dissertation (1999)

-- 15 Jan 1999 P. Manis

Added export of start states for some variables to do perturbation tests
These start values replace the "inf" values used in the initialization procedure
Note that if the start variable is set to a value less than 0, 
then the default initialization will be done. Typically I use a value of -1 for this flagging
Note also that it is possible to set the initial values > 1 but this is meaningless in terms of
the present equations. 
-- 5 Feb 1999 P. Manis

ENDCOMMENT

UNITS {
       (mA) = (milliamp)
       (mV) = (millivolt)
}

NEURON {
    THREADSAFE
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
    kd_avh = -40 (mV)
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
    
BREAKPOINT {
	SOLVE states METHOD cnexp
	gk = gbar*n*n
	ik = gk*(v - ek)
}
    
INITIAL {
	rates(v)
	n = ninf
}

DERIVATIVE states {  
	rates(v)
	n' = (ninf - n) / ntau
}
    
LOCAL q10

PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
	                      :Call once from HOC to initialize inf at resting v.
	LOCAL  alpha, beta, sum
	TABLE ninf, ntau DEPEND celsius FROM -200 TO 100 WITH 400

	q10 = 3^((celsius - 22)/10 (degC))
	: "n" potassium activation system
       ninf = kd_m(v)
}
    
: Make these as functions so we can view them from hoc, although this 
: may slow things down a bit
FUNCTION kd_m(x (mV)) { : potassium activation
	kd_m = 1/(1+exp(-(x-kd_avh)/(3 (mV))))		: flat time constants
}

