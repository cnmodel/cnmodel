:
:  ichanWT2005.mod 
:
:   Alan Goldin Lab, University of California, Irvine
:   Jay Lickfett - Last Modified: 6 July 2005
:
:  This file is the Nav1.1 wild-type channel model described in:
:
:		Barela et al. An Epilepsy Mutation in the Sodium Channel SCN1A That Decreases
:	    Channel Excitability.  J. Neurosci. 26(10): p. 2714-2723 
:
:
:   The model is derived from the one described in:
:
:    	Spampanato et al. (2004a) Increased Neuronal Firing in Computer Simulations 
:		of Sodium Channel Mutations that Cause Generalized Epilepsy with Febrile Seizures Plus.
:		Journal of Neurophysiology 91:2040-2050
:
:	and
:
:	 	Spampanato et al. (2004b) A Novel Epilepsy Mutation 
:   	in the Sodium Channel SCN1A Identifies a Cytoplasmic Domain for 
:		Beta Subunit Interaction. J. Neurosci. 24(44):10022-10034
:
 
: delayed rectifier removed (p.b.manis 2/22/2009)

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (uF) = (microfarad)
    (molar) = (1/liter)
    (nA) = (nanoamp)
    (mM) = (millimolar)
    (um) = (micron)
    (S) = (siemens)
    FARADAY = 96520 (coul)
    R = 8.3134  (joule/degC)

}

 
NEURON { 
    SUFFIX inav11 
    USEION nat READ enat WRITE inat VALENCE 1
    RANGE gnat
    RANGE gnatbar
    RANGE minf, mtau, hinf, htau, sinf, stau, inat, m, h, s
	RANGE vsna : voltage shift parameter
}

 
INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}

 
PARAMETER {
    
    celsius = 6.3 (degC)
    dt (ms) 
    enat  (mV)
    gnatbar (mho/cm2)   
}


ASSIGNED {
      
    v (mV) 
    gnat (mho/cm2) 
    inat (mA/cm2)
    minf hinf sinf
    mtau (ms) htau (ms) stau (ms)
    mexp hexp sexp
	vsna (mV)
} 


STATE {
    m h s
}
 

BREAKPOINT {
    SOLVE states
    gnat = gnatbar*m*m*m*h*s  
    inat = gnat*(v - enat)
}

 
UNITSOFF

 
INITIAL {

    trates(v)
    
    m = minf
    h = hinf
    s = sinf
    
    VERBATIM
    return 0;
    ENDVERBATIM
}


PROCEDURE states() {        : Computes state variables m, h, s and n 
                            : at the current v and dt.        
    trates(v)           

    m = m + mexp*(minf-m)
    h = h + hexp*(hinf-h)
    s = s + sexp*(sinf-s)    
    VERBATIM
    return 0;
    ENDVERBATIM
}
 

LOCAL q10


PROCEDURE rates(v (mV)) {   :Computes rate and other constants at current v.
                            :Call once from HOC to initialize inf at resting v.

    LOCAL  alpha, beta, sum
    q10 = 3^((celsius - 6.3)/10)
:    q10 = 3^((celsius - 22)/10)
    
    
    :"m" sodium activation system
    minf = 1/(1+exp(-(v+27.4+vsna)*4.7*0.03937))		
    mtau = 0.15

    :"h" sodium fast inactivation system
    hinf = 1/(1+exp((v+41.9+vsna)/6.7))				
    htau = 23.12*exp(-0.5*((v+77.58+vsna)/43.92)^2) 	
       
    :"s" sodium slow inactivation system
    sinf = 1/(1+exp((v+46.0+vsna)/6.6))						
    stau = 1000*(140.4*exp(-0.5*((v+71.3+vsna)/30.9)^2))   	


}
 

PROCEDURE trates(v (mV)) {  :Build table with rate and other constants at current v.
                            :Call once from HOC to initialize inf at resting v.
    LOCAL tinc
 
    TABLE minf, mexp, hinf, hexp, sinf, sexp, mtau, htau, stau
       DEPEND dt, celsius FROM -100 TO 100 WITH 200
                           
    rates(v)    : not consistently executed from here if usetable_hh == 1
                : so don't expect the tau values to be tracking along with
                : the inf values in hoc

    tinc = -dt * q10
    mexp = 1 - exp(tinc/mtau)
    hexp = 1 - exp(tinc/htau)
    sexp = 1 - exp(tinc/stau)
}

 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.

    if (fabs(x/y) < 1e-6) {
        vtrap = y*(1 - x/y/2)
    }else{  
        vtrap = x/(exp(x/y) - 1)
    }
}
 

UNITSON

