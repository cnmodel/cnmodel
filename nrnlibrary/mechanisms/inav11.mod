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
    THREADSAFE
    SUFFIX nav11
    USEION na READ ena WRITE ina VALENCE 1
    RANGE gna
    RANGE gbar
    RANGE minf, mtau, hinf, htau, sinf, stau, inat, m, h, s
	RANGE vsna : voltage shift parameter
}

 
INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}

 
PARAMETER {
    vsna = 4.3 (mV)
    celsius (degC)
    dt (ms) 
    ena (mV)
    :enat = 50  (mV)
    gbar = 0.1 (mho/cm2)
    q10 = 3.0 (1)
}


ASSIGNED {
      
    v (mV) 
    gna (mho/cm2)
    ina (mA/cm2)
    minf hinf sinf
    mtau (ms) htau (ms) stau (ms)
    mexp hexp sexp
:	vsna (mV)
} 


STATE {
    m h s
}
 

BREAKPOINT {
    SOLVE states METHOD cnexp
    gna = gbar*m*m*m*h*s
    ina = gna*(v - ena)
}

 
UNITSOFF

 
INITIAL {

    trates(v)
    
    m = minf
    h = hinf
    s = sinf

}


DERIVATIVE states {        : Computes state variables m, h, s and n
                            : at the current v and dt.        
    rates(v)
    m' = (minf - m)/mtau
    h' = (hinf - h)/htau
    s' = (sinf - s)/stau

}
 

LOCAL qt


PROCEDURE rates(v (mV)) {   :Computes rate and other constants at current v.
                            :Call once from HOC to initialize inf at resting v.

    LOCAL  alpha, beta, sum
    qt = q10^((celsius - 22)/10) : original recordings in Barela et al made at "room temperature"
    
    
    : "m" sodium activation system
    minf = f_minf(v)
    mtau = f_mtau(v)/qt

    : "h" sodium fast inactivation system
    hinf = f_hinf(v)
    htau = f_htau(v)/qt

    : "s" sodium slow inactivation system
    sinf = f_sinf(v)
    stau = f_stau(v)/qt

}
 

PROCEDURE trates(v (mV)) {  :Build table with rate and other constants at current v.
                            :Call once from HOC to initialize inf at resting v.
    LOCAL tinc
 
    TABLE minf, mexp, hinf, hexp, sinf, sexp, mtau, htau, stau
       DEPEND dt, celsius FROM -100 TO 100 WITH 200
                           
    rates(v)    : not consistently executed from here if usetable_hh == 1
                : so don't expect the tau values to be tracking along with
                : the inf values in hoc

    tinc = -dt  :  * q10 q10 is handled in rates, above
    mexp = 1 - exp(tinc/mtau)
    hexp = 1 - exp(tinc/htau)
    sexp = 1 - exp(tinc/stau)
}

FUNCTION f_minf(v (mV)) {
        f_minf = 1/(1+exp(-(v+27.4+vsna)*4.7*0.03937))

        }
FUNCTION f_mtau(v (mV)) {
        f_mtau = 0.15
        }

FUNCTION f_hinf(v (mV)) {
        f_hinf = 1/(1+exp((v+41.9+vsna)/6.7))
    }

FUNCTION f_htau(v (mV)) {
        f_htau = 23.12*exp(-0.5*((v+77.58+vsna)/43.92)^2)
        }


FUNCTION f_sinf(v (mV)) {
    f_sinf = 1/(1+exp((v+46.0+vsna)/6.6))
        }

FUNCTION f_stau(v (mV)) {
        f_stau = 1000*140.4*exp(-0.5*((v+71.3+vsna)/30.9)^2)
    }


UNITSON

