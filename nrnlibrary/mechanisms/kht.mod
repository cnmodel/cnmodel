TITLE kht.mod  The high threshold conductance of cochlear nucleus neurons

COMMENT

NEURON implementation of Jason Rothman's measurements of VCN conductances.

This file implements the high threshold potassium current found in several brainstem
 nuclei of the auditory system, including the spherical and globular bushy cells
  (Manis and Marx, 1991; Rothman and Manis, 2003a,b) and multipolar (stellate) 
  cells of the ventral cochlear nucleus, principal cells of the medial 
  nucleus of the trapzoid body (Brew and Forsythe, 1995, Wang and Kaczmarek, 
  1997) and neurons of the medial superior olive. The current is likely mediated by 
  Kv3.1  potassium channel subunits. The specific 
  implementation is described in Rothman and Manis, J. Neurophysiol. 2003, in the 
  appendix. Measurements were made from isolated neurons from adult guinea pig, 
  under reasonably stringent voltage clamp conditions. The measured current is 
  sensitive to 4-aminopyridine and TEA, but is spared by mamba snake toxi
  dendrotoxin I.


Similar conductrances are found in the homologous neurons of the avian auditory 
system (Reyes and Rubel; Zhang and Trussell; Rathouz and Trussell), and the 
conductance described here, in the absence of more detailed kinetic measurements
, is probably suitable for use in modeling that system.


Original implementation by Paul B. Manis, April (JHU) and Sept, (UNC)1999.

File split implementation, February 28, 2004.

Contact: pmanis@med.unc.edu

ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
        SUFFIX kht
        USEION k READ ek WRITE ik
        RANGE gkhtbar, gkht, ik
        GLOBAL ninf, pinf, ntau, ptau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius = 22 (degC)  : model is defined on measurements made at room temp in Baltimore
        dt (ms)
        ek = -77 (mV)
        gkhtbar = 0.01592 (mho/cm2) <0,1e9>
		nf = 0.85 <0,1> :proportion of n vs p kinetics
}

STATE {
        n p
}

ASSIGNED {
    ik (mA/cm) 
    gkht (mho/cm2)
    pinf ninf
    ptau (ms) ntau (ms)
    }

LOCAL nexp, pexp

BREAKPOINT {
	SOLVE states
    
	gkht = gkhtbar*(nf*(n^2) + (1-nf)*p)
    ik = gkht*(v - ek)

}

UNITSOFF

INITIAL {
    trates(v)
    p = pinf
    n = ninf
}

PROCEDURE states() {  :Computes state variables m, h, and n
	trates(v)      :             at the current v and dt.
	n = n + nexp*(ninf-n)
	p = p + pexp*(pinf-p)
VERBATIM
	return 0;
ENDVERBATIM
}

LOCAL q10

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

	q10 = 3^((celsius - 22)/10) : if you don't like room temp, it can be changed!

    ninf =   (1 + exp(-(v + 15) / 5))^-0.5
    pinf =  1 / (1 + exp(-(v + 23) / 6))

	ntau =  (100 / (11*exp((v+60) / 24) + 21*exp(-(v+60) / 23))) + 0.7
    ptau = (100 / (4*exp((v+60) / 32) + 5*exp(-(v+60) / 22))) + 5
}

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc
	TABLE ninf, nexp, pinf, pexp
	DEPEND dt, celsius FROM -150 TO 150 WITH 300

    rates(v)    : not consistently executed from here if usetable_hh == 1
        : so don't expect the tau values to be tracking along with
        : the inf values in hoc

	tinc = -dt * q10
	nexp = 1 - exp(tinc/ntau)
	pexp = 1 - exp(tinc/ptau)
	}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

UNITSON
