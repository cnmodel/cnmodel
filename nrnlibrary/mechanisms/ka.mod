TITLE klt.mod  The low threshold conductance of cochlear nucleus neurons

COMMENT

NEURON implementation of Jason Rothman's measurements of VCN conductances.

This file implements the transient potassium current found in ventral cochlear
nucleus "Type I" cells, which are largely "stellate" or "multipolar" cells  (Manis and
Marx, 1991; Rothman and Manis, 2003a,b; Manis et al, 1996). The current is likely
 mediated by Kv4.2 potassium channel subunits, but this has not been directly
demonstrated. The specific implementation is described in Rothman and Manis, J.
Neurophysiol. 2003, in the appendix. Measurements were made from isolated 
neurons from adult guinea pig, under reasonably stringent voltage clamp conditions.
 The measured current is sensitive to 4-aminopyridine. 
Original implementation by Paul B. Manis, April (JHU) and Sept, (UNC)1999.

File split implementaiton, April 1, 2004.

Contact: pmanis@med.unc.edu

ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
        SUFFIX ka
        USEION k READ ek WRITE ik
        RANGE gkabar, gka, ik
        GLOBAL ainf, binf, cinf, atau, btau, ctau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius = 22 (degC)  : model is defined on measurements made at room temp in Baltimore
        dt (ms)
        ek = -77 (mV)
        gkabar = 0.00477 (mho/cm2) <0,1e9>
}

STATE {
        a b c
}

ASSIGNED {
    ik (mA/cm2) 
    gka (mho/cm2)
    ainf binf cinf
    atau (ms) btau (ms) ctau (ms)
    }

LOCAL aexp, bexp, cexp

BREAKPOINT {
	SOLVE states
    
	gka = gkabar*(a^4)*b*c
    ik = gka*(v - ek)

}

UNITSOFF

INITIAL {
    trates(v)
    a = ainf
    b = binf
    c = cinf
}

PROCEDURE states() {  :Computes state variables m, h, and n
	trates(v)      :             at the current v and dt.
	a = a + aexp*(ainf-a)
	b = b + bexp*(binf-b)
	c = c + cexp*(cinf-c)
VERBATIM
	return 0;
ENDVERBATIM
}

LOCAL q10

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

	q10 = 3^((celsius - 22)/10) : if you don't like room temp, it can be changed!

    ainf = (1 / (1 + exp(-1*(v + 31) / 6)))^0.25
    binf = 1 / (1 + exp((v + 66) / 7))^0.5
    cinf = 1 / (1 + exp((v + 66) / 7))^0.5

    atau =  (100 / (7*exp((v+60) / 14) + 29*exp(-(v+60) / 24))) + 0.1
    btau =  (1000 / (14*exp((v+60) / 27) + 29*exp(-(v+60) / 24))) + 1
    ctau = (90 / (1 + exp((-66-v) / 17))) + 10
}

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc
	TABLE ainf, aexp, binf, bexp, cinf, cexp
	DEPEND dt, celsius FROM -150 TO 150 WITH 300

    rates(v)    : not consistently executed from here if usetable_hh == 1
        : so don't expect the tau values to be tracking along with
        : the inf values in hoc

	tinc = -dt * q10
	aexp = 1 - exp(tinc/atau)
	bexp = 1 - exp(tinc/btau)
	cexp = 1 - exp(tinc/ctau)
	}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

UNITSON
