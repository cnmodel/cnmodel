TITLE jsr.mod  VCN conductances

COMMENT
Ih for VCN neurons - average from several studies in auditory neurons


Implementation by Paul B. Manis, April (JHU) and Sept, (UNC)1999.
revised 2/28/04 pbm
Naming convention modified 2/10/2014
pmanis@med.unc.edu

ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
THREADSAFE
        SUFFIX ih
        NONSPECIFIC_CURRENT i
        RANGE gbar, gh, ih, eh
        GLOBAL rinf, rtau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius (degC) : data collected at 22
        dt (ms)
        gbar = 0.00318 (mho/cm2) <0,1e9>
        eh = -43 (mV)
        q10 = 3.0 (1)
}

STATE {
        r
}

ASSIGNED {
	gh (mho/cm2)
	i (mA/cm2)
	rinf
    rtau (ms)
}

LOCAL rexp

BREAKPOINT {
	SOLVE states METHOD derivimplicit
    
	gh = gbar*r
    i = gh*(v - eh)
    }

UNITSOFF

INITIAL {
    trates(v)
    r = rinf
}

DERIVATIVE states {  :Computes state variables m, h, and n
	rates(v)      :             at the current v and dt.
	r' = (rinf - r)/rtau
        :r = r + rexp*(rinf-r)
:VERBATIM
:	return 0;
:ENDVERBATIM
}

LOCAL qt
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

	qt = q10^((celsius - 22.0)/10.0)
    rinf = 1.0 / (1+exp((v + 76.0) / 7.0))
    rtau = (100000.0 / (237.0*exp((v+60.0) / 12.0) + 17.0*exp(-(v+60.0) / 14.0))) + 25.0
    rtau = rtau/qt

}

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc
	TABLE rinf, rexp
	DEPEND dt, celsius FROM -200 TO 150 WITH 350

    rates(v)    : not consistently executed from here if usetable_hh == 1
        : so don't expect the tau values to be tracking along with
        : the inf values in hoc

	tinc = -dt  :  * q10
	rexp = 1 - exp(tinc/rtau)
}


UNITSON
