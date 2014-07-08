TITLE jsr.mod  VCN conductances

COMMENT
Ih for VCN neurons - average from several studies in auditory neurons


Implementation by Paul B. Manis, April (JHU) and Sept, (UNC)1999.
revised 2/28/04 pbm

pmanis@med.unc.edu

ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
        THREADSAFE
        SUFFIX ihvcn
        NONSPECIFIC_CURRENT i
        RANGE gbar, gh, ih, eh, vshift
        RANGE vh, k
        RANGE vtau, taumin
        RANGE tausc1, tausc2
        GLOBAL rinf, rtau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
:        celsius = 22 (degC)
        dt (ms)
        gbar = 0.00318 (mho/cm2) <0,1e9>
        eh = -43 (mV)
        vh = 76 (mV)
        k = 7 (mV)
        vtau = 60 (mV)
        taufac = 100000.0 (1)
        taumin = 25 (ms)
        tausc1 = 237 (ms)
        tausc2 = 17 (ms)
        vshift = 0 (mV)
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
	celsius (degC)
}

LOCAL rexp

BREAKPOINT {
	SOLVE states METHOD cnexp

	gh = gbar*r
    i = gh*(v - eh)
    }

UNITSOFF

INITIAL {
    trates(v)
    r = rinf
}

DERIVATIVE states {  :Computes state variables m, h, and n
	trates(v)      :             at the current v and dt.
	r' = (rinf - r)/rtau
:	r = r + rexp*(rinf-r)
:VERBATIM
:	return 0;
:ENDVERBATIM
}

LOCAL qt
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

	qt = q10^((celsius - 22.0)/10.0)
    rinf = 1.0 / (1+exp((v + vh + vshift) / k))
    rtau = (taufac / (tausc1*exp((v + vtau + vshift) / 12.0) + tausc2*exp(-(v + vtau + vshift) / 14.0)))
    rtau = rtau + taumin
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

	tinc = -dt : * q10 is handled in rates
	rexp = 1.0 - exp(tinc/rtau)
}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

UNITSON
