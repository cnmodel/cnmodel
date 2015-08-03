TITLE nacn.mod  A sodium conductance for a ventral cochlear nucleus neuron model

COMMENT

NEURON implementation of Jason Rothman's measurements of VCN conductances.

This file implements the average brain sodium current used in the Rothman model.
In the absence of direct measurements in the VCN, this is a fair assumption.
The model differs from the one used in Rothman et al, (1993) in that the steep
voltage dependence of recovery from inactivation in that model is missing. This
may affect the refractory period. To use the other model, use najsr.mod instead.

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
        SUFFIX nacn
        USEION na READ ena WRITE ina
        RANGE gbar, gna, ina
        GLOBAL hinf, minf, htau, mtau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius (degC) : 22 (degC) model is defined on measurements made at room temp in Baltimore
        dt (ms)
        ena (mV)
        gbar =  0.07958 (mho/cm2) <0,1e9>
        q10 = 3.0 : q10 for rates
}

STATE {
        m h
}

ASSIGNED {
    ina (mA/cm2) 
    gna (mho/cm2)
    minf hinf
    mtau (ms) htau (ms)
    }

LOCAL mexp, hexp

BREAKPOINT {
	SOLVE states
    
    gna = gbar*(m^3)*h
    ina = gna*(v - ena)

}

UNITSOFF

INITIAL {
    trates(v)
    m = minf
    h = hinf
}

PROCEDURE states() {  :Computes state variables m, h, and n
	trates(v)      :             at the current v and dt.
	m = m + mexp*(minf-m)
	h = h + hexp*(hinf-h)
VERBATIM
	return 0;
ENDVERBATIM
}

LOCAL qt

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

	qt = q10^((celsius - 22)/10) : if you don't like room temp, it can be changed!

: average sodium channel
    minf = 1 / (1+exp(-(v + 38) / 7))
    hinf = 1 / (1+exp((v + 65) / 6))

    mtau =  (10 / (5*exp((v+60) / 18) + 36*exp(-(v+60) / 25))) + 0.04
    mtau = mtau/qt
    htau =  (100 / (7*exp((v+60) / 11) + 10*exp(-(v+60) / 25))) + 0.6
    htau = htau/qt
}

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc
	TABLE minf, mexp, hinf, hexp
	DEPEND dt, celsius FROM -150 TO 150 WITH 300

    rates(v)    : not consistently executed from here if usetable_hh == 1
        : so don't expect the tau values to be tracking along with
        : the inf values in hoc

	tinc = -dt :  * qt (note q10 is handled in mtau/htau calculation above
	mexp = 1 - exp(tinc/mtau)
	hexp = 1 - exp(tinc/htau)
	}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

UNITSON
