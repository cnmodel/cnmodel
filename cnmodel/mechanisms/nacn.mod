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

 Does not pass modlunit.
Should work at 22C and scales by Rothman and Manis, 2003c for temperature

Contact: pmanis@med.unc.edu

ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
        THREADSAFE
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
        q10tau = 3.0 : q10 for rates
        q10g = 2.0 : q10 for conductances
}

STATE {
        m h
}

ASSIGNED {
    ina (mA/cm2) 
    gna (mho/cm2)
    minf hinf
    mtau (ms) htau (ms)
    qg ()  : computed q10 for gnabar based on q10g
    q10 ()
    
    }

LOCAL mexp, hexp

BREAKPOINT {
	SOLVE states METHOD cnexp
    
    gna = qg*gbar*(m^3)*h
    ina = gna*(v - ena)
}

UNITSOFF

INITIAL {
    rates(v)
    m = minf
    h = hinf
}

DERIVATIVE states {  :Computes state variables m, h, and n
	rates(v)      :             at the current v and dt.
	m' = (minf - m)/mtau : m = m + mexp*(minf-m)
	h' = (hinf - h)/htau : h = h + hexp*(hinf-h)

}

LOCAL qt

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

    qg = q10g^((celsius-22)/10 (degC))
	q10 = q10tau^((celsius - 22)/10) : if you don't like room temp, it can be changed!

: average sodium channel
    minf = 1 / (1+exp(-(v + 38) / 7))
    hinf = 1 / (1+exp((v + 65) / 6))

    mtau =  (10 / (5*exp((v+60) / 18) + 36*exp(-(v+60) / 25))) + 0.04
    mtau = mtau/q10
    htau =  (100 / (7*exp((v+60) / 11) + 10*exp(-(v+60) / 25))) + 0.6
    htau = htau/q10
}


UNITSON
