TITLE klt.mod  The low threshold conductance of cochlear nucleus neurons

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

Modifed implementation; includes all temperature scaling, passes modlunit
7/10/2014 pbm

ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
        THREADSAFE
        SUFFIX na
        USEION na READ ena WRITE ina
        RANGE gbar, gna, ina
        GLOBAL hinf, minf, htau, mtau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        dt (ms)
        ena (mV)
        gbar =  0.07958 (mho/cm2) <0,1e9>
        q10tau = 3.0
        q10g = 2.0

}

STATE {
        m h
}

ASSIGNED {
    celsius (degC)  : model is defined on measurements made at room temp in Baltimore
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

INITIAL {
    qg = q10g^((celsius-22)/10 (degC))
    q10 = q10tau^((celsius - 22)/10 (degC)) : if you don't like room temp, it can be changed!
    rates(v)
    m = minf
    h = hinf
}

DERIVATIVE states {  :Computes state variables m, h, and n
    rates(v)      :             at the current v and dt.
    m' = (minf - m)/mtau
    h' = (hinf - h)/htau
}

PROCEDURE rates(v (mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

: average sodium channel
    minf = 1 / (1+exp(-(v + 38) / 7 (mV)))
    hinf = 1 / (1+exp((v + 65) / 6 (mV)))

    mtau =  (10 (ms)/ (5*exp((v+60) / 18 (mV)) + 36*exp(-(v+60) / 25 (mV)))) + 0.04
    mtau = mtau/q10
    htau =  (100 (ms)/ (7*exp((v+60) / 11 (mV)) + 10*exp(-(v+60) / 25 (mV)))) + 0.6
    htau = htau/q10
}

