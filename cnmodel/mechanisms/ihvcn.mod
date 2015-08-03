TITLE jsr.mod  VCN conductances

COMMENT
Ih for VCN neurons - average from several studies in auditory neurons

Implementation by Paul B. Manis, April (JHU) and Sept, (UNC)1999.
revised 2/28/04 pbm

pmanis@med.unc.edu

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
    SUFFIX ihvcn
    NONSPECIFIC_CURRENT i
    RANGE gbar, gh, ih, eh
    GLOBAL rinf, rtau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
    v (mV)
    dt (ms)
    gbar = 0.00318 (mho/cm2) <0,1e9>
    q10tau = 3.0
    q10g = 2.0
}

STATE {
    r
}

ASSIGNED {
    celsius (degC)
    gh (mho/cm2)
    eh (mV)
    i (mA/cm2)
    rinf
    rtau (ms)
    qg ()  : computed q10 for gnabar based on q10g
    q10 ()
}

BREAKPOINT {
    SOLVE states METHOD cnexp

    gh = qg*gbar*r
    i = gh*(v - eh)
}

INITIAL {
    qg = q10g^((celsius-22)/10 (degC))
    q10 = q10tau^((celsius - 22)/10 (degC)) : if you don't like room temp, it can be changed!
    rates(v)
    r = rinf
}

DERIVATIVE states {  :Computes state variables m, h, and n
    rates(v)      :         at the current v and dt.
    r' = (rinf - r)/rtau
}

PROCEDURE rates(v (mV)) {  :Computes rate and other constants at current v.
              :Call once from HOC to initialize inf at resting v.

    rinf = 1 / (1+exp((v + 76) / 7 (mV)))
    rtau = (100000 (ms)/ (237*exp((v+60) / 12 (mV)) + 17*exp(-(v+60) / 14 (mV)))) + 25
    rtau = rtau/q10

}

