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
        THREADSAFE
        SUFFIX ka
        USEION k READ ek WRITE ik
        RANGE gbar, gka, ik
        GLOBAL ainf, binf, cinf, atau, btau, ctau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        dt (ms)
        gbar = 0.00477 (mho/cm2) <0,1e9>
        q10tau = 3.0
        q10g = 2.0
}

STATE {
        a b c
}

ASSIGNED {
    celsius  (degC)  : model is defined on measurements made at room temp in Baltimore
    ik (mA/cm2) 
    ek (mV)
    gka (mho/cm2)
    ainf binf cinf
    atau (ms) btau (ms) ctau (ms)
    qg ()  : computed q10 for gnabar based on q10g
    q10 ()
}

LOCAL aexp, bexp, cexp

BREAKPOINT {
    SOLVE states METHOD cnexp
    
    gka = gbar*(a^4)*b*c
    ik = gka*(v - ek)

}


INITIAL {
    qg = q10g^((celsius-22)/10 (degC))
    q10 = q10tau^((celsius - 22)/10 (degC)) : if you don't like room temp, it can be changed!
    rates(v)
    a = ainf
    b = binf
    c = cinf
}

DERIVATIVE states {  :Computes state variables m, h, and n
    rates(v)      :             at the current v and dt.
    a' = (ainf - a)/atau
    b' = (binf - b)/btau
    c' = (cinf - c)/ctau
}


PROCEDURE rates(v (mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

    ainf = (1 / (1 + exp(-1*(v + 31) / 6 (mV))))^0.25
    binf = 1 / (1 + exp((v + 66) / 7 (mV)))^0.5
    cinf = 1 / (1 + exp((v + 66) / 7 (mV)))^0.5

    atau =  (100 (ms)/ (7*exp((v+60) / 14 (mV)) + 29*exp(-(v+60) / 24 (mV)))) + 0.1
    atau = atau/q10
    btau =  (1000 (ms) / (14*exp((v+60) / 27 (mV)) + 29*exp(-(v+60) / 24 (mV)))) + 1
    btau = btau/q10
    ctau = (90 (ms)/ (1 + exp((-66-v) / 17 (mV)))) + 10
    ctau = ctau/q10
}

