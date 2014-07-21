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
        THREADSAFE
        SUFFIX kht
        USEION k READ ek WRITE ik
        RANGE gbar, gkht, ik
        GLOBAL ninf, pinf, ntau, ptau
}

:INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
ASSIGNED {
    celsius (degC)  : model is defined on measurements made at room temp in Baltimore: 22 degC
    ik (mA/cm2)
    ek (mV)
    gkht (mho/cm2)
    pinf ninf
    ptau (ms)
    ntau (ms)
    qg ()  : computed q10 for gnabar based on q10g
    q10 ()
}
    
PARAMETER {
        v (mV)
        dt (ms)
        gbar = 0.01592 (mho/cm2) <0,1e9>
        nf = 0.85 <0,1> :proportion of n vs p kinetics
        q10tau = 3.0
        q10g = 2.0
}

STATE {
        n p
}

LOCAL nexp, pexp

BREAKPOINT {
    SOLVE states METHOD cnexp
    
    gkht = qg*gbar*(nf*(n^2) + (1-nf)*p)
    ik = gkht*(v - ek)
}

INITIAL {
    qg = q10g^((celsius-22)/10 (degC))
    q10 = q10tau^((celsius - 22)/10 (degC)) : if you don't like room temp, it can be changed!
    rates(v)
    p = pinf
    n = ninf
}

DERIVATIVE states {  :Computes state variables m, h, and n
    rates(v)      :             at the current v and dt.
    n' = (ninf - n)/ntau
    p' = (pinf - p)/ptau
}

PROCEDURE rates(v (mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

    ninf =   (1 + exp(-(v + 15) / 5 (mV)))^-0.5
    pinf =  1 / (1 + exp(-(v + 23) / 6 (mV)))

    ntau =  (100 (ms)/ (11*exp((v+60) / 24 (mV)) + 21*exp(-(v+60) / 23 (mV)))) + 0.7
    ntau = ntau/q10
    ptau = (100 (ms)/ (4*exp((v+60) / 32 (mV)) + 5*exp(-(v+60) / 22 (mV)))) + 5
    ptau = ptau/q10
    
}


