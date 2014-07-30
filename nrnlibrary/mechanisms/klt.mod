TITLE klt.mod  The low threshold conductance of cochlear nucleus neurons

COMMENT

NEURON implementation of Jason Rothman's measurements of VCN conductances.

This file implements the low threshold potassium current found in several brainstem
 nuclei of the auditory system, including the spherical and globular bushy cells
  (Manis and Marx, 1991; Rothman and Manis, 2003a,b) and octopus cells (Bal and
  Oertel, 2000) of the ventral cochlear nucleus, principal cells of the medial 
  nucleus of the trapzoid body (Brew and Forsythe, 1995, Wang and Kaczmarek, 
  1997) and neurons of the medial superior olive. The current is likely mediated by 
  heteromultimers of Kv1.1 and Kv1.2 potassium channel subunits. The specific 
  implementation is described in Rothman and Manis, J. Neurophysiol. 2003, in the 
  appendix. Measurements were made from isolated neurons from adult guinea pig, 
  under reasonably stringent voltage clamp conditions. The measured current is 
  sensitive to the mamba snake toxin dendrotoxin-I.


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
    SUFFIX klt
    USEION k READ ek WRITE ik
    RANGE gbar, gklt, ik
    GLOBAL winf, zinf, wtau, ztau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
    v (mV)
    dt (ms)
    gbar = 0.01592 (mho/cm2) <0,1e9>
    zss = 0.5   <0,1>   : steady state inactivation of glt
    q10tau = 3.0
    q10g = 2.0
}

STATE {
    w z
}

ASSIGNED {
    celsius (degC)  : model is defined on measurements made at room temp in Baltimore
    ik (mA/cm2) 
    ek (mV)
    gklt (mho/cm2)
    winf zinf
    wtau (ms) ztau (ms)
    qg ()  : computed q10 for gnabar based on q10g
    q10 ()
}

LOCAL wexp, zexp

BREAKPOINT {
    SOLVE states METHOD cnexp
    
    gklt = qg*gbar*(w^4)*z
    ik = gklt*(v - ek)
}

INITIAL {
    qg = q10g^((celsius-22)/10 (degC))
    q10 = q10tau^((celsius - 22)/10 (degC)) : if you don't like room temp, it can be changed!
    rates(v)
    w = winf
    z = zinf
}

DERIVATIVE states {  :Computes state variables m, h, and n
    rates(v)      :         at the current v and dt.
    w' = (winf - w)/wtau
    z' = (zinf - z)/ztau
}


PROCEDURE rates(v (mV)) {  :Computes rate and other constants at current v.
              :Call once from HOC to initialize inf at resting v.

    winf = (1 / (1 + exp(-(v + 48) / 6 (mV))))^0.25
    zinf = zss + ((1-zss) / (1 + exp((v + 71) / 10 (mV))))

    wtau =  (100 (ms)/ (6*exp((v+60) / 6 (mV)) + 16*exp(-(v+60) / 45 (mV)))) + 1.5
    wtau = wtau/q10
    ztau =  (1000 (ms)/ (exp((v+60) / 20 (mV)) + exp(-(v+60) / 8 (mV)))) + 50
    ztau = ztau/q10
}

