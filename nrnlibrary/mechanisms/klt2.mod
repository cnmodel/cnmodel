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


Original implementation by Paul B. Manis, April (JHU) and Sept, (UNC) 1999.

File split implementation, February 28, 2004.
Updated to contemporary NMODL needs, 2012.

Contact: pmanis@med.unc.edu

ENDCOMMENT



NEURON {
    SUFFIX klt2
    THREADSAFE
    NONSPECIFIC_CURRENT i
    USEION k READ ek WRITE ik
    RANGE gkltbar, gklt, iklt
    RANGE winf, zinf, wtau, ztau
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (nA) = (nanoamp)
}

PARAMETER {
    gkltbar = 0.01592 (mho/cm2) <0,1e9>
    zss = 0.5   <0,1>   : steady state inactivation of glt
    Q10 = 3.0 (1)
    KineticMeasureTemp = 22 (degC)
}

ASSIGNED {
    v (mV)
    ik (mA/cm2)
    i (mA/cm2)
    ek (mV)
    gklt (mho/cm2)
    iklt (mA/cm2)
    winf (1) zinf (1)
    wtau (ms) ztau (ms)
    celsius (degC)  : model is defined on measurements made at room 
    qten (1)
}

STATE {
        w z
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    gklt = gkltbar*(w^4)*z
    iklt = gklt*(v - ek)
    i = iklt
    ik = i
}



INITIAL {
    qten = Q10^((celsius - KineticMeasureTemp)/10 (degC))  : if you don't like room
                                                    :temp, it can be changed!
    trates(v)
    w = winf
    z = zinf
}

DERIVATIVE states {  : Update state variables w and z
    rates(v)
    w' = (winf - w)/wtau
    z' = (zinf - z)/ztau
}

PROCEDURE rates(v (mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

    winf = (1 / (1 + exp(-(v + 48 (mV)) / 6 (mV))))^0.25
    wtau =  (100 (ms)/ (6*exp((v+60 (mV)) / 6 (mV)) + 16*exp(-(v+60 (mV)) / 45 (mV)))) + 1.5(ms)
    wtau = wtau/qten

    zinf = zss + ((1-zss) / (1 + exp((v + 71 (mV)) / 10 (mV))))
    ztau =  (1000(ms) / (exp((v+60 (mV)) / 20 (mV)) + exp(-(v+60 (mV)) / 8 (mV)))) + 50(ms)
    ztau = ztau/qten

}

PROCEDURE trates(v (mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
    TABLE winf, zinf
    DEPEND celsius FROM -250 TO 200 WITH 901

    rates(v)    : not consistently executed from here if usetable_hh == 1
        : so don't expect the tau values to be tracking along with
        : the inf values in hoc
}


