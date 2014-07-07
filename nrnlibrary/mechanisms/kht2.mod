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
    SUFFIX kht2
    THREADSAFE
    NONSPECIFIC_CURRENT i
    USEION k READ ek WRITE ik
    RANGE gkhtbar, gkht
    RANGE ninf, pinf, ntau, ptau
}


PARAMETER {
    dt (ms)
    gkhtbar = 0.01592 (mho/cm2) <0,1e9>
    nf = 0.85 <0,1> :proportion of n vs p kinetics
    Q10 = 3.0 (1)
    KineticMeasureTemp = 22 (degC) : model is defined on measurements made
                                    : at room temp in Baltimore
}

STATE {
    n p
}

ASSIGNED {
    v (mV)
    ik (mA/cm2) 
    i (mA/cm2)
    gkht (mho/cm2)
    pinf (1) ninf (1)
    ptau (ms) ntau (ms)
    ek (mV)
    celsius (degC)
    qten (1)
}


BREAKPOINT {
    SOLVE states METHOD cnexp
    gkht = gkhtbar*(nf*(n^2) + (1-nf)*p)
    i = gkht*(v - ek)
    ik = i
}


INITIAL {
    qten = Q10^((celsius - KineticMeasureTemp)/10 (degC)) : if you don't like room temp, it can be changed!
    rates(v)
    p = pinf
    n = ninf
}

DERIVATIVE states {  :Computes state variables m, h, and n
    trates(v)      :             at the current v and dt.
    n' = (ninf - n)/ntau 
    p' = (pinf - p)/ptau

}

PROCEDURE rates(v (mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

    ninf =   (1 + exp(-(v + 15 (mV)) / 5 (mV)))^-0.5
    ntau =  (100 (ms) / (11*exp((v + 60(mV)) / 24 (mV)) + 21*exp(-(v+60 (mV)) / 23 (mV)))) + 0.7 (ms)
    ntau = ntau/qten
    
    pinf =  1 / (1 + exp(-(v + 23 (mV)) / 6 (mV)))
    ptau = (100 (ms)/ (4*exp((v+60 (mV)) / 32 (mV)) + 5*exp(-(v+60 (mV)) / 22 (mV)))) + 5 (ms)
    ptau = ptau/qten

}

PROCEDURE trates(v (mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
    TABLE ninf, pinf
    DEPEND celsius FROM -250 TO 200 WITH 901

    rates(v)    : not consistently executed from here if usetable_hh == 1
        : so don't expect the tau values to be tracking along with
        : the inf values in hoc

}
