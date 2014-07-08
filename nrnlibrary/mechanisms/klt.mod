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
:        celsius = 22 (degC)  : model is defined on measurements made at room temp in Baltimore
        dt (ms)
:        ek = -77 (mV)
        gbar = 0.01592 (mho/cm2) <0,1e9>
        zss = 0.5   <0,1>   : steady state inactivation of glt
        celsius (degC)
        q10 = 3.0 (1)
        ek (mV)
}

STATE {
        w z
}

ASSIGNED {
    ik (mA/cm2) 
    gklt (mho/cm2)
    winf zinf
    wtau (ms) ztau (ms)

}

LOCAL wexp, zexp

BREAKPOINT {
	SOLVE states METHOD cnexp
    
	gklt = gbar*(w^4)*z
    ik = gklt*(v - ek)

}

:UNITSOFF

INITIAL {
    trates(v)
    w = winf
    z = zinf
}

DERIVATIVE states {  :Computes state variables m, h, and n
	rates(v)      :             at the current v and dt.
	w' = (winf - w)/wtau
    z' = (zinf - z)/ztau
:        w = w + wexp*(winf-w)
:	z = z + zexp*(zinf-z)
:VERBATIM
:	return 0;
:ENDVERBATIM
}

LOCAL qt

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

	qt = q10^((celsius - 22)/10) : if you don't like room temp, it can be changed!

    winf = (1 / (1 + exp(-(v + 48) / 6)))^0.25
    zinf = zss + ((1-zss) / (1 + exp((v + 71) / 10)))

    wtau =  (100 / (6*exp((v+60) / 6) + 16*exp(-(v+60) / 45))) + 1.5
    wtau = wtau/q10
    ztau =  (1000 / (exp((v+60) / 20) + exp(-(v+60) / 8))) + 50
    ztau = ztau/q10
}

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc
	TABLE winf, wexp, zinf, zexp
	DEPEND dt, celsius FROM -150 TO 150 WITH 300

    rates(v)    : not consistently executed from here if usetable_hh == 1
        : so don't expect the tau values to be tracking along with
        : the inf values in hoc

:	tinc = -dt  :  * qt (note qt is handled in rates)
:	wexp = 1 - exp(tinc/wtau)
:	zexp = 1 - exp(tinc/ztau)
	}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

:UNITSON
