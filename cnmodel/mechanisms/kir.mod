TITLE kir.mod 
 
COMMENT

Inwardly rectifying potassium (Kir) currents
based on mathematical model from Leao et al (2012)
 Ceballos CC, Li S, Roque AC, Tzounopoulos T, Leão RM (2016) Ih Equalizes Membrane Input Resistance in a Heterogeneous Population of Fusiform Neurons in the Dorsal Cochlear Nucleus. Front Cell Neurosci 10:249 [PubMed]
ModelDB Accession: 206252
Reduced from fusiform.mod from Leao et al. 2012, 

ENDCOMMENT
 
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}

 
NEURON {
    SUFFIX kir
    USEION k READ ek WRITE ik

	RANGE gk, gbar 	: Persistent Na channels and KIR channels 
	GLOBAL ntau, ninf	: time constants for Na channels and K channels 

}

 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
    v (mV)
    celsius (degC)
    dt (ms)
	ek (mV)
    ena (mV)
	ekir (mV)
    gbar = 0.0005 (mho/cm2)    <0,1e9>
	ntau = 0.5 (ms) <0,100>
}


STATE {
    nir
}
 
ASSIGNED {
    gk (mho/cm2)

    ik (mA/cm2)

    ninf 
}


LOCAL nexp
 
? currents
BREAKPOINT {
    SOLVE states METHOD cnexp

    gk = gbar*nir
	ik = gk*(v - ek)
}
? currents

UNITSOFF 
 

INITIAL {
	rates(v)
	nir = ninf
}


? states
DERIVATIVE states {  
    rates(v)
    nir' = (ninf - nir) / ntau
}
 
LOCAL q10

? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
	                      :Call once from HOC to initialize inf at resting v.
	TABLE ninf DEPEND celsius FROM -200 TO 100 WITH 400

UNITSOFF
	q10 = 3^((celsius - 32)/10)

	:"n" potassium activation system
        ninf = kird_m(v)
}

FUNCTION kird_m(x) { : KIR activation
	kird_m = 1/(1+exp((x+85.48)/12)) 
}

UNITSON