TITLE KIR channel
COMMENT
Reference: Steephen JE, Manchanda R (2009) Differences in biophysical
properties of nucleus accumbens medium spiny neurons emerging from
inactivation of inward rectifying potassium currents. J Comput Neurosci [PubMed] 

Found on ModelDB, 1/21/2013 PBManis


ENDCOMMENT

NEURON {
	SUFFIX KIR
	USEION  k READ ek WRITE ik
	RANGE g, ik, gmax
	GLOBAL minf, mtau
}

UNITS {
	(mA) = (milliamp)
	(uA) = (microamp)
	(mV) = (millivolt)
	(mS) = (millimho)
}

PARAMETER {
	celsius (degC)
	ek			(mV)
	gmax = 1.4e-4	(mho/cm2)	<0,1e9>
	m_vh = -82	(mV)	: half activation
	m_ve = 13		(mV)	: slope
}

ASSIGNED {
	v	(mV)
	g	(mho/cm2)
	ik	(mA/cm2)
	minf	(1)
	mtau	(ms)
	qt (1)
}

STATE {
	m
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gmax*m
	ik = g*(v - ek)
}

INITIAL {
	qt = 3^((celsius-35)/10)
	rates(v)
	m = minf
}

DERIVATIVE states { 
	rates(v)
	m' = (minf-m)/mtau
}

FUNCTION_TABLE tabmtau(v(mV)) (ms)

: rates() computes rate and other constants at present v
: call once from hoc to initialize inf at resting v

PROCEDURE rates(v(mV)) {
:	mtau = tabmtau(v)
	mtau = tabmtau(v)/qt
	minf = 1/(1 + exp((v - m_vh)/m_ve))
}
 
COMMENT
/* TABLES
The tables here are built as vecs in hoc, and then loaded into the mod file arrays before use...

 * 
 * Steephen, J. E., & Manchanda, R. (2009). Differences in biophysical properties of nucleus accumbens medium spiny neurons emerging from inactivation of inward rectifying potassium currents. J Comput Neurosci, 
 * doi:10.1007/s10827-009-0161-7
 */
 
//KIR
objref vecmtau_KIR, vecv_KIR
vecmtau_KIR = new Vector()
vecv_KIR = new Vector()
vecv_KIR.indgen(-120, 0, 10)
vecmtau_KIR.append(7.465, 7.465, 7.465, 8, 9.435, 10.755, 12.12, 13.795, 15.385, 14.285, 11.765, 8.89, 8) // At 35 deg C
table_tabmtau_KIR(&vecmtau_KIR.x[0], vecv_KIR.size, &vecv_KIR.x[0]) 

//inKIR
objref vecv_inKIR, vechinf_inKIR, vechtau_inKIR, vecv_tau_inKIR
vecv_tau_inKIR = new Vector()
vechtau_inKIR= new Vector()
vecv_tau_inKIR.append(-120,-90, -50)
vechtau_inKIR.append(7.767, 15, 25.333) // At 35 deg C
table_tabhtau_inKIR(&vechtau_inKIR.x[0],vecv_tau_inKIR.size, &vecv_tau_inKIR.x[0])

vecv_inKIR = new Vector()
vechinf_inKIR = new Vector()
vecv_inKIR.append(-120,-90, -50)
vechinf_inKIR.append(0, 0.13, 1)
table_tabhinf_inKIR(&vechinf_inKIR.x[0], vecv_inKIR.size, &vecv_inKIR.x[0])
table_tabmtau_inKIR(&vecmtau_KIR.x[0], vecv_KIR.size, &vecv_KIR.x[0]) 

ENDCOMMENT

