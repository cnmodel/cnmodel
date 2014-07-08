NEURON {
	SUFFIX capmp
	USEION ca READ cao, ica, cai WRITE cai, ica
	RANGE tau, width, cabulk, ica, pump0
}
	
UNITS {
	(um)	=	(micron)
	(molar) =	(1/liter)
	(mM)	=	(millimolar)
	(uM)	= 	(micromolar)
	(mA)	=	(milliamp)
	(mol)	=	(1)
	FARADAY = 	(faraday)	(coulomb)
}

PARAMETER {
	width = 0.1 (um)
	tau = 1 (ms)
	k1 = 5e8	(/mM-s)
	k2 = 0.25e6	(/s)
	k3 = 0.5e3	(/s)
	k4 = 5e0	(/mM-s)
	cabulk = 0.1 (uM)
	pump0 = 3e-14 (mol/cm2)
}

ASSIGNED {
	cao (mM) : 2
	cai (mM) : 100e-6
	ica (mA/cm2)
	ica_pmp (mA/cm2)
	ica_pmp_last (mA/cm2)
}

STATE {
	cam (uM)	<1e-6>
	pump (mol/cm2)  <1e-16>
	capump (mol/cm2) <1e-16>
}

INITIAL {
	ica = 0
	ica_pmp = 0
	ica_pmp_last = 0
	SOLVE pmp STEADYSTATE sparse
}

BREAKPOINT {
	SOLVE pmp METHOD sparse
	ica_pmp_last = ica_pmp
	ica = ica_pmp
}

KINETIC pmp {
	~ cabulk <-> cam (width/tau, width/tau)
	~ cam + pump <-> capump ((1e7)*k1, (1e10)*k2)
	~ capump <-> cao + pump ((1e10)*k3, (1e10)*k4)
	ica_pmp = (1e-7)*2*FARADAY*(f_flux - b_flux)

	: ica_pmp_last vs ica_pmp needed because of STEADYSTATE calculation
	~ cam << (-(ica - ica_pmp_last)/(2*FARADAY)*(1e7))

	CONSERVE pump + capump = (1e13)*pump0
	COMPARTMENT width {cam}	: volume has dimensions of um
	COMPARTMENT (1e13) {pump capump}	: area is dimensionless
	COMPARTMENT 1(um) {cabulk}
	COMPARTMENT (1e3)*1(um) {cao}

	cai = (0.001)*cam
}
