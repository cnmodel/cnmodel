TITLE kinetic NMDA receptor model

COMMENT
-----------------------------------------------------------------------------

	Kinetic model of NMDA receptors
	===============================

	10-state gating model:
	Kampa et al. (2004) J Physiol
  
	  U -- Cl  --  O
         \   | \	    \
          \  |  \      \
         UMg --  ClMg - OMg
		 |	|
		D1	|
		 | \	|
		D2  \	|
		   \	D1Mg
		    \	|
			D2Mg
-----------------------------------------------------------------------------

  Based on voltage-clamp recordings of NMDA receptor-mediated currents in 
  nucleated patches of  rat neocortical layer 5 pyramidal neurons (Kampa 2004), 
  this model was fit with AxoGraph directly to experimental recordings in 
  order to obtain the optimal values for the parameters.

-----------------------------------------------------------------------------

  This mod file does not include mechanisms for the release and time course
  of transmitter; it should to be used in conjunction with a sepearate mechanism
  to describe the release of transmitter and tiemcourse of the concentration
  of transmitter in the synaptic cleft (to be connected to pointer XMTR here).

-----------------------------------------------------------------------------

  See details of NEURON kinetic models in:

  Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  Kinetic models of 
  synaptic transmission.  In: Methods in Neuronal Modeling (2nd edition; 
  edited by Koch, C. and Segev, I.), MIT press, Cambridge, 1996.


  Written by Bjoern Kampa in 2004 
  Lightly modified, Paul Manis 2010.
	Note that data were taken at 23 deg C
	Q10 was taken from native receptors:
	Korinek M, Sedlacek M, Cais O, Dittert I, Vyklicky L Jr. Temperature
dependence of N-methyl-D-aspartate receptor channels and N-methyl-D-aspartate
receptor excitatory postsynaptic currents. Neuroscience. 2010 Feb
3;165(3):736-48. Epub 2009 Oct 31. PubMed PMID: 19883737.
	
-----------------------------------------------------------------------------
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
THREADSAFE
THREADSAFE
	POINT_PROCESS NMDA_Kampa
	POINTER XMTR
	RANGE U, Cl, D1, D2, Open, MaxOpen, UMg, ClMg, D1Mg, D2Mg, OMg
	RANGE g, gmax, vshift, Erev, rb, rmb, rmu, rbMg,rmc1b,rmc1u,rmc2b,rmc2u
	GLOBAL Erev, mg, Rb, Ru, Rd1, Rr1, Rd2, Rr2, Ro, Rc, Rmb, Rmu
	GLOBAL RbMg, RuMg, Rd1Mg, Rr1Mg, Rd2Mg, Rr2Mg, RoMg, RcMg
	GLOBAL Rmd1b,Rmd1u,Rmd2b,Rmd2u,rmd1b,rmd1u,rmd2b,rmd2u
	GLOBAL Rmc1b,Rmc1u,Rmc2b,Rmc2u
	GLOBAL vmin, vmax, valence, memb_fraction
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(umho) = (micromho)
	(mM) = (milli/liter)
	(uM) = (micro/liter)
}

PARAMETER {

	Erev	= 5    	(mV)	: reversal potential
	gmax	= 500  	(pS)	: maximal conductance
	mg	= 1  	(mM)	: external magnesium concentration
	vmin 	= -120	(mV)
	vmax 	= 100	(mV)
	valence = -2		: parameters of voltage-dependent Mg block
	memb_fraction = 0.8
	vshift = 0.0 (mV)
	Q10 = 2.0 : temperature sensitivity (see above)
	
    : Maximum open probability with Mode=0 (no rectification). 
    : This is determined empirically by holding XMTR at a large
    : value and v=40mV for 100 timesteps and measuring the 
    : maximum value of Open.
    MaxOpen = 0.01988893957 (1) 

: Rates

	Rb		= 10e-3   	(/uM /ms)	: binding 		
	Ru		= 5.6e-3  	(/ms)	: unbinding		
	Ro		= 10e-3   	(/ms)	: opening
	Rc		= 273e-3   	(/ms)	: closing 
:	Rd1		= 2.2e-3   	(/ms)	: fast desensitisation
	Rd1		= 0.1		(/ms) 	: fast desensitisation
	Rr1		= 1.6e-3   	(/ms)	: fast resensitisation
:	Rd2 	= 0.43e-3 	(/ms)	: slow desensitisation
	Rd2 	= 1e-4 		(/ms) 	: slow desensitisation
	Rr2 	= 0.5e-3	(/ms)	: slow resensitisation
	Rmb		= 0.05e-3	(/uM /ms)	: Mg binding Open
	Rmu		= 12800e-3	(/ms)	: Mg unbinding Open
	Rmc1b		= 0.00005e-3	(/uM /ms)	: Mg binding Closed
	Rmc1u		= 2.438312e-3	(/ms)	: Mg unbinding Closed
	Rmc2b		= 0.00005e-3	(/uM /ms)	: Mg binding Closed2
	Rmc2u		= 5.041915e-3	(/ms)	: Mg unbinding Closed2
	Rmd1b		= 0.00005e-3	(/uM /ms)	: Mg binding Desens1
	Rmd1u		= 2.98874e-3	(/ms)	: Mg unbinding Desens1
	Rmd2b		= 0.00005e-3	(/uM /ms)	: Mg binding Desens2
	Rmd2u		= 2.953408e-3	(/ms)	: Mg unbinding Desens2
	RbMg		= 10e-3		(/uM /ms)	: binding with Mg
	RuMg		= 17.1e-3	(/ms)	: unbinding with Mg
	RoMg		= 10e-3		(/ms)	: opening with Mg
	RcMg		= 548e-3	(/ms)	: closing with Mg
	Rd1Mg		= 2.1e-3	(/ms)	: fast desensitisation with Mg
	Rr1Mg		= 0.87e-3	(/ms)	: fast resensitisation with Mg
	Rd2Mg		= 0.26e-3	(/ms)	: slow desensitisation with Mg
	Rr2Mg		= 0.42e-3	(/ms)	: slow resensitisation with Mg
}

ASSIGNED {
	v		(mV)	: postsynaptic voltage
	i 		(nA)	: current = g*(v - Erev)
	g 		(pS)	: conductance
	XMTR 		(mM)	: pointer to glutamate concentration

	rb		(/ms)   : binding, [glu] dependent
	rmb		(/ms)	: blocking V and [Mg] dependent
	rmu		(/ms)	: unblocking V and [Mg] dependent
	rbMg		(/ms)	: binding, [glu] dependent
	rmc1b		(/ms)	: blocking V and [Mg] dependent
	rmc1u		(/ms)	: unblocking V and [Mg] dependent
	rmc2b		(/ms)	: blocking V and [Mg] dependent
	rmc2u		(/ms)	: unblocking V and [Mg] dependent
	rmd1b		(/ms)	: blocking V and [Mg] dependent
	rmd1u		(/ms)	: unblocking V and [Mg] dependent
	rmd2b		(/ms)	: blocking V and [Mg] dependent
	rmd2u		(/ms)	: unblocking V and [Mg] dependent

	qfac  : Q10
    celsius (degC)
}

STATE {
	: Channel states (all fractions)
	U		: unbound
	Cl		: closed
	D1		: desensitised 1
	D2		: desensitised 2
	Open		: open
	UMg		: unbound with Mg
	ClMg		: closed with Mg
	D1Mg		: desensitised 1 with Mg
	D2Mg		: desensitised 2 with Mg
	OMg		: open with Mg
}

INITIAL {
	U = 1
	qfac = Q10^((celsius-23)/10 (degC))}

BREAKPOINT {
	SOLVE kstates METHOD sparse

	g = gmax * Open / MaxOpen
	i = (1e-6) * g * (v - Erev)
}

KINETIC kstates {

	rb 	= Rb 	* (1e3) * XMTR
	rbMg 	= RbMg 	* (1e3) * XMTR
	rmb 	= Rmb 	* mg * (1e3) * exp((v-40+vshift) * valence * memb_fraction /25 (mV))
	rmu 	= Rmu 	* exp((-1)*(v-40+vshift) * valence * (1-memb_fraction) /25 (mV))
	rmc1b 	= Rmc1b * mg * (1e3) * exp((v-40+vshift) * valence * memb_fraction /25 (mV))
	rmc1u 	= Rmc1u * exp((-1)*(v-40+vshift) * valence * (1-memb_fraction) /25 (mV))
	rmc2b 	= Rmc2b * mg * (1e3) * exp((v-40+vshift) * valence * memb_fraction /25 (mV))
	rmc2u 	= Rmc2u * exp((-1)*(v-40+vshift) * valence * (1-memb_fraction) /25 (mV))
	rmd1b 	= Rmd1b * mg * (1e3) * exp((v-40+vshift) * valence * memb_fraction /25 (mV))
	rmd1u 	= Rmd1u * exp((-1)*(v-40+vshift) * valence * (1-memb_fraction) /25 (mV))
	rmd2b 	= Rmd2b * mg * (1e3) * exp((v-40+vshift) * valence * memb_fraction /25 (mV))
	rmd2u 	= Rmd2u * exp((-1)*(v-40+vshift) * valence * (1-memb_fraction) /25 (mV))

	~ U <-> Cl	(rb*qfac,Ru*qfac)
	~ Cl <-> Open	(Ro*qfac,Rc*qfac)
	~ Cl <-> D1	(Rd1*qfac,Rr1*qfac)
	~ D1 <-> D2	(Rd2*qfac,Rr2*qfac)
	~ Open <-> OMg	(rmb*qfac,rmu*qfac)
	~ UMg <-> ClMg 	(rbMg*qfac,RuMg*qfac)
	~ ClMg <-> OMg 	(RoMg*qfac,RcMg*qfac)
	~ ClMg <-> D1Mg (Rd1Mg*qfac,Rr1Mg*qfac)
	~ D1Mg <-> D2Mg (Rd2Mg*qfac,Rr2Mg*qfac)
	~ U <-> UMg     (rmc1b*qfac,rmc1u*qfac)
	~ Cl <-> ClMg	(rmc2b*qfac,rmc2u*qfac)
	~ D1 <-> D1Mg	(rmd1b*qfac,rmd1u*qfac)
	~ D2 <-> D2Mg	(rmd2b*qfac,rmd2u*qfac)

	CONSERVE U+Cl+D1+D2+Open+UMg+ClMg+D1Mg+D2Mg+OMg = 1
}
