TITLE NaKPump adapted from Kyoto's model

COMMENT

Author: Fabio M Simoes de Souza

The equations that simulated the Na+/K+ ATPase are described in Table S10 of Takeuchi et al. (2006)

Reference:
Takeuchi A, Tatsumi S, Sarai N, Terashima K, Matsuoka S, Noma A (2006) Ionic mechanisms of cardiac cell swelling induced by blocking Na+/K+ pump as revealed by experiments and simulation. J. G. Physiol. 128: 495-507.

ENDCOMMENT

:********************************************/

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

NEURON {
	SUFFIX nakpump
	USEION na READ nai, nao WRITE ina
	USEION k READ ki,ko WRITE ik
	RANGE Nai_inf, Ki_inf, ATPi 
	RANGE imax, nai, ki
	RANGE y,inak, inakmax
}

PARAMETER {
	celsius   (degC)
	inakmax=1 (mA/cm2)
	R=8.314472 (j/(K.molar)
	F=96.4853 (C/mol)	
	nao=145   (mM)
	ko=5 	  (mM)
        Nai_inf=5 (mM)
	Ki_inf=14 (mM)
        ouabain=1 (mM)
	KdNao=69.8 
	KdKo=0.258 
	KdNai=4.05 
	KdKi=32.88 
	ATPi=6.97 (mM)
	k2=0.04    
	k3=0.01    
	k4=0.165   
	y0=0
}

ASSIGNED {
	v	(mV)
	inak	(mA/cm2)
	ina	(mA/cm2)
	ik	(mA/cm2)
	nai 	(mM)
	ki	(mM)
	pE1Na	
	pE1K	
	pE2Na	
	pE2K	
	alfay	
	betay	
	Naeff	(mM)
        k1	
	drugblock 
	kelvin	(degK)
}

STATE {
	y
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	settables(nai, ki, v, celsius, ATPi,ouabain)
	inak=inakmax*0.001*11.5*(k1*pE1Na*y-k2*pE2Na*(1-y))*drugblock
    ina=3*inak
	ik =-2*inak
}

INITIAL {
	settables(Nai_inf, Ki_inf, v, celsius, ATPi,ouabain)
	y=.6224
}

DERIVATIVE states {
	settables(nai,ki,v, celsius, ATPi,ouabain)
	y'=alfay*(1-y)-betay*y
}

UNITSOFF

PROCEDURE settables(nai, ki, v, celsius, ATPi) {
	drugblock=1/(1+ouabain/0.0006)
        k1=0.37*(1/(1+0.094/ATPi))  
	kelvin=celsius+273.16  
	Naeff=nao*exp(-0.82*F*v/R/kelvin)
	alfay=k2*pE2Na+k4*pE2K
	betay=k1*pE1Na+k3*pE1K	
	pE1Na=1/(1+(KdNai/nai)^1.06*(1+(ki/KdKi)^1.12))
	pE1K=1/(1+(KdKi/ki)^1.12*(1+(nai/KdNai)^1.06))
	pE2Na=1/(1+(KdNao/Naeff)^1.06*(1+(ko/KdKo)^1.12))
	pE2K=1/(1+(KdKo/ko)^1.12*(1+(Naeff/KdNao)^1.06))
}

UNITSON
