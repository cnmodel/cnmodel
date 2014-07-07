TITLE Cerebellum Golgi Cell Model

COMMENT
        Sodium first order kinetics
   
	Author: F. Simoes de Souza
	Revised: 15.07.09
ENDCOMMENT

NEURON {
        SUFFIX na_conc
        USEION na READ ina, nao, nai WRITE nai, nao
        RANGE d, beta, nai0, nao0
}

UNITS {
        (mV)    = (millivolt)
        (mA)    = (milliamp)
	(um)    = (micron)
	(molar) = (1/liter)
        (mM)    = (millimolar)
   	FARADAY      = (faraday) (coulomb)
}

PARAMETER {
        celsius 	 (degC)
        d = .2           (um)
        nao0 = 145        (mM)         
        nai0 = 5         (mM)         
        beta = 1.3       (/ms)
	ina              (mA/cm2)
}

ASSIGNED {

}
STATE {
	nai (mM)
	nao (mM)
}

INITIAL {
        nai = nai0 
        nao = nao0 
}

BREAKPOINT {
       SOLVE conc METHOD derivimplicit
}

DERIVATIVE conc {    
        nai' =  -(ina)/(2*FARADAY*d)*(1e4) - beta*(nai-nai0)
        nao' =  (ina)/(2*FARADAY*d)*(1e4) - beta*(nao-nao0)
	
}
