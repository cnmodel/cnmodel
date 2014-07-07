: Na diffusion in a cell
: based on cadiff from: Created 8/15/02 - nwg

NEURON {
       SUFFIX nadiff
       USEION na READ ina, nai WRITE nai
       RANGE na
       GLOBAL depth, beta
}

UNITS {
      (mV) = (millivolt)
      (mA) = (milliamp)
      (mM) = (milli/liter)
      (um) = (micron)
}

CONSTANT {
      F = 9.6485e4 (coul)
}

PARAMETER {
          nai      (mM)
          dt       (ms)

          depth  = 0.1  (um)
          beta = 1 (/ms)
}

ASSIGNED {
         ina       (mA/cm2)
}

STATE {
      na           (mM)
}

INITIAL {
        na = 0.004
}

BREAKPOINT {
        na = na + (10000.0) * dt * ( ( -1/(2*F)*ina / (depth)) - (0.0001) * beta * na )

        if ( na < 1e-5 ) {: minimum 10 nM Na
           na = 1e-4
        }

        nai = na
}
