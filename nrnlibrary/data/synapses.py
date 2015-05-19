# -*- encoding: utf-8 -*-
from ._db import add_table_data

add_table_data('sgc_synapse', row_key='field', col_key='post_type', 
               species='mouse', data=u"""

AMPA_gmax and NMDA_gmax are the estimated average peak conductances (in nS) 
resulting from an action potential in a single auditory nerve terminal, under 
conditions that minimize the effects of short-term plasticity.

Ro1, Ro2, Rc1, Rc2, and PA are kinetic constants affecting the AMPA receptor
mechanism.

-------------------------------------------------------------------------------
             bushy        tstellate   dstellate   pyramidal    octopus         tuberculoventral

AMPA_gmax    15.0±6.5 [1] 2.2±1.5 [2] 2.2±1.5 [7]              0.87±0.23 [3]   2.2±1.5 [7]
NMDA_gmax    10.8±4.6 [1] 2.4±1.6 [2] 2.4±1.6 [7]              0.17±0.046 [3]  2.4±1.6 [7]
EPSC_cv      0.17 [8]     0.20 [9]    0.27 [9]


Ro1          107.85 [4]   39.25 [4]   39.25 [7]                                39.25 [7]
Ro2          0.6193 [4]   4.40 [4]    4.40 [7]                                 4.40 [7] 
Rc1          3.678 [4]    0.667 [4]   0.667 [7]                                0.667 [7]
Rc2          0.3212 [4]   0.237 [4]   0.237 [7]                                0.237 [7]
PA                        0.1 [4]     0.1 [7]                                  0.1 [7]

n_rsites     100 [5]      50 [6]      15 [8]                                   15 [8]

-------------------------------------------------------------------------------

[1] Derived from Cao, X. & Oertel, D. (2010). Single-terminal conductance was
    reported as 21.5±15.4 nS. The ratio of NMDA current to 
    total current is 0.3, so AMPA and NMDA currents are:
       AMPA_gmax = 21.5±15.4 nS * 0.7 = 15±6.5 nS
       NMDA_gmax = 21.5±15.4 nS * 0.3 = 10.8±4.6 nS
    Age>p17, Temperature=33C, [Mg2+]=1.3mM, [Ca2+]=2.4mM

[2] Derived from Cao, X. & Oertel, D. (2010). Single-terminal conductance was
    estimated as 4.6±3.1 nS. The ratio of NMDA current to 
    total current is 0.53, so AMPA and NMDA currents are:
       AMPA_gmax = 4.6±3.1 nS * 0.47 = 2.2±1.5 nS
       NMDA_gmax = 4.6±3.1 nS * 0.53 = 2.4±1.6 nS
    Age>p17, Temperature=33C, [Mg2+]=1.3mM, [Ca2+]=2.4mM

[3] Derived from Cao, X. & Oertel, D. (2010). Single-terminal conductance was
    estimated as 52±14 nS / 60 = 0.87±0.23 nS. The ratio of NMDA current to 
    total current is 0.2, so AMPA and NMDA currents are:
       AMPA_gmax = 0.87±0.23 nS * 0.8 = 0.70±0.18 nS
       NMDA_gmax = 0.87±0.23 nS * 0.2 = 0.17±0.046 nS
    Age>p17, Temperature=33C, [Mg2+]=1.3mM, [Ca2+]=2.4mM

[4] Xie & Manis (2013)

[5] Oleskevich & Walmsley ~2002, Wang & Manis 2005

[6] Value was chosen to satisfy the CV of EPSC amplitude determined in [9]

[7] Data copied from t-stellate column (no literature on these cells)

[8] Thin air.

[9] Reanalysis of evoked EPSCs in stellate cells

[10] Lu, Harris, & Rubel 2007

""")
