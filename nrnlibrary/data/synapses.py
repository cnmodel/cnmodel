# -*- encoding: utf-8 -*-
from ._db import add_table_data

add_table_data('sgc_synapse', row_key='field', col_key='post_type', 
               species='mouse', data=u"""

AMPA_gmax and NMDA_gmax are the estimated average peak conductances (in nS) 
resulting from an action potential in a single auditory nerve terminal, under 
conditions that minimize the effects of short-term plasticity.

Ro1, Ro2, Rc1, and Rc2 are kinetic constants affecting the AMPA receptor
mechanism.

-------------------------------------------------------------------------------
             bushy        tstellate   dstellate  pyramidal    octopus         tuberculoventral
AMPA_gmax    15±6.5 [1]   2.2±1.5 [2]                         0.87±0.23 [3]
NMDA_gmax    10.8±4.6 [1] 2.4±1.6 [2]                         0.17±0.046 [3]
Ro1          107.85 [4]
Ro2          0.6193 [4]
Rc1          3.678 [4]
Rc2          0.3212 [4]



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

""")
