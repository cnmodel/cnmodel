# -*- encoding: utf-8 -*-
from ._db import add_table_data

add_table_data('sgc_synapse', row_key='field', col_key='post_type', 
               species='mouse', data=u"""

AMPA_gmax and NMDA_gmax are the estimated average peak conductances (in nS) 
resulting from an action potential in a single auditory nerve terminal, under 
conditions that minimize the effects of short-term plasticity.

n_rsites is the number of release sites per SGC terminal.

------------------------------------------------------------------------------------------------
             bushy              tstellate        dstellate        pyramidal    octopus         tuberculoventral
                                                                  
AMPA_gmax    15.0±6.5 [1]       2.2±1.5 [2]      2.2±1.5 [7]                   0.87±0.23 [3]   2.2±1.5 [7]
NMDA_gmax    10.8±4.6 [1]       2.4±1.6 [2]      2.4±1.6 [7]                   0.17±0.046 [3]  2.4±1.6 [7]
EPSC_cv      0.12 [8]           0.20 [9]         0.27 [9]         
                                                                  
n_rsites     100 [5]            45 [6]           20 [6]                                        15 [8]

------------------------------------------------------------------------------------------------

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

[5] Oleskevich & Walmsley ~2002, Wang & Manis 2005

[6] Value was chosen to satisfy the CV of EPSC amplitude determined in [9]

[7] Data copied from t-stellate column (no literature on these cells)

[8] Thin air.

[9] Reanalysis of evoked EPSCs in stellate cells

""")



add_table_data('sgc_ampa_kinetics', row_key='field', col_key='post_type', 
               species='mouse', data=u"""
AMPA receptor kinetic values obtained by fitting the model of Raman and 
Trussell (1992) to measured EPSCs in the mouse VCN.

Ro1, Ro2, Rc1, Rc2, and PA are kinetic constants affecting the AMPA receptor
mechanism. tau_g and A affect the speed and amplitude of transmitter release
(implemented in the presynaptic release mechanism).
These parameters were selected to fit the model output to known EPSC shapes.

PA is a polyamine block parameter ued in the AMPAR mechanism.

------------------------------------------------------------------------------------------------
             bushy              tstellate        dstellate        pyramidal    octopus         tuberculoventral
                                                                  
Ro1          107.85 [4]         39.25 [4]        39.25 [7]                                     39.25 [7]
Ro2          0.6193 [4]         4.40 [4]         4.40 [7]                                      4.40 [7] 
Rc1          3.678 [4]          0.667 [4]        0.667 [7]                                     0.667 [7]
Rc2          0.3212 [4]         0.237 [4]        0.237 [7]                                     0.237 [7]
tau_g        0.10 [4]           0.25 [4]         0.25 [7]         
amp_g        0.770 [4]          1.56625 [4]      1.56625 [7]      
                                                                  
PA           45 [12]            0.1 [12]         0.1 [7]                                       0.1 [7]

------------------------------------------------------------------------------------------------

[4] Xie & Manis 2013, Table 2

[7] Data copied from t-stellate column (no literature on these cells)

[12] Wang & Manis (unpublished)

""")



add_table_data('sgc_epsp_kinetics', row_key='field', col_key='post_type', 
               species='mouse', data=u"""

EPSC shape parameters obtained from fits of Xie & Manis 2013 Equation 3 to measured EPSCs.

------------------------------------------------------------------------------------------------
             bushy              tstellate        dstellate        pyramidal    octopus         tuberculoventral
                                                                  
tau_r        0.253 [11]         0.19 [11]                          
tau_f        0.16 [11]          1.073 [11]                         
tau_s        0.765 [11]         3.3082 [11]                        
F            0.984 [11]         0.917 [11]                        
                                                 
------------------------------------------------------------------------------------------------

[11] Xie & Manis 2013, Table 3

""")



add_table_data('sgc_release_dynamics', row_key='field', col_key='post_type', 
               species='mouse', data=u"""

Kinetic parameters correspond to variables as described by Dittman et al. 
(2000), their Table 1.

F: Resting release probability

------------------------------------------------------------------------------------------------
             bushy             tstellate         dstellate        pyramidal    octopus         tuberculoventral
                                                                  
F            0.29366 [1]       0.43435 [1]       0.43435 [2]    
k0           0.52313e-3 [1]    0.06717e-3 [1]    0.06717e-3 [2] 
kmax         19.33805e-3 [1]   52.82713e-3 [1]   52.82713e-3 [2]
kd           0.11283 [1]       0.08209 [1]       0.08209 [2]    
ks           11.531 [1]        14.24460 [1]      14.24460 [2]   
kf           17.78 [1]         18.16292 [1]      18.16292 [2]   
taud         15.16 [1]         3.98 [1]          3.98 [2]       
taus         17912.2 [1]       16917.120 [1]     16917.120 [2]  
tauf         9.75 [1]          11.38 [1]         11.38 [2]      
dD           0.57771 [1]       2.46535 [1]       2.46535 [2]    
dF           0.60364 [1]       1.44543 [1]       1.44543 [2]

------------------------------------------------------------------------------------------------

[1] Xie & Manis 2013, Table 1

[2] Data copied from t-stellate column (no literature on these cells)

""")



add_table_data('gly_kinetics', row_key='field', col_key='post_type', 
               species='mouse', data=u"""

Kinetic parameters for glycine receptor mechanisms.

These are currently used for both DS and TV synapses, but should probably be 
separated in the future.

KV, KU, and XMax are kinetic parameters for the cleft transmitter mechanism.


------------------------------------------------------------------------------------------------
             bushy              tstellate        dstellate        pyramidal    octopus         tuberculoventral
                                                                  
KV           1e9                531.0            531.0
KU           4.46               4.17             4.17
XMax         0.733              0.731            0.731

------------------------------------------------------------------------------------------------

""")
