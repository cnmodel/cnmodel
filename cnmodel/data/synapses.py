# -*- encoding: utf-8 -*-
from ._db import add_table_data

add_table_data('sgc_synapse', row_key='field', col_key='post_type', 
               species='mouse', data=u"""

AMPA_gmax and NMDA_gmax are the estimated average peak conductances (in nS) 
resulting from an action potential in a single auditory nerve terminal, under 
conditions that minimize the effects of short-term plasticity.
AMPA_gmax are from values measured at -65 mV (or -70mV), and represent SINGLE TERMINAL 
conductances
AMPAR_gmax are the individual synapse postsynaptic conductance
NMDA_gmax values are taken as the fraction of the current that is NMDAR dependent
at +40 mV (see below)

n_rsites is the number of release sites per SGC terminal.

-----------------------------------------------------------------------------------------------------------------------------------
             bushy             tstellate          dstellate         octopus        pyramidal      tuberculoventral
                                                           
AMPA_gmax    21.05±15.4 [1]    4.6±3.1 [2]        0.49±0.29 [7]     0.87±0.23 [3]                 2.2±1.5 [7]
AMPAR_gmax   4.6516398 [10]    4.632848  [10]     1.7587450 [10]    16.975147 [10]
NMDA_gmax    10.8±4.6 [1]      2.4±1.6 [2]        0.552±0.322 [7]   0.17±0.046 [3]                2.4±1.6 [7]
NMDAR_gmax   0.4531933 [10]    1.2127097 [10]     0.9960820 [10]    0.6562702 [10]
EPSC_cv      0.12 [8]          0.499759 [9]       0.886406 [9]      1.393382 [9]
Pr           1.000 [11]        1.000 [11]         1.000 [11]        1.000 [11]
n_rsites     100 [5]           4 [6]              1 [4]             1 [4]                         15 [8]

-----------------------------------------------------------------------------------------------------------------------------------

[1] Derived from Cao, X. & Oertel, D. (2010). Single-terminal conductance was
    reported as 21.5±15.4 nS (1.4±1.0 nA at -65 mV). The ratio of NMDA current to 
    total current is 0.3, so AMPA and NMDA currents are:
       AMPA_gmax = 21.5±15.4 nS (measured at -65 mV)
       NMDA_gmax = 21.5±15.4 nS * 0.3 = 10.8±4.6 nS
    Age>p17, Temperature=33C, [Mg2+]=1.3mM, [Ca2+]=2.4mM
    Units are nS.

[2] Derived from Cao, X. & Oertel, D. (2010). Single-terminal conductance was
    estimated as 4.6±3.1 nS. The ratio of NMDA current to 
    total current is 0.53, so AMPA and NMDA currents are:
       AMPA_gmax = 4.6±3.1 nS
       NMDA_gmax = 4.6±3.1 nS * 0.53 = 2.4±1.6 nS
    Estimated number of inputs per AN fiber:
        0.3 nA step, 0.08 nA mini size = ~ 4 inputs per AN fiber
    Age>p17, Temperature=33C, [Mg2+]=1.3mM, [Ca2+]=2.4mM
    Units are nS

[3] Derived from Cao, X. & Oertel, D. (2010). Single-terminal conductance was
    estimated as 52±14 nS / 60 = 0.87±0.23 nS. The ratio of NMDA current to 
    total current is 0.2, so AMPA and NMDA currents are:
       AMPA_gmax = 0.87±0.23 nS
       NMDA_gmax = 0.87±0.23 nS * 0.2 = 0.17±0.046 nS
    Age>p17, Temperature=33C, [Mg2+]=1.3mM, [Ca2+]=2.4mM
    Units are nS

[4] Assumption based on mini size and lack of discernable EPSC step (guess).
    Should be verified.

[5] Oleskevich & Walmsley ~2002, Wang & Manis 2005. Units are nS

[6] A value of 45 would be chosen to satisfy the CV of EPSC amplitude determined in [9].
    However, those measures are for simultaneous stimulation of multiple AN fibers.
    A value of 4 is included here to correspond to measures in Cao and Oertel (2010)
    (see note [2])

[7] (Xie and Manis, unpublished, 2017):
    Measurements from CBA/CaJ mouse "radiate" multipolar cells in the AVCN.
    Single terminal conductance = (1.2 ± 0.70 nA/70 mV)/ 35 inputs = 0.490 ± 0.286 nS
    (see connections.py) 
    Single terminal conductance from mini = 34 pA/70 mV = 0.486 nS (single mini)
    Assume same AMPA/NMDA ratio as tstellate cells, but measures made where NMDA = 0
    (at negative V):
       AMPA_gmax = 0.490±0.286 nS
       NMDA_gmax = 0.490±0.286 nS * 0.53/0.47 = 0.552±0.322 nS
    Age > P35, Temperature=34C, [Mg2+]=1.5mM, [Ca2+]=2.5mM

[8] Thin air.

[9] Reanalysis of evoked EPSCs in stellate cells (Manis/Xie, 2014)

[10]  Maximum AMPA open conductance per synaptic site (units are pS). 
      These values are calculated by running python cnmodel/synapses/tests/test_psd.py
      for a specific cell type (if the cell uses the receptor mechanisms; this is 
      not necessary for simple exp2syn style mechanisms)
      to ensure that maximum AMPA conductance during PSG matches [1, 2 or 3]
      For a bushy cell, the original default values  (bushy cell) were:
          AMPAR_gmax   3.314707700918133
          NMDAR_gmax   0.4531929783503451
      These values will also depend on the number of release sites per
      synapse (the total conductance is produce of site gmax and nsites).
      
      A note on the precision of these values: This precision is only
      required for the tests of the model, as a way of ensuring numerical
      equivalency after potential modifications of the code. The precision
      of the value is in no way intended to specificy biological precision.

      For example, a change in the rate constants in the AMPA_Trussell AMPA
      receptor model could (and probably would) change the open probability,
      and therefore the maximal conductance of an EPSC. However, as this is
      only a representation of the EPSC, the "receptor" conductance should
      be scaled so that the computed EPSC has the same maximal conductance
      as prior to the kinetic modifications. Because the receptor model is
      numerically computed (and not analytically tractable without
      additional knowledge of the ligand time course), a numerical solution
      is required.

[11]  Pr is the initial release probability. The value can be computed by
      setting Pr to 1 in this file, and running the cnmodel test_synapses.py
      with the appropriate presynaptic source and postsynaptic target,
      once all other parameters are set. The Pr is used to rescale
      the AMPAR_gmax so that the total current matches the data in 
      AMPA_gmax in the table (on average).
""")



add_table_data('sgc_ampa_kinetics', row_key='field', col_key='post_type', 
               species='mouse', data=u"""
AMPA receptor kinetic values obtained by fitting the model of Raman and 
Trussell (1992) to measured EPSCs in the mouse VCN.

Ro1, Ro2, Rc1, Rc2, and PA are kinetic constants affecting the AMPA receptor
mechanism. tau_g and A affect the speed and amplitude of transmitter release
(implemented in the presynaptic release mechanism).
These parameters were selected to fit the model output to known EPSC shapes.

PA is a polyamine block parameter ued in the AMPAR mechanism (concentration in micromolar).

------------------------------------------------------------------------------------------------
             bushy              tstellate        dstellate        pyramidal    octopus         tuberculoventral
                                                                  
Ro1          107.85 [4]         39.25 [4]        39.25 [7]                     107.85 [5]      39.25 [7]
Ro2          0.6193 [4]         4.40 [4]         4.40 [7]                      0.6193 [5]      4.40 [7] 
Rc1          3.678 [4]          0.667 [4]        0.667 [7]                     3.678 [5]       0.667 [7]
Rc2          0.3212 [4]         0.237 [4]        0.237 [7]                     0.3212 [5]      0.237 [7]
tau_g        0.10 [4]           0.25 [4]         0.25 [7]                      0.10 [5]  
amp_g        0.770 [4]          1.56625 [4]      1.56625 [7]                   0.770 [5] 
                                                                                         
PA           45 [12]            0.1 [12]         0.1 [7]                       45 [5]          0.1 [7]

------------------------------------------------------------------------------------------------

[4] Xie & Manis 2013, Table 2

[5] copied from bushy cells; no direct data.

[7] Data copied from t-stellate column (no literature on these cells). Unpublished data suggests these
    should be slightly different, but is complicated by electrotonically distant synaptic sites that
    preclude accurate measurement of kinetics.

[12] Wang & Manis (unpublished)

""")



add_table_data('sgc_epsp_kinetics', row_key='field', col_key='post_type', 
               species='mouse', data=u"""

EPSC shape parameters obtained from fits of Xie & Manis 2013 Equation 3 to measured EPSCs.

------------------------------------------------------------------------------------------------
             bushy         tstellate        dstellate        pyramidal    octopus         tuberculoventral
                                                                          0.253 [13]
tau_r        0.253 [11]    0.19 [11]                                      0.16 [13]
tau_f        0.16 [11]     1.073 [11]                                     0.765 [13]
tau_s        0.765 [11]    3.3082 [11]                                    0.984 [13]
F            0.984 [11]    0.917 [11]                        
                                            
------------------------------------------------------------------------------------------------

[11] Xie & Manis 2013, Table 3
[13] Copied from bushy cells; no direct data

""")



add_table_data('sgc_release_dynamics', row_key='field', col_key='post_type', 
               species='mouse', data=u"""

Kinetic parameters correspond to variables as described by Dittman et al. 
(2000), their Table 1.

F: ~ Resting release probability

------------------------------------------------------------------------------------------------
             bushy             tstellate         dstellate        pyramidal    octopus         tuberculoventral
                                                                               
F            0.29366 [1]       0.43435 [1]       0.43435 [2]                   0.29366 [14] 
k0           0.52313 [1]       0.06717 [1]       0.06717 [2]                   0.52313 [14] 
kmax         19.33805 [1]      52.82713 [1]      52.82713 [2]                  19.33805 [14]
kd           0.11283 [1]       0.08209 [1]       0.08209 [2]                   0.11283 [14] 
ks           11.531 [1]        14.24460 [1]      14.24460 [2]                  11.531 [14]  
kf           17.78 [1]         18.16292 [1]      18.16292 [2]                  17.78 [14]   
taud         15.16 [1]         3.98 [1]          3.98 [2]                      15.16 [14]   
taus         17912.2 [1]       16917.120 [1]     16917.120 [2]                 17912.2 [14] 
tauf         9.75 [1]          11.38 [1]         11.38 [2]                     9.75 [14]    
dD           0.57771 [1]       2.46535 [1]       2.46535 [2]                   0.57771 [14] 
dF           0.60364 [1]       1.44543 [1]       1.44543 [2]                   0.60364 [14] 

------------------------------------------------------------------------------------------------

[1] Xie & Manis 2013, Table 1. Although independently measured in > P30 CBA/CaJ mice,
    the values are similar to the measurements from Yang and Xu-Friedman, 2008
    in P14-P21 CBA/CaJ mice.

[2] Data copied from t-stellate column (no literature on these cells)

[14] Data copied from bushy cell column (no literature on these cells)
""")


add_table_data('gly_kinetics', row_key='field', col_key='post_type', 
               species='mouse', data=u"""

Kinetic parameters for glycine receptor mechanisms.

These are currently used for both DS and TV synapses, but should probably be 
separated in the future.

KV, KU, and XMax are kinetic parameters for the cleft transmitter mechanism.


------------------------------------------------------------------------------------------------
             bushy        tstellate        dstellate     pyramidal         tuberculoventral
                                                         
KV           1e9 [1]      531.0 [1]        531.0 [1]     531.0 [2]
KU           4.46 [1]     4.17 [1]         4.17 [1]      4.17 [2] 
XMax         0.733 [1]    0.731 [1]        0.731 [1]     0.731 [2]

------------------------------------------------------------------------------------------------

[1] Xie & Manis 2013

[2] Copied from tstellate data (Kuo et al., J. Neurophysiol. indicate glycinergic IPSCs in TV
    and pyramidal cells are fast, with a decay time constant similar to that seen in tstellate
    cells). In pyramidal cells, this is consistent with the brief cross-correlation tip (Voigt
    and Young, 1980) and brief somatic current source (Manis and Brownell, 1983).


""")


# Mouse data
# TV conductance onto pyr cells: 2.1 nS SD 2.9 nS (Kuo et al., 2012)
# TV conductance onto TV cells: 1.8 ns SD 2.3 nS.
#

