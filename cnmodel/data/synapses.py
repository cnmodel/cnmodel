# -*- encoding: utf-8 -*-
from ._db import add_table_data

# sgc old weights:
             # bushy             tstellate          dstellate         octopus         pyramidal      tuberculoventral
# weight       0.027 [12]        0.006 [12]         0.00064 [12]      0.0011 [12]     0.0023 [12]    0.0029 [12]
# tau1         0.1 [5]           0.1 [5]            0.2 [5]           0.1 [5]         0.1 [5]        0.1 [5]
# tau2         0.3 [5]           0.3 [5]            0.5 [5]           0.3 [5]         0.3 [5]        0.3 [5]
# erev         0   [5]           0   [5]            0   [5]           0   [5]         0   [5]        0   [5]

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

---------------------------------------------------------------------------------------------------------------------------------------
                  bushy             tstellate         dstellate         octopus           pyramidal         tuberculoventral  cartwheel                                                                                                                    
                                                                                                            
AMPA_gmax         21.05±15.4 [1]    4.6±3.1 [2]       0.49±0.29 [7]     0.87±0.23 [3]     0.6±0.3 [8]       2.2±1.5 [8]       0
AMPAR_gmax        4.6516398 [10]    4.632848  [10]    1.7587450 [10]    16.975147 [10]    0.9 [8]           2.2  [8]          0
NMDA_gmax         10.8±4.6 [1]      2.4±1.6 [2]       0.552±0.322 [7]   0.17±0.046 [3]    0.4±0.33 [8]      2.4±1.6 [8]       0
NMDAR_gmax        0.4531933 [10]    1.2127097 [10]    0.9960820 [10]    0.6562702 [10]    0.2 [8]           1.2127097 [8]     0
NMDAR_vshift      0.0   [12]        0.0   [12]        0.0   [12]        0.0   [12]        0.0   [12]        0.0   [12]        0
EPSC_cv           0.12 [8]          0.499759 [9]      0.886406 [9]      1.393382 [9]      0.499 [8]         0.499 [8]         0
Pr                1.000 [11]        1.000 [11]        1.000 [11]        1.000 [11]        1.000 [8]         1.000 [8]         0
n_rsites          100 [5]           4 [6]             1 [4]             1 [4]             2 [8]             2 [8]             0
delay             0.600             0.600             0.600             0.600             0.600             0.600             0
weight            0.020377          0.003679          0.000457          0.001311          0.000327          0.000808          0
tau1              0.158             0.174             0.152             0.125             0.167             0.157             0
tau2              0.246             1.501             1.652             0.251             1.489             1.641             0
erev              0.0               0.0               0.0               0.0               0.0               0.0               0
----------------------------------------------------------------------------------------------------------------------------------------

[1] Derived from Cao, X. & Oertel, D. (2010). Single-terminal conductance was
    reported as 21.5±15.4 nS (1.4±1.0 nA at -65 mV). The ratio of NMDA current to 
    total current is 0.3, so AMPA and NMDA currents are:
       AMPA_gmax = 21.5±15.4 nS (measured at -65 mV)
       NMDA_gmax = 21.5±15.4 nS * 0.3 = 10.8±4.6 nS
    Age>p17, Temperature=33C, [Mg2+]=1.3mM, [Ca2+]=2.4mM
    Units are nS.
    See also Pliss et al., J. Neurophys., 2009 (and note [12])

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

[7] (Xie and Manis, Frontiers in Neural Circuits, 2017):
    Measurements from CBA/CaJ mouse "radiate" multipolar cells in the AVCN.
    Single terminal conductance = (1.2 ± 0.70 nA/70 mV)/ 35 inputs = 0.490 ± 0.286 nS
    (see connections.py) 
    Single terminal conductance from mini = 34 pA/70 mV = 0.486 nS (single mini)
    Assume same AMPA/NMDA ratio as tstellate cells, but measures made where NMDA = 0
    (at negative V):
       AMPA_gmax = 0.490±0.286 nS
       NMDA_gmax = 0.490±0.286 nS * 0.53/0.47 = 0.552±0.322 nS
    Age > P35, Temperature=34C, [Mg2+]=1.5mM, [Ca2+]=2.5mM

[8] Thin air.  These are for testing the software, not necessarily for performing
    real simulations.  Note: Pyramidal cell strength has been reduced 
    because of large convergence and high input resistance of the reference cell model.
    Release 1 (Nov 2017):
     pyramidal    
                  
     0.6 ±1.05 [8]
     1.8 [8]      
     0.8±0.66 [8] 
     0.4 [8]      
     -15.0 [12]   
     0.499 [8]    
     1.000 [8]    
     2 [8]        
     
     
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

[12]  NMDA_vshift is the voltage shift for the activation of the NMDAR's, relative
      to 0 (standard in the NMDA_Kampa model). A negative value shifts the voltage
      dependence to the right (depolarizing).
      The value of the shift (0 or -15 mV) was chosen based on an exploration
      of fitting functions against the NMDA-Kampa IV curve in an SGC-bushy cell
      model, and comparing them against data. The functions were the modified
      Woodhull function and a Boltzmann function, yielding values of 1.19 mM for 
      k0 and 0.78 for delta (tau decay at +40 mV of 16.4 ms), and Vr -3 mV, Vh
      16 mV for the Boltzmann fit. These are close to the values reported in 
      for NMDA currents in p14-p26 CBA/CaJ mice in Pliss et al. (J. Neurophys. 
      102, 2627, 2009). Note: Pliss et al. agree with Cao and Oertel regarding
      an approximate 10-fold difference between AMPA and NMDA conductance in
      mouse bushy cells. An exact fit was not obtained, but no other parameters
      of the NMDA_Kampa model were changed. 
      ***
      Removed the follwing line 4/30/2019, as it was confusing.
      These were the -15 mV shifts. They cause the sgc->busy psd (for example) to fail, 
      because that was computed with a 0 mV scaling.
      NMDAR_vsh         -15.0 [12]        -15.0 [12]        -15.0 [12]        -15.0 [12]        -15.0 [12]        -15.0 [12]        0
      
   
[13]  weight is the weight to use in a netcon object (NEURON) for "simple"
      synapses based on the exp2syn mechanism.
      Parameters Weight, tau1, tau2, delay and erev from comare_simple_multisynapses 
      run and curve fitting (all cells)

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
             bushy              tstellate        dstellate        pyramidal    octopus         tuberculoventral       mso     
                                                                                                                                
Ro1          107.85 [4]         39.25 [4]        39.25 [7]        39.25 [4]    107.85 [5]      39.25 [7]              107.85 [4]
Ro2          0.6193 [4]         4.40 [4]         4.40 [7]         4.40 [4]     0.6193 [5]      4.40 [7]               0.6193 [4]
Rc1          3.678 [4]          0.667 [4]        0.667 [7]        0.667 [4]    3.678 [5]       0.667 [7]              3.678 [4] 
Rc2          0.3212 [4]         0.237 [4]        0.237 [7]        0.237 [4]    0.3212 [5]      0.237 [7]              0.3212 [4]
tau_g        0.10 [4]           0.25 [4]         0.25 [7]         0.25 [4]     0.10 [5]        0.25 [4]               0.10 [4]  
amp_g        0.770 [4]          1.56625 [4]      1.56625 [7]      1.56625 [4]  0.770 [5]       1.56625 [4]            0.770 [4] 
                                                                                                                                
PA           45 [12]            0.1 [12]         0.1 [7]          0.1 [12]     45 [5]          0.1 [7]                45 [12]   

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
                                                                          
tau_r        0.253 [11]    0.19 [11]                                      0.253 [13]
tau_f        0.16 [11]     1.073 [11]                                     0.16 [13]
tau_s        0.765 [11]    3.3082 [11]                                    0.765 [13]
F            0.984 [11]    0.917 [11]                                     0.984 [13]
                                            
------------------------------------------------------------------------------------------------

[11] Xie & Manis 2013, Table 3
[13] Copied from bushy cells; no direct data

""")


add_table_data('sgc_release_dynamics', row_key='field', col_key='post_type', 
               species='mouse', data=u"""

Kinetic parameters correspond to variables as described by Dittman et al. 
(2000), their Table 1.

F: ~ Resting release probability

---------------------------------------------------------------------------------------------------------------
             bushy             tstellate         dstellate        pyramidal     octopus        tuberculoventral
                                                                               
F            0.29366 [1]       0.43435 [1]       0.43435 [2]      0.43435 [1]   0.29366 [14]   0.43435 [1]   
k0           0.52313 [1]       0.06717 [1]       0.06717 [2]      0.06717 [1]   0.52313 [14]   0.06717 [1]   
kmax         19.33805 [1]      52.82713 [1]      52.82713 [2]     52.82713 [1]  19.33805 [14]  52.82713 [1]  
kd           0.11283 [1]       0.08209 [1]       0.08209 [2]      0.08209 [1]   0.11283 [14]   0.08209 [1]   
ks           11.531 [1]        14.24460 [1]      14.24460 [2]     14.24460 [1]  11.531 [14]    14.24460 [1]  
kf           17.78 [1]         18.16292 [1]      18.16292 [2]     18.16292 [1]  17.78 [14]     18.16292 [1]  
taud         15.16 [1]         3.98 [1]          3.98 [2]         3.98 [1]      15.16 [14]     3.98 [1]      
taus         17912.2 [1]       16917.120 [1]     16917.120 [2]    16917.120 [1] 17912.2 [14]   16917.120 [1] 
tauf         9.75 [1]          11.38 [1]         11.38 [2]        11.38 [1]     9.75 [14]      11.38 [1]     
dD           0.57771 [1]       2.46535 [1]       2.46535 [2]      2.46535 [1]   0.57771 [14]   2.46535 [1]   
dF           0.60364 [1]       1.44543 [1]       1.44543 [2]      1.44543 [1]   0.60364 [14]   1.44543 [1]   

---------------------------------------------------------------------------------------------------------------

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
                                                                            
KV           1e9 [1]      531.0 [1]        531.0 [1]     531.0 [2]         531.0 [2]
KU           4.46 [1]     4.17 [1]         4.17 [1]      4.17 [2]          4.17 [2] 
XMax         0.733 [1]    0.731 [1]        0.731 [1]     0.731 [2]         0.731 [2]

------------------------------------------------------------------------------------------------

[1] Xie & Manis 2013

[2] Copied from tstellate data (Kuo et al., J. Neurophysiol. indicate glycinergic IPSCs in TV
    and pyramidal cells are fast, with a decay time constant similar to that seen in tstellate
    cells). In pyramidal cells, this is consistent with the brief cross-correlation tip (Voigt
    and Young, 1980) and brief somatic current source (Manis and Brownell, 1983).


""")

add_table_data('dstellate_synapse', row_key='field', col_key='post_type', 
                species='mouse', data=u"""

DStellate Synapse values
gly_gmax is the default value in the program (scaled by Po for the receptors). See synapses/gly_psd.py
IPSC_cv is the coefficient of variation of the IPSC. (Not currently used in the model)
Pr is the release probabilty (not currently used); built into release mechanism for multisite synapses.
n_rsites is the number of release sites per dstellate terminal.

---------------------------------------------------------------------------------------------------------------------------------------
                  bushy             tstellate         dstellate         octopus           pyramidal         tuberculoventral  cartwheel
                                                                                                                              
gly_gmax          2.5 [1]           1.0 [5]           1.0 [2]           0. [2]            2.0 [3]           2.0 [3]           0±0 [2]
IPSC_cv           0.3 [3]           0.3 [3]           0.3 [3]           0.3 [3]           0.3 [3]           0.3 [3]           0.3 [3]    
Pr                1.000 [4]         1.000 [4]         1.000 [4]         1.000 [4]         1.000 [4]         1.000 [4]         1.000 [4]  
n_rsites          10 [5]            5 [5]             5 [5]             0 [2]             5 [5]             25 [5]            0 [2]      
delay             0.000             0.000             0.000             0                 0.000             0.000             0
weight            0.004131          0.004455          0.0               0                 0.002228          0.012097          0
tau1              0.187             0.152             0.152             0                 0.152             0.152             0
tau2              7.953             1.247             1.247             0                 1.247             1.247             0
erev              -70.0             -70.0             -70.0             0                 -70.0             -70.0             0
---------------------------------------------------------------------------------------------------------------------------------------

[1] Estimate

[2] No evidence for dstellate inputs to other d stellate cells or cartwheel cells.
    Octopus cells do not get inhibitory input
    
[3] Guess

[4] Default value

[5] Guess *educated* DS->TS from Xie and Manis, 2013. 99 pA mini @ 50 mV driving ~ 2 nS

[6] delay from pre to post; default is 0

[7] Parameters Weight, tau1, tau2, delay and erev from comare_simple_multisynapses run and curve fitting (all cells)

""")



add_table_data('tuberculoventral_synapse', row_key='field', col_key='post_type', 
                species='mouse', data=u"""

Tuberculventral Synapse values
gly_gmax is the default value in the program (scaled by Po for the receptors). See synapses/gly_psd.py
IPSC_cv is the coefficient of variation of the IPSC. (Not currently used in the model)
Pr is the release probabilty (not currently used)
n_rsites is the number of release sites per tuberculoventral terminal.

-----------------------------------------------------------------------------------------------------------------------------------
                  bushy             tstellate         dstellate         octopus           pyramidal         tuberculoventral  cartwheel
                                                                                                                              
gly_gmax          5.0 [3]           3.0 [3]           3.0 [3]           0. [2]            2.1±2.9 [6]       1.8±2.3 [6]       0±0 [6]
IPSC_cv           0.3 [3]           0.3 [3]           0.3 [3]           0.3 [3]           1.0 [3]           0.3 [3]           0.3 [3]    
Pr                1.000 [4]         1.000 [4]         1.000 [4]         1.000 [4]         1.000 [4]         1.000 [4]         1.000 [4]  
n_rsites          6 [5]             6 [5]             0 [1]             0 [2]             6 [5]             6 [5]             6 [5]      
delay             0.600             0.600             0                 0                 0.600             0.600             0
weight            0.002371          0.008114          0                 0                 0.002705          0.002705          0
tau1              0.190             0.149             0                 0                 0.149             0.149             0
tau2              7.952             1.250             0                 0                 1.250             1.250             0
erev              -70.0             -70.0             0                 0                 -70.0             -70.0             0
-----------------------------------------------------------------------------------------------------------------------------------

[1] Default value from GlyPSD

[2] No evidence for tuberculo inputs to other d stellate cells or cartwheel cells.
    Octopus cells do not get inhibitory input
    
[3] Guess

[4] Default value

[5] Guess

[6] Mouse data
    TV conductance onto pyr cells: 2.1 nS SD 2.9 nS (Kuo et al., 2012)
    TV conductance onto TV cells: 1.8 ns SD 2.3 nS.

[7] Parameters Weight, tau1, tau2, delay and erev from comare_simple_multisynapses run and curve fitting (all cells)
    Fitting done against 200 rep average for bushy, 500 rep average for all others.

""")

add_table_data('cartwheel_synapse', row_key='field', col_key='post_type', 
                species='mouse', data=u"""

Cartwheel cell synapse values
gly_gmax is the default value in the program (scaled by Po for the receptors). See synapses/gly_psd.py
IPSC_cv is the coefficient of variation of the IPSC. (Not currently used in the model)
Pr is the release probabilty (not currently used)
n_rsites is the number of release sites per cartwheel cell terminal.

-----------------------------------------------------------------------------------------------------------------------------------
             bushy             tstellate          dstellate         octopus         pyramidal      tuberculoventral  cartwheel
                                                                                                                
gly_gmax     0.0 [3]           0.0 [3]            0.0 [3]           0. [2]          2.1±2.9 [6]    0±0 [6]           1.8±2.3 [6]    
IPSC_cv      0.3 [3]           0.3 [3]            0.3 [3]           0.3 [3]         1.0 [3]        0.3 [3]           0.3 [3]       
Pr           1.000 [4]         1.000 [4]          1.000 [4]         1.000 [4]       1.000 [4]      1.000 [4]         1.000 [4]       
n_rsites     6 [5]             6 [5]              0 [1]             0 [2]           6 [5]          6 [5]             6 [5]            
delay        0 [7]             0                  0                 0               0              0                 0
weight       0.01              0.01               0.01              0.0             0.01           0.01              0.01 
tau1         0.3 [5]           0.3 [5]            0.3 [5]           0.3 [5]         0.3 [5]        0.3 [5]           0.3 [5]
tau2         2.0 [5]           2.0 [5]            2.0 [5]           2.0 [5]         2.0 [5]        2.0 [5]           2.0 [5]
erev         -70 [5]           -70 [5]            -70 [5]           -70 [5]         -70 [5]        -70 [5]           -70 [5]
-----------------------------------------------------------------------------------------------------------------------------------

[1] Default value from GlyPSD

[2] No evidence for cartwheel inputs to Dstellate, bushy or tstellate cells.
    Octopus cells do not get inhibitory input
    
[3] Guess

[4] Default value

[5] Guess

[6] Mouse data
    TV conductance onto pyr cells: 2.1 nS SD 2.9 nS (Kuo et al., 2012)
    TV conductance onto TV cells: 1.8 ns SD 2.3 nS.

[7] delay from pre to post; default is just 0

""")

add_table_data('bushy_synapse', row_key='field', col_key='post_type', 
               species='mouse', data=u"""

AMPA_gmax and NMDA_gmax are the estimated average peak conductances (in nS) 
resulting from an action potential in a single presynaptic terminal under 
conditions that minimize the effects of short-term plasticity.
AMPA_gmax are from values measured at -65 mV (or -70mV), and represent SINGLE TERMINAL 
conductances
AMPAR_gmax are the individual synapse postsynaptic conductance
NMDA_gmax values are taken as the fraction of the current that is NMDAR dependent
at +40 mV (see below)

n_rsites is the number of release sites per terminal.

-----------------------------------------------------------------------------------------------------------------------------------
             mso           
                          
AMPA_gmax    21.05±15.4 [1]
AMPAR_gmax   4.6516398 [2]
NMDA_gmax    0 [3]  
NMDAR_gmax   0 [3]
EPSC_cv      0.12 [4]      
Pr           1.000 [5]    
n_rsites     36 [6]         
weight       0.01
delay        0
-----------------------------------------------------------------------------------------------------------------------------------

[1] Taken from the mouse bushy cell model.
    Units are nS.
    
[2] See note [10] for the SGC-bushy synapse

[3] Assume no NMDA receptors at this synapse

[4] See SGC-bushy synapse

[5] Just to scale with the multisite synapse model

[6] This is a guess.

""")

