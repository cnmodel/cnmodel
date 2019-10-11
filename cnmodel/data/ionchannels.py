# -*- encoding: utf-8 -*-
from ._db import add_table_data
"""
Ion channel density tables
All of the ion channel densities for the models implemented in cnmodel
are (or should be) stated here, and should not be modified in the
cnmodel code itself. 

"""

add_table_data('RM03_channels', row_key='field', col_key='model_type', 
               species='guineapig', data=u"""

This table describes the ion channel densities (and voltage shifts if necessary)
for different cell types in the original Rothman Manis 2003 model.
Data from Table 1, except for "octopus" cells, which is modified (see note 3)
map to cell:    bushy-II      bushy-II-I    tstellate     tstellate-t   bushy-I-II    octopus
-----------------------------------------------------------------------------------------------------------------------------------
               II            II-I          I-c           I-t           I-II          II-o

nacn_gbar      1000. [1]     1000. [1]     1000. [1]     1000. [1]     1000. [2]     0000. [3]
jsrna_gbar     0000. [1]     0000. [1]     0000. [1]     0000. [1]     0000. [2]     1000. [3]
kht_gbar       150.0 [1]     150.0 [1]     150.0 [1]     80.0  [1]     150.0 [2]     150.0 [3] 
klt_gbar       200.0 [1]     35.0  [1]     0.0   [1]     0.0   [1]     20.0  [2]     1000. [3] 
ka_gbar        0.0   [1]     0.0   [1]     0.0   [1]     65.0  [1]     0.0   [2]     0.0   [3]
ih_gbar        20.0  [1]     3.5   [1]     0.5   [1]     0.5   [1]     2.0   [2]     30.0  [3]
leak_gbar      2.0   [1]     2.0   [1]     2.0   [1]     2.0   [1]     2.0   [2]     2.0   [3]
leak_erev      -65   [1]     -65   [1]     -65   [1]     -65   [1]     -65   [2]     -65   [3]
na_type        nacn  [1]     nacn  [1]     nacn  [1]     nacn  [1]     nacn  [2]     jsrna [3]
ih_type        ihvcn [1]     ihvcn [1]     ihvcn [1]     ihvcn [1]     ihvcn [2]     ihvcn [3]
soma_Cap       12.0  [1]     12.0  [1]     12.0  [1]     12.0  [1]     12.0  [2]     25.0  [3]
e_k            -84   [1]     -84   [1]     -84   [1]     -84   [2]     -84   [2]     -84   [2] 
e_na           50.   [1]     50.   [1]     50.   [1]     50.   [2]     50.   [2]     50.   [2] 
ih_eh          -43   [1]     -43   [1]     -43   [1]     -43   [2]     -43   [2]     -43   [2] 

-----------------------------------------------------------------------------------------------------------------------------------

[1] Rothman and Manis, 2003
    Age "adult", Temperature=22C
    Units are nS.

[2] Rothman and manis, 2003, model I-II
    Some low-voltage K current, based on observations of
    a single spike near threshold and regular firing for higher
    currents (Xie and Manis, 2017)
    
[3] Derived from Rothman and Manis, 2003, model II
    Large amounts of low-voltage K current, and elevated HCN. Conductances
    based on Rothman and Manis, 2003; concept from Cao and Oertel

[4] Designation for elevated LTK and Ih for octopus cells

""")

add_table_data('XM13_channels', row_key='field', col_key='model_type', 
               species='mouse', data=u"""

This table describes the REFERENCE ion channel densities (and voltage shifts if necessary)
for different cell types based on the Xie and Manis 2013 models for mouse.

The REFERENCE values are applied to "point" models, and to the soma of
compartmental models.
The names of the mechanisms must match a channel mechanism (Neuron .mod files)
and the following _(gbar, vshift, etc) must match an attribute of that channel
that can be accessed.

-----------------------------------------------------------------------------------------------------------------------------------
               II             II-I           I-c           I-II          I-t       
                                                                                   
nav11_gbar     1000.  [4]     0000.  [4]     800.   [4]    800.   [4]    1000.  [4] 
nacn_gbar      2300.  [1]     1000.  [1]     3000.  [1]    0000.  [2]    0000.  [1] 
na_gbar        0000.  [1]     0000.  [1]     3000.  [1]    1800.  [2]    0000.  [1] 
kht_gbar       58.0   [1]     58.0   [1]     500.0  [1]    150.0  [2]    500.0  [1] 
klt_gbar       80.0   [1]     20.0   [1]     0.0    [1]    14.0   [3]    0.0    [1] 
ka_gbar        0.0    [1]     0.0    [1]     0.0    [1]    0.0    [2]    125.0  [1] 
ihvcn_gbar     30.0   [1]     30.0   [1]     18.0   [1]    2.0    [2]    18.0   [1] 
leak_gbar      2.0    [1]     2.0    [1]     8.0    [1]    2.0    [2]    8.0    [1] 
leak_erev      -65    [1]     -65    [1]     -65    [1]    -65    [2]    -65    [1] 
na_type        nacn   [1]     nacn   [1]     nacn   [1]    na     [3]    nav11  [1] 
ih_type        ihvcn  [1]     ihvcn  [1]     ihvcn  [1]    ihvcn  [2]    ihvcn  [1] 
soma_Cap       26.0   [1]     26.0   [1]     25.0   [1]    25.0   [2]    25.0   [1] 
nav11_vshift   4.3    [1]     4.3    [1]     4.3    [1]    4.3    [1]    4.3    [1]
e_k            -84    [1]     -84    [1]     -84    [1]    -70    [3]    -84    [1] 
e_na           50.    [1]     50.    [1]     50.    [1]    55.    [3]    50.    [1] 
ih_eh          -43    [1]     -43    [1]     -43    [1]    -43    [2]    -43    [1] 

-----------------------------------------------------------------------------------------------------------------------------------

[1] Uses channels from Rothman and Manis, 2003
    Conductances are for Mouse bushy cells
    Xie and Manis, 2013
    Age "adult", Temperature=34C
    Units are nS.

[2] Rothman and manis, 2003, model I-II
    Some low-voltage K current, based on observations of
    a single spike near threshold and regular firing for higher
    currents (Xie and Manis, 2017)

[3] These values for the I-II (dstellate) are from the original checkpoint test
    for cnmodel 12/2017. 

[4] nav11 channels were used in original Xie and Manis (2013) ms, but are not
    used for mice in the master distribution of cnmodel, which used only the nacn
    channels. The channel type can be overridden however.

""")

add_table_data('XM13_channels_compartments', row_key='parameter', col_key='compartment', 
               species='mouse', model_type='II', data=u"""

This table describes the ion channel densities relative to somatic densities,
e.g., relative to REFERENCE densities in the table XM13_channels.
and voltage shifts, for different compartments of the specified neuron,
Conductances will be calculated from the Model derived from Xie and Manis 2013 for mouse
(data table: mGVC_channels).

------------------------------------------------------------------------------------------------------------------------------------------------------------------
               axon       unmyelinatedaxon     myelinatedaxon     initialsegment    hillock     soma        dendrite         primarydendrite    secondarydendrite
                                                                                                                                                                                              
nav11_gbar     3.0 [1]    3.0 [1]              0.0 [1]            5.0 [1]           5.0 [1]     1.0 [1]     0.5 [1]          0.50 [1]           0.25 [1]       
kht_gbar       1.0 [1]    2.0 [1]              0.01 [1]           2.0 [1]           2.0 [1]     1.0 [1]     0.5 [1]          0.5 [1]            0.25 [1]       
klt_gbar       1.0 [1]    1.0 [1]              0.01 [1]           1.0 [1]           1.0 [1]     1.0 [1]     0.5 [1]          0.5 [1]            0.25 [1]       
ihvcn_gbar     0.0 [1]    0.0 [1]              0.0 [1]            0.5 [1]           0.0 [1]     1.0 [1]     0.5 [1]          0.5 [1]            0.5 [1]       
leak_gbar      1.0 [1]    0.25 [1]             0.25e-3 [1]        1.0 [1]           1.0 [1]     1.0 [1]     0.5 [1]          0.5 [1]            0.5 [1]       
leak_erev      -65. [1]   -65. [1]             -65. [1]           -65. [1]          -65. [1]    -65. [1]    -65. [1]         -65. [1]           -65. [1]      
nav11_vshift   4.3  [1]   4.3  [1]             0.0 [1]            4.3  [1]          4.3  [1]    0.0 [1]     0.0 [1]          0.0 [1]            0.0 [1]       
na_type        nav11      nav11                nav11              nav11             nav11       nav11       nav11            nav11              nav11
ih_type        ihvcn      ihvcn                ihvcn              ihvcn             ihvcn       ihvcn       ihvcn            ihvcn              ihvcn                            
-------------------------------------------------------------------------------------------------------------------------------------------------------------------

[1] Scaling is relative to soma scaling. Numbers are estimates based on general distribution from literature on cortical neurons.


""")



# ***** BEGINNING OF XM13_Channels for nacncoop version of model


add_table_data('XM13nacncoop_channels', row_key='field', col_key='model_type', 
               species='mouse', data=u"""

This table describes the REFERENCE ion channel densities (and voltage shifts if necessary)
for different cell types based on the Xie and Manis 2013 models for mouse, but using
the nacncoop mechanism (coooperative sodium channels)

!!!!!!!!!!!! USAGE OF THIS TABLE SHOULD BE CONSIDERED EXPERIMENTAL !!!!!!!!!!!!!!

The REFERENCE values are applied to "point" models, and to the soma of
compartmental models.
The names of the mechanisms must match a channel mechanism (Neuron .mod files)
and the following _(gbar, vshift, etc) must match an attribute of that channel
that can be accessed.

-----------------------------------------------------------------------------------------------------------------------------------
                 II             II-I           I-c           I-II          I-t       
                                                                                     
nacncoop_gbar    3000.  [4]     1000.  [4]     1000.  [4]    1000.  [4]    1000.  [4] 
kht_gbar         58.0   [1]     58.0   [1]     500.0  [1]    150.0  [2]    500.0  [1] 
klt_gbar         80.0   [1]     20.0   [1]     0.0    [1]    14.0   [3]    0.0    [1] 
ka_gbar          0.0    [1]     0.0    [1]     0.0    [1]    0.0    [2]    125.0  [1] 
ihvcn_gbar       30.0   [1]     30.0   [1]     18.0   [1]    2.0    [2]    18.0   [1] 
leak_gbar        2.0    [1]     2.0    [1]     8.0    [1]    2.0    [2]    8.0    [1] 
leak_erev        -65    [1]     -65    [1]     -65    [1]    -65    [2]    -65    [1] 
na_type          nacncoop [4]   nacncoop  [4]  nacncoop [4]  nacncoop [3]  nacncoop [4] 
ih_type          ihvcn  [1]     ihvcn  [1]     ihvcn  [1]    ihvcn  [2]    ihvcn  [1] 
soma_Cap         26.0   [1]     26.0   [1]     25.0   [1]    25.0   [2]    25.0   [1] 
nacncoop_vshift  0.0    [1]     0.0    [1]     0.0    [1]    0.0    [1]    0.0    [1]
e_k              -84    [1]     -84    [1]     -84    [1]    -70    [3]    -84    [1] 
e_na             50.    [1]     50.    [1]     50.    [1]    55.    [3]    50.    [1] 
ih_eh            -43    [1]     -43    [1]     -43    [1]    -43    [2]    -43    [1] 

-----------------------------------------------------------------------------------------------------------------------------------

[1] Uses channels from Xie and Manis, 2013
    Age "adult", Temperature=34C
    Units are nS.

[2] Rothman and manis, 2003, model I-II
    Some low-voltage K current, based on observations of
    a single spike near threshold and regular firing for higher
    currents (Xie and Manis, 2017)

[3] These values for the I-II (dstellate) are from the original checkpoint test
    for cnmodel 12/2017. 

[4] nav11 channels were used in original Xie and Manis (2013) ms, 
    However, this version uses cooperative na channels for faster activation

""")

add_table_data('XM13nacncooop_channels_compartments', row_key='parameter', col_key='compartment', 
               species='mouse', model_type='II', data=u"""

!!!!!!!!!!!! USAGE OF THIS TABLE SHOULD BE CONSIDERED EXPERIMENTAL !!!!!!!!!!!!!!

This table describes the ion channel densities relative to somatic densities,
e.g., relative to REFERENCE densities in the table XM13_nacncoop_channels.
and voltage shifts, for different compartments of the specified neuron,
Conductances will be calculated from the Model derived from Xie and Manis 2013 for mouse

------------------------------------------------------------------------------------------------------------------------------------------------------------------
                   axon           unmyelinatedaxon     myelinatedaxon     initialsegment    hillock     soma        dendrite         primarydendrite    secondarydendrite
                                                                                                                                                                                                      
nacncoop_gbar      3.0 [1]        3.0 [1]              0.0 [1]            5.0 [1]           5.0 [1]     1.0 [1]     0.5 [1]          0.50 [1]           0.25 [1]       
kht_gbar           1.0 [1]        2.0 [1]              0.01 [1]           2.0 [1]           2.0 [1]     1.0 [1]     0.5 [1]          0.5 [1]            0.25 [1]       
klt_gbar           1.0 [1]        1.0 [1]              0.01 [1]           1.0 [1]           1.0 [1]     1.0 [1]     0.5 [1]          0.5 [1]            0.25 [1]       
ihvcn_gbar         0.0 [1]        0.0 [1]              0.0 [1]            0.5 [1]           0.0 [1]     1.0 [1]     0.5 [1]          0.5 [1]            0.5 [1]       
leak_gbar          1.0 [1]        0.25 [1]             0.25e-3 [1]        1.0 [1]           1.0 [1]     1.0 [1]     0.5 [1]          0.5 [1]            0.5 [1]       
leak_erev          -65. [1]       -65. [1]             -65. [1]           -65. [1]          -65. [1]    -65. [1]    -65. [1]         -65. [1]           -65. [1]      
nacncoop_vshift    0.0  [1]       0.0  [1]             0.0 [1]            0.0  [1]          0.0  [1]    0.0 [1]     0.0 [1]          0.0 [1]            0.0 [1]       
na_type            nacncoop       nacncoop             nacncoop           nacncoop          nacncoop    nacncoop    nacncoop            nacncoop              nacncoop
ih_type            ihvcn          ihvcn                ihvcn              ihvcn             ihvcn       ihvcn       ihvcn            ihvcn              ihvcn                            
--------------------------------------------------------------------------------------------------------------------------------------------------------------------

[1] Scaling is relative to soma scaling. Numbers are estimates based on general distribution from literature on cortical neurons.


""")

# ***** END OF XM13_Channels for nacncoop version of model

add_table_data('mGBC_channels', row_key='field', col_key='model_type', 
               species='mouse', data=u"""

This table describes the REFERENCE ion channel densities (and voltage shifts if necessary)
for different cell types based on the Xie and Manis 2013 models for mouse.

The REFERENCE values are applied to "point" models, and to the soma of
compartmental models.
The names of the mechanisms must match a channel mechanism (Neuron .mod files)
and the following _(gbar, vshift, etc) must match an attribute of that channel
that can be accessed.

-----------------------------------------------------------------------------------------------------------------------------------
               II             II-I           I-c           I-II          I-t       
                                                                                   
nav11_gbar     1600.  [1]     1600.  [1]     3000.  [1]    1600.  [2]    3000.  [1] 
kht_gbar       58.0   [1]     58.0   [1]     500.0  [1]    150.0  [2]    500.0  [1] 
klt_gbar       80.0   [1]     14.0   [1]     0.0    [1]    20.0   [2]    0.0    [1] 
ka_gbar        0.0    [1]     0.0    [1]     0.0    [1]    0.0    [2]    125.0    [1] 
ihvcn_gbar     30.0   [1]     30.0   [1]     18.0   [1]    2.0    [2]    18.0   [1] 
leak_gbar      2.0    [1]     2.0    [1]     8.0    [1]    2.0    [2]    8.0    [1] 
leak_erev      -65    [1]     -65    [1]     -65    [1]    -65    [2]    -65    [1] 
na_type        nav11  [1]     nav11  [1]     nav11  [1]    nav11  [1]    nav11  [1] 
ih_type        ihvcn  [1]     ihvcn  [1]     ihvcn  [1]    ihvcn  [2]    ihvcn  [1] 
soma_Cap       26.0   [1]     26.0   [1]     25.0   [1]    26.0   [2]    25.0   [1] 
nav11_vshift   4.3    [1]     4.3    [1]     4.3    [1]    4.3    [1]    4.3    [1]
e_k            -84    [1]     -84    [1]     -84    [1]    -84    [2]    -84    [1] 
e_na           50.    [1]     50.    [1]     50.    [1]    50.    [2]    50.    [1] 
ih_eh          -43    [1]     -43    [1]     -43    [1]    -43    [2]    -43    [1] 

-----------------------------------------------------------------------------------------------------------------------------------

[1] Uses channels from Rothman and Manis, 2003, except for Na channels
    Conductances are for Mouse bushy cells
    Xie and Manis, 2013
    Age "adult", Temperature=34C
    Units are nS.

[2] Rothman and Manis, 2003, model I-II
    Some low-voltage K current, based on observations of
    a single spike near threshold and regular firing for higher
    currents (Xie and Manis, 2017)


""")



add_table_data('mGBC_channels_compartments', row_key='parameter', col_key='compartment', 
               species='mouse', model_type='II', data=u"""

This table describes the ion channel densities relative to somatic densities,
e.g., relative to REFERENCE densities in the table XM13_channels.
and voltage shifts, for different compartments of the specified neuron,
Conductances will be calculated from the Model for Xie and Manis 2013 for mouse
(data table: XM13_channels).

------------------------------------------------------------------------------------------------------------------------------------------------------------------
               axon       unmyelinatedaxon     myelinatedaxon     initialsegment    hillock     soma        dendrite         primarydendrite    secondarydendrite
                                                                                                                                                                                              
nav11_gbar     3.0 [1]    3.0 [1]              0.0 [1]            3.0 [1]           2.0 [1]     1.0 [1]     0.25 [1]         0.25 [1]           0.25 [1]       
kht_gbar       1.0 [1]    2.0 [1]              0.01 [1]           2.0 [1]           2.0 [1]     1.0 [1]     0.5 [1]          0.5 [1]            0.25 [1]       
klt_gbar       1.0 [1]    1.0 [1]              0.01 [1]           1.0 [1]           1.0 [1]     1.0 [1]     0.5 [1]          0.5 [1]            0.25 [1]       
ihvcn_gbar     0.0 [1]    0.0 [1]              0.0 [1]            0.5 [1]           0.0 [1]     1.0 [1]     0.5 [1]          0.5 [1]            0.5 [1]       
leak_gbar      1.0 [1]    0.25 [1]             0.25e-3 [1]        1.0 [1]           1.0 [1]     1.0 [1]     0.5 [1]          0.5 [1]            0.5 [1]       
leak_erev      -65. [1]   -65. [1]             -65. [1]           -65. [1]          -65. [1]    -65. [1]    -65. [1]         -65. [1]           -65. [1]      
nav11_vshift   4.3  [1]   4.3  [1]             0.0 [1]            4.3  [1]          4.3  [1]    0.0 [1]     0.0 [1]          0.0 [1]            0.0 [1]       
na_type        nav11      nav11                nav11              nav11             nav11       nav11       nav11            nav11              nav11
ih_type        ihvcn      ihvcn                ihvcn              ihvcn             ihvcn       ihvcn       ihvcn            ihvcn              ihvcn                            
-------------------------------------------------------------------------------------------------------------------------------------------------------------------

[1] Scaling is relative to soma scaling. Numbers are estimates based on general distribution from literature on cortical neurons.


""")


add_table_data('POK_channels', row_key='field', col_key='model_type', 
               species='rat', data=u"""

This table describes the ion channel densities and voltage shifts for rat DCN pyramidal cells,
from Kanold and Manis, 2001

------------------------------------------------------------------------------------------------------------------------------------------
                          pyramidal   
                                   
soma_napyr_gbar           350.0  [1]      
soma_nap_gbar             0.
soma_cap_pcabar           0.     [3]
soma_kdpyr_gbar           80.0   [1]
soma_kcnq_gbar            0.     [3]
soma_kpksk_gbar           0.     [3]
soma_kir_gbar             0.     [3]
soma_kif_gbar             150.0  [1]
soma_kis_gbar             40.0   [1]
soma_ihpyr_gbar           2.8    [1]     
soma_leak_gbar            2.8    [1]
soma_leak_erev            -62.0  [1]
soma_e_na                 50.    [1]
soma_e_k                  -81.5  [1]
soma_e_h                  -43.0  [1]
soma_natype               napyr
soma_Cap                  22.0   [1]
------------------------------------------------------------------------------------------------------------------------------------------

[1] Kanold and Manis, 1999, 2001, 2005
    Age P11-14, Temperature=22C
    Units are nS.
    Default cap is 12 pF.
[2] Adjustable q10 added for fitting
    soma_ihpyr_adj_q10        1.0    [2]      (removed for testing)

[3] for implementing the additional channels from Li et al., and Leao et al. Default remains
    original model set to 0; also see Ceballo et al. 2016.

""")

add_table_data('CW_channels', row_key='field', col_key='model_type', 
               species='mouse', data=u"""

This table describes the ion channel densities and voltage shifts
for a mouse carthweel cell model.
Ad-hoc model, based on a Purkinje cell model (ref [1]).


-----------------------------------------------------------------------------------------------------------------------------------
                   cartwheel   
                           
soma_narsg_gbar    500.0  [1]     
soma_bkpkj_gbar    2.0
soma_kpkj_gbar     100.   [1]
soma_kpkj2_gbar    50.
soma_kpkjslow_gbar 150    [1]
soma_kpksk_gbar    25.0   [1]
soma_lkpkj_gbar    5.0    [1]     
soma_hpkj_gbar     5.0    [1]
soma_e_na          50.    [1]
soma_e_k           -80.0  [1]
soma_hpkj_eh       -43.0  [1]
soma_lkpkj_e       -65.0  [1]
soma_e_ca          50.
soma_na_type       narsg
soma_pcabar        0.00015 [1]
soma_Dia           18  

-----------------------------------------------------------------------------------------------------------------------------------

[1] Channels from Khaliq, Gouwens and Raman, J. Neurosci. 2003
    Conductance levels modified. 

""")

add_table_data('TV_channels', row_key='field', col_key='model_type', 
               species='mouse', data=u"""

This table describes the ion channel densities and voltage shifts
for a mouse tuberculoventral cell model.
Ad-hoc model, based on the t-stellate cell model, but adjusted
to match the data from Kuo and Trussell.

-----------------------------------------------------------------------------------------------------------------------------------
                     TVmouse   
                              
soma_nacncoop_gbar   5800.0   [2]      
soma_kht_gbar        400.0    [1]
soma_ihvcn_gbar      2.5      [2]
soma_ka_gbar         65.0     [1]
soma_leak_gbar       4.5      [1]
soma_leak_erev       -72.0    [1]
soma_e_na            50.      [1]
soma_e_k             -81.5    [1]
soma_ihvcn_eh        -43.0    [1]
soma_na_type         nacncoop [2]
soma_Cap             35       [1]

-----------------------------------------------------------------------------------------------------------------------------------

[1] Values obtained from brute force runs and comparision to 
    FI curve from Kuo, Lu and Trussell, J Neurophysiol. 2012 Aug 15; 
    108(4): 1186–1198.

[2] Cooperative sodium channel model, based on (see the mechanisms folder)
    concepts and implementation similar to Oz et al. J.Comp. Neurosci. 39: 63, 2015,
    and Huang et al., PloSOne 7:e37729, 2012.


""")

add_table_data('sgc_mouse_channels', row_key='field', col_key='model_type', 
               species='mouse', data=u"""

This table describes the ion channel densities (and voltage shifts if necessary)
for SGC cells, based on 

-----------------------------------------------------------------------------------------------------------------------------------
                    sgc-a         sgc-bm   
                                        
sgc_name            a             bm     
soma_na_gbar        350.  [2]     350.  [2]
soma_kht_gbar       58.0  [1]     58.0  [1]
soma_klt_gbar       80.0  [1]     80.0  [1]
soma_ihap_gbar      3.0   [3]     0.0   [1]
soma_ihap_eh        -41.0 [3]     -41.0 [3]
soma_ihbm_gbar      0.0   [3]     3.0   [3]
soma_ihbm_eh        -41.0 [3]     -41.0 [3]
soma_leak_gbar      2.0   [1]     2.0   [1]
soma_leak_erev      -65   [1]     -65   [1]
soma_na_type        jsrna [2]     jsrna [2]
soma_Cap            12.0  [1]     12.0  [1]
soma_e_k            -84   [1]     -84   [1]
soma_e_na           50.   [1]     50.   [1]

-----------------------------------------------------------------------------------------------------------------------------------

[1] Model is based on the mouse bushy cell model (XM13, above),
    but with a fast sodium channel from Rothman et al, 1993. and Ih currents
    from Liu et al. 2014
    
[2] Sodium channel from Rothman, Young and Manis, J Neurophysiol. 1993 Dec;70(6):2562-83. 

[3] Ih Currents from Liu, Manis, Davis, J Assoc Res Otolaryngol. 2014 Aug;15(4):585-99.
    doi: 10.1007/s10162-014-0446-z. Epub 2014 Feb 21.
    Age "P10" (cultured SGC cells), Original data temperature=22C.
    Units are nS.

""")


add_table_data('sgc_guineapig_channels', row_key='field', col_key='model_type', 
               species='guineapig', data=u"""

This table describes the ion channel densities (and voltage shifts if necessary)
for a model SGC cell, which is based on a bushy cell with a different Na channel.

-----------------------------------------------------------------------------------------------------------------------------------
                    sgc-a         sgc-bm   
                                        
sgc_name            a             bm     
soma_na_gbar        1000. [2]     1000. [2]
soma_kht_gbar       150.0 [1]     150.0 [1]
soma_klt_gbar       200.0 [1]     200.0 [1]
soma_ihap_gbar      3.0   [3]     0.0   [3]
soma_ihap_eh        -41.0 [3]     -41.0 [3]
soma_ihbm_gbar      0.0   [3]     3.0   [3]
soma_ihbm_eh        -41.0 [3]     -41.0 [3]
soma_leak_gbar      2.0   [1]     2.0   [1]
soma_leak_erev      -65   [1]     -65   [1]
soma_na_type        jsrna [2]     jsrna [2]
soma_Cap            12.0  [1]     12.0  [1]
soma_e_k            -84   [1]     -84   [1]
soma_e_na           50.   [1]     50.   [1]

-----------------------------------------------------------------------------------------------------------------------------------

[1] Model is based on the guinea pig bushy cell model (RM03, above),
    but with a fast sodium channel from Rothman et al, 1993. and Ih currents
    from Liu et al. 2014
    
[2] Sodium channel from Rothman, Young and Manis, J Neurophysiol. 1993 Dec;70(6):2562-83. 

[3] Ih Currents from Liu, Manis, Davis, J Assoc Res Otolaryngol. 2014 Aug;15(4):585-99.
    doi: 10.1007/s10162-014-0446-z. Epub 2014 Feb 21.
    Age "P10" (cultured SGC cells), Temperature=22C.
    Units are nS.

""")

add_table_data('MSO_principal_channels', row_key='field', col_key='model_type', 
               species='guineapig', data=u"""

This table describes the ion channel densities
for a putative MSO principal neuron based on the original Rothman Manis 2003 model for bushy cells.

-----------------------------------------------------------------------------------------------------------------------------------
                    MSO-principal   
                             
MSO_name            Principal       
na_gbar             1000. [1]
soma_kht_gbar       150.0 [1]
soma_klt_gbar       200.0 [1]
soma_ka_gbar        0.0   [1]
soma_ih_gbar        20.0  [1]
soma_leak_gbar      2.0   [1]
soma_leak_erev      -65   [1]
soma_na_type        nacn  [1]
soma_ih_type        ihvcn [1]
soma_Cap            12.0  [1]
soma_e_k            -84   [1]
soma_e_na           50.   [1]
soma_ih_eh          -43   [1]

-----------------------------------------------------------------------------------------------------------------------------------

[1] This MSO neuron model is basied on Rothman and Manis, 2003 bushy cell, type II
    Age "adult", Temperature=22C
    Units are nS.


""")

