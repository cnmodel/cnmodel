# -*- encoding: utf-8 -*-
from ._db import add_table_data

mouse_convergence = u"""

Convergence defines the average number of presynaptic cells of a particular
type (rows) that synapse onto a single postsynaptic cell of a particular
type (columns).
This connectivity matrix is currently incomplete.
Note: Bushy and pyramidal cells are known to have no (or very few)
collaterals within the CN, and so they are not listed as presynaptic cells in
this table. Octopus cells have collaterals (including in granule cell domains),
and should be added to this table when more data are available (Golding et al.,
J. Neurosci. 15: 3138, 1995)

----------------------------------------------------------------------------------------------
                  bushy       tstellate   dstellate   octopus     pyramidal    tuberculoventral
sgc               3.3±0.6 [2] 6.5±1.0 [2] 35±0 [3]    60±0 [2]    48±0 [5]     24±0 [5]
dstellate         7 [1]       20 [1]      3 [1]       0 [4]       15 [5]       15 [5]
tstellate         0 [6]       0 [6]       0 [6]       0 [6]       0 [6]        0 [6]
tuberculoventral  6           6           0           0 [4]       21 [5]       0 [7]
pyramidal         0           0           0           0           0            0    
----------------------------------------------------------------------------------------------

[1] Guesses based on Campagnola & Manis 2014

[2] Cao, X. & Oertel, D. (2010). Auditory nerve fibers excite targets through
    synapses that vary in convergence, strength, and short-term plasticity. 
    Journal of Neurophysiology, 104(5), 2308–20.
    Xie and Manis (unpublished): max EPSC = 3.4 ± 1.5 nA with ~0.3 nA steps
    (Cao and Oertel, 2010) = ~11 AN inputs. However neither we nor Cao and Oertel
    see that many clear steps in the responses, so use lower bound.
    
[3] Lower bound based on estimates from unpublished data Xie and Manis (2017)
    Assumptions: No discernable step sizes when increasing shock intensity 
    at ANFs in radiate multipolars (dstellate)
     Measured: 0.034 ± 15 nA sEPSC @ -70 mV
     Measured: Maximal current from AN stim = 1.2 ± 0.7 nA @ -70 mV
     Assuming that each AN provides 1 input, then N = ~35
     
[4] Octopus cells are devoid of inhibitory input (Golding et al., J. Neurosci., 1995)

[5] Convergence from Hancock and Voigt, Ann. Biomed. Eng. 27, 1999 and Zheng and Voigt,
    Ann. Biomed. Eng., 34, 2006.  Numbers are based on models for cat and gerbil,
    respectively. Adjusted to 1/2 to avoid overexciting TV cells in network model.

[6] tstellate cells have collaterals within the CN. It has been proposed that they
    provide auditory-driven input to the DCN (Oertel and Young, ), and also synapse
    within the VCN (Oertel, SFN abstract). These parameters may need to be adjusted
    once the convergence and strength is known.

[7] In the models of Hancock and Voigt (1999) and Zheng and Voigt (2006), the TV cells
    have no connections with each other. However, Kuo et al. (J. Neurophysiol., 2015)
    did see connections between pairs of TV cells in the mouse.

"""

add_table_data('convergence', row_key='pre_type', col_key='post_type', 
               species='mouse', data=mouse_convergence)



mouse_convergence_range = u"""

The convergence range table describes, for each type of connection from 
presynaptic (rows) to postsynaptic (columns), the variance in frequency of
presynaptic cells relative to the postsynaptic cell.

All values are expressed as the sigma for a lognormal distribution scaled to
the CF of the postsynaptic cell. 

----------------------------------------------------------------------------------------------
                  bushy       tstellate   dstellate   octopus     pyramidal    tuberculoventral
sgc               0.05 [1]    0.1 [1]     0.4 [1]     0.5 [5]     0.1 [1]      0.1 [1]
dstellate         0.208 [2]   0.347 [2]   0.5 [1]     0           0.2 [1]      0.2 [1]      
tstellate         0.1 [4]     0.1 [4]     0           0           0            0    
tuberculoventral  0.069 [3]   0.111 [3]   0           0           0.15 [1]     0    
pyramidal         0           0           0           0           0            0    
----------------------------------------------------------------------------------------------

[1] Guess based on axonal / dendritic morphology.

[2] Calculated from Campagnola & Manis 2014 fig. 7C
    Distribution widths are given in stdev(octaves), so we multiply by ln(2) to
    get the sigma for a lognormal distribution.
        DS->Bushy:     ln(2) * 0.3 = 0.208
        DS->TStellate: ln(2) * 0.5 = 0.347

[3] Calculated from Campagnola & Manis 2014 fig. 9C
    Distribution widths are given in stdev(octaves), so we multiply by ln(2) to
    get the sigma for a lognormal distribution.
        TV->Bushy:     ln(2) * 0.10 = 0.069
        TV->TStellate: ln(2) * 0.16 = 0.111

[4] Guess based on very limited information in Campagnola & Manis 2014 fig. 12

[5] Octopus cells get a wide range of ANF input (but weak on a per input basis)
    For example, see McGinley et al., 2012 or Spencer et al., 2012.


"""

add_table_data('convergence_range', row_key='pre_type', col_key='post_type', 
               species='mouse', data=mouse_convergence_range)


