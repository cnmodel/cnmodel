# -*- encoding: utf-8 -*-
from ._db import add_table_data

mouse_convergence = u"""

Convergence defines the average number of presynaptic cells of a particular
type (rows) that synapse onto a single postsynaptic cell of a particular
type (columns).
Note: Bushy and pyramidal are known to have no (or very few)
collaterals within the CN, and so they are not listed as presynaptic cells in
this table. Octopus cells have collaterals (including in granule cell domains),
and should be added to this table when data are available (Golding et al.,
J. Neurosci. 15: 3138, 1995)

-------------------------------------------------------------------------------
                  bushy       tstellate   dstellate   octopus     pyramidal    tuberculoventral
sgc               3.3±0.6 [2] 6.5±1.0 [2] 35±0 [3]    60±0 [2]    48±0 [5]     48±0 [5]
dstellate         7 [1]       20 [1]      3 [1]       0 [4]       15 [5]       15 [5]
tstellate         0 [6]       0 [6]       0 [6]       0 [6]       0 [6]        0 [6]
tuberculoventral  6           6           0           0 [4]       21 [5]       0 [7]
-------------------------------------------------------------------------------

[1] Guesses based on Campagnola & Manis 2013

[2] Cao, X. & Oertel, D. (2010). Auditory nerve fibers excite targets through
    synapses that vary in convergence, strength, and short-term plasticity. 
    Journal of Neurophysiology, 104(5), 2308–20.
    Xie and Manis (unpublished): max EPSC = 3.4 ± 1.5 nA with ~0.3 nA steps
    (Cao and Oertel, 2010) = ~11 AN inputs. However neither we nor Cao and Oertel
    see that many clear steps in the responses, so use lower bound.
    
[3] Lower bound based on estimates from unpublished data Xie and Manis (2017)
    Assumptions: No discernable step sizes in radiate multipolars (dstellate)
     Measured: 0.034 ± 15 pA sEPSC @ -70 mV ()
     Measured: Maximal current from AN stim = 1.2 ± 0.7 nA @ -70 mV
     Assuming that each AN provides 1 input, then N = ~35
     
[4] Octopus cells are devoid of inhibitory input (Golding et al., J. Neurosci., 1995)

[5] Convergence from Hancock and Voigt, Ann. Biomed. Eng. 27, 1999 and Zheng and Voigt,
    Ann. Biomed. Eng., 34, 2006.  Numbers are based on models for cat and gerbil,
    respectively.

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

