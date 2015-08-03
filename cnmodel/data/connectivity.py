# -*- encoding: utf-8 -*-
from ._db import add_table_data

mouse_convergence = u"""

Convergence defines the average number of presynaptic cells of a particular
type (rows) that synapse onto a single postsynaptic cell of a particular
type (columns). 

-------------------------------------------------------------------------------
                  bushy       tstellate   dstellate  pyramidal    octopus   tuberculoventral
sgc               3.3±0.6 [2] 6.5±1.0 [2] 10                      60 [2]
dstellate         7 [1]       20 [1]      3 [1]
tuberculoventral  6           6           0
-------------------------------------------------------------------------------

[1] Wild guesses based on Campagnola & Manis 2013
[2] Cao, X. & Oertel, D. (2010). Auditory nerve fibers excite targets through
    synapses that vary in convergence, strength, and short-term plasticity. 
    Journal of Neurophysiology, 104(5), 2308–20.
"""

add_table_data('convergence', row_key='pre_type', col_key='post_type', 
               species='mouse', data=mouse_convergence)

