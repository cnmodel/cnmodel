from ._db import add_table_data

mouse_convergence = """

Convergence defines the average number of presynaptic cells of a particular
type (rows) that synapse onto a single postsynaptic cell of a particular
type (columns). 

-------------------------------------------------------------------------------
                  bushy  tstellate   dstellate  pyramidal    octopus   tuberculoventral
sgc               2      6           10
dstellate         7 [1]  20 [1]      3 [1]
tuberculoventral  6      6           0
-------------------------------------------------------------------------------

[1] Wild guesses based on Campagnola & Manis 2013
"""

add_table_data('convergence', row_key='pre_type', col_key='post_type', 
               species='mouse', data=mouse_convergence)

