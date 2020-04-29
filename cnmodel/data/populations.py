# -*- encoding: utf-8 -*-
from ._db import add_table_data

add_table_data('populations', row_key='field', col_key='cell_type', 
               species='mouse', data=u"""

---------------------------------------------------------------------------------------------------------------------------
             sgc        bushy      tstellate    dstellate     octopus     pyramidal   tuberculoventral   pyramidal_ceballos  
                                                                                                                    
n_cells      10000 [1]  6500 [2]   6500 [2]     650 [3]       5000        3000        5000               3000       
cf_min       500        500        500          500           500         500         500                500        
cf_max       64000      64000      64000        64000         64000       64000       64000              64000      
--------------------------------------------------------------------------------------------------------------------------

[1] ?

[2] Rough estimate from allen brain atlas data:
    Volume of VCN is 0.377 mm^3, by counting voxels with 'VCO' (101) label in Common Coordinate Framework atlas.
        753370 voxels * 0.5 * 10e-6**3 m^3/vox = 0.377 mm^3
    Counted Slc17a7 (pan-excitatory) cell bodies in a 500x500 um chunk of VCN
        http://mouse.brain-map.org/experiment/siv?id=69014470&imageId=68856767&initImage=ish&coordSystem=pixel&x=7616.5&y=4144.5&z=1
        266 cells in 500x500 um = 34707 cells / mm^2
        34707**3/2 * 0.377 mm^3 = 13084 cells total
        Assume half are bushy, half are T-stellate
        
[3] Rough estimate from allen brain atlas data:
    Similar to [2], using Gad1 inhibitory marker
    http://mouse.brain-map.org/experiment/siv?id=75492764&imageId=75405134&initImage=ish&coordSystem=pixel&x=5320.5&y=3232.5&z=1
    36 cells in 500x500 um = 144e6 / m^2  ~= 1728 / mm^2
    = 651 cells total  (VCN, unilateral)

""")

add_table_data('populations', row_key='field', col_key='cell_type', 
               species='rat', data=u"""

---------------------------------------------------------------------------------------------------------------------------
             sgc        bushy      tstellate    dstellate     octopus     pyramidal   tuberculoventral   pyramidal_ceballos  
                                                                                                                    
n_cells      15800 [1]  6500 [2]   6500 [2]     650 [3]       5000        3000        5000               3000       
cf_min       500        500        500          500           500         500         500                500        
cf_max       64000      64000      64000        64000         64000       64000       64000              64000      
--------------------------------------------------------------------------------------------------------------------------

[1] Kiethly and Feldman J. Comp Neuol., 1979 (1-2 month old rat)

[2] All others are same as mice, until we have a better set of measurements
    Rough estimate from allen brain atlas data:
    Volume of VCN is 0.377 mm^3, by counting voxels with 'VCO' (101) label in Common Coordinate Framework atlas.
        753370 voxels * 0.5 * 10e-6**3 m^3/vox = 0.377 mm^3
    Counted Slc17a7 (pan-excitatory) cell bodies in a 500x500 um chunk of VCN
        http://mouse.brain-map.org/experiment/siv?id=69014470&imageId=68856767&initImage=ish&coordSystem=pixel&x=7616.5&y=4144.5&z=1
        266 cells in 500x500 um = 34707 cells / mm^2
        34707**3/2 * 0.377 mm^3 = 13084 cells total
        Assume half are bushy, half are T-stellate
        
[3] Rough estimate from allen brain atlas data:
    Similar to [2], using Gad1 inhibitory marker
    http://mouse.brain-map.org/experiment/siv?id=75492764&imageId=75405134&initImage=ish&coordSystem=pixel&x=5320.5&y=3232.5&z=1
    36 cells in 500x500 um = 144e6 / m^2  ~= 1728 / mm^2
    = 651 cells total  (VCN, unilateral)

""")
