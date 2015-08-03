"""
Cell definitions for models.

This file includes a number of different cell definitions and default
conductances for point models. Most are models from the lab for neurons
of the cochlear nucleus.
Pyramidal cell from the DCN: pyr (Kanold and Manis, 1999, 2001, 2005)
Bushy cell from the VCN: bushy (Rothman and Manis, 2003abc) Type II
T-stellate cell from the VCN: tstellate (Rothman and Manis, 2003abc) Type I-c and I-t

Paul B. Manis, Ph.D. 2009 (August - November 2009)
"""

from .bushy import *
from .tstellate import *
from .dstellate import *
from .cartwheel import *
from .pyramidal import *
from .sgc import *
from .octopus import *
from .hh import *

from .cell import Cell

def cell_from_section(sec):
    return Cell.from_section(sec)

