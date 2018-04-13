"""
The cnmodel.data package contains information about ion channel densities, 
connectivity, synaptic properties, and population distributions. These values
are used by the Cell, Synapse, Population, and related classes to determine
all model construction parameters.

Values are stored in python strings that contain human-readable tables with
provenance documentation.
"""

from ._db import get, get_source, add_table_data
from . import connectivity
from . import synapses
from . import populations
from . import ionchannels

