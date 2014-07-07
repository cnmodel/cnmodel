"""
Cell definitions for models.

This file includes a number of different cell definitions and default
conductances for point models. Most are models from the lab for neurons
of the cochlear nucleus.
Pyramidal cell from the DCN: pyr (Kanold and Manis, 1999, 2001, 2005)
Bushy cell from the VCN: bushy (Rothman and Manis, 2003abc) Type II
T-stellate cell from teh VCN: tstellate (Rothman and Manis, 2003abc) Type I

Paul B. Manis, Ph.D. 2009 (August - November 2009)
"""

from .bushy import *

def getcelltypes():
    """
    Return the types of cells that area generally avaialble.
    This is sorely outdated for the files here.
    """
    return({'pyramidal': (pyr, 'pyrprc.txt',
                'PYR-PRC.py output file for pyramidal'),
            'cartwheel': (cartwheel, 'cwprc.txt',
                'CW-PRC.py output file for Cartwheel'),
            'bushy': (bushy, 'bushyprc.txt',
                'BU-PRC.py output file for Bushy'),
            'tstellate': (tstellate_rothman, 'tstellateprc.txt',
                'ST-PRC.py output file for T-Stellate')})
