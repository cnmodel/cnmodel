"""
Test connection between two cell populations.

Usage: python test_populations.py <pre_celltype> <post_celltype>

This script: 

1. Creates two cell populations (pop1 and pop2)
2. Connects pop1 => pop2
3. Instantiates a single cell in pop2
4. Automatically generates presynaptic cells and synapses from pop1
5. Stimulates presynaptic cells and records postsynaptically

This is a high-level approach to generating networks in that the supporting
cells (those in pop1) are created automatically based on expected patterns
of connectivity in the cochlear nucleus. A lower-level approach is demonstrated
in test_synapses.py, in which the individual pre- and postsynaptic cells are
manually created and connected.
"""
import sys
import argparse

from cnmodel import populations
from cnmodel.protocols import PopulationTest
import pyqtgraph as pg


def testpopulation():
    parser = argparse.ArgumentParser(description="test simple synapses")
    parser.add_argument(
        "cells",
        type=str,
        nargs=2,
        choices=[
            "sgc",
            "bushy",
            "tstellate",
            "dstellate",
            "tuberculoventral",
            "pyramidal",
        ],
        help="Specify source and target cells.",
    )
    args = parser.parse_args()

    pop_types = {
        "sgc": populations.SGC,
        "bushy": populations.Bushy,
        "tstellate": populations.TStellate,
        "dstellate": populations.DStellate,
        "pyramidal": populations.Pyramidal,
        "tuberculoventral": populations.Tuberculoventral,
    }

    pops = []
    for cell_type in args.cells:
        if cell_type not in pop_types:
            print(
                '\nUnsupported cell type: "%s". Options are %s'
                % (cell_type, list(pop_types.keys()))
            )
            sys.exit(-1)
        pops.append(pop_types[cell_type]())

    pt = PopulationTest()
    pt.run(pops)
    pt.show()

    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()


if __name__ == "__main__":
    testpopulation()
