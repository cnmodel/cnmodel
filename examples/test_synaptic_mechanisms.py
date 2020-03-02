import argparse
from pathlib import Path
import numpy as np
import pickle

from neuron import h
from cnmodel.util import reset

"""
Test routine to look carefully at NMDA_KAMPA states during a run.
"""


def test_max_open_probability(fn=None):
    print("Running sim for NMDARs")
    reset(
        raiseError=False
    )  # reset() fails as unable to remove all neuron objects, unless we ignore the error
    sec = h.Section()

    # Create AMPA and NMDA mechanisms
    # AMPA uses mode=0; no rectification
    apsd = h.AMPATRUSSELL(0.5, sec=sec)
    # For NMDA we will hold the cell at +40 mV
    npsd = h.NMDA_Kampa(0.5, sec=sec)

    # And a presynaptic terminal to provide XMTR input
    term = h.MultiSiteSynapse(0.5, sec=sec)
    term.nZones = 1
    h.setpointer(term._ref_XMTR[0], "XMTR", apsd)
    h.setpointer(term._ref_XMTR[0], "XMTR", npsd)
    print("h.dt: ", h.dt)
    h.celsius = 34.0
    h.finitialize()
    # h.dt = 0.025
    tstop = 50.0
    npts = int(tstop / h.dt) + 1
    states = [
        npsd._ref_U,
        npsd._ref_Cl,
        npsd._ref_D1,
        npsd._ref_D2,
        npsd._ref_Open,
        npsd._ref_MaxOpen,
        npsd._ref_UMg,
        npsd._ref_ClMg,
        npsd._ref_D1Mg,
        npsd._ref_D2Mg,
        npsd._ref_OMg,
        npsd._ref_XMTR,
    ]

    statenames = [
        "U",
        "Cl",
        "D1",
        "D2",
        "Open",
        "MaxOpen",
        "UMg",
        "ClMg",
        "D1Mg",
        "D2Mg",
        "OMg",
        "Glut",
    ]

    # NMDAS = dict.fromkeys(statenames, np.zeros(npts))
    # NMDAS['t'] = np.zeros(npts)

    NSM = {}  # {k: states[i] for i, k in enumerate(statenames)}
    for i, k in enumerate(statenames):
        NSM[k] = h.Vector()  # Membrane potential vector
        if k not in ['Glut']:
            NSM[k].record(states[i])
        else:
            NSM[k].record(states[i])
    NSM["t"] = h.Vector()
    NSM["t"].record(h._ref_t)

    h.finitialize()
    h.t = 0.0
    h.dt = 0.025
    while h.t < tstop:
        # force very high transmitter concentration for every timestep
        if h.t < 5.0:
            term.XMTR[0] = 0.
        else:
            term.XMTR[0] = 1.0
        # term.XMTR[0] = 10000.0
        sec.v = 40.0
        h.fadvance()

    for s in NSM.keys():
        NSM[s] = np.array(NSM[s])  # no longer vector objects
    # assert np.allclose(max(op1), apsd.MaxOpen)
    # assert np.allclose(max(op2), npsd.MaxOpen)
    with open(fn, "wb") as fh:
        pickle.dump(NSM, fh)
    return


def plot(fn=None):
    with open(fn, "rb") as fh:
        NMDAS = pickle.load(fh)
    P = PH.regular_grid(4, 3)
    axs = P.axarr.ravel()
    t = np.array(NMDAS["t"])
    for i, s in enumerate(list(NMDAS.keys())):
        if s == "t":
            continue
        axs[i].plot(t, np.array(NMDAS[s]))
        axs[i].set_title(s, fontsize=9)
    mpl.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="test states in NMDA receptor model")
    parser.add_argument(
        "-p", "--plot", action="store_true", dest="plot", help="enable plot",
    )
    parser.add_argument(
        "-r",
        "--retrieve",
        action="store_true",
        dest="retrieve",
        help="retrieve plot without running simulation",
    )
    parser.add_argument(
        "-f",
        "--filename",
        default="test_syn_mech.p",
        help="set filename"
    )
    args = parser.parse_args()
    if not args.retrieve:
        test_max_open_probability(args.filename)
    if args.plot or args.retrieve:
        import pylibrary.plotting.plothelpers as PH
        import matplotlib.pyplot as mpl
        plot(args.filename)
