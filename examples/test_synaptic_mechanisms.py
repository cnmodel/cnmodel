import argparse
from pathlib import Path
import numpy as np
import pickle
import platform
import datetime

print(platform.node())

from neuron import h
from cnmodel.util import reset

"""
Test routine to look carefully at NMDA_KAMPA nmda_states during a run.
"""


def test_max_open_probability(fn=None, savetype=False):
    print("Running sim for NMDARs; ")
    reset(
        raiseError=False
    )  # reset() fails as unable to remove all neuron objects, unless we ignore the error
    postsec = h.Section(0.5)

    # Create AMPA and NMDA mechanisms
    # AMPA uses mode=0; no rectification
    apsd = h.AMPATRUSSELL(postsec(0.5), sec=postsec)
    # For NMDA we will hold the cell at +40 mV
    npsd = h.NMDA_Kampa(postsec(0.5), sec=postsec)
    # npsd.loc(sec(0.5))
    # And a presynaptic terminal to provide XMTR input
    presec = h.Section(0.5)
    term = h.MultiSiteSynapse(presec(0.5), sec=presec)
    term.nZones = 2
    term.multisite = 0
    azone = 0
    nzone = 0
    h.setpointer(term._ref_XMTR[azone], "XMTR", apsd)
    h.setpointer(term._ref_XMTR[nzone], "XMTR", npsd)
    print("h.dt: ", h.dt)
    print("terminal nzones: ", term.nZones)
    h.celsius = 34.0
    h.finitialize()
    # h.dt = 0.025
    print("a")
    tstop = 10.0
    npts = int(tstop / h.dt) + 1
    nmda_states = [
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
        term._ref_XMTR[azone],
    ]

    nmda_statenames = [
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

    ampa_states = [
        apsd._ref_C0,
        apsd._ref_C1,
        apsd._ref_C2,
        apsd._ref_D,
        apsd._ref_O1,
        apsd._ref_O2,
    ]

    ampa_statenames = [
        "C0",
        "C1",
        "C2",
        "D",
        "O1",
        "O2",
    ]
    # NMDAS = dict.fromkeys(nmda_statenames, np.zeros(npts))
    # NMDAS['t'] = np.zeros(npts)
    xmtrs = {"nmda": h.Vector(), "ampa": h.Vector()}
    print(dir(term))
    xmtrs["nmda"].record(h.ref(term.XMTR[nzone]))
    xmtrs["ampa"].record(h.ref(term.XMTR[azone]))
    NSM = {}  # {k: nmda_states[i] for i, k in enumerate(nmda_statenames)}
    for i, k in enumerate(nmda_statenames):
        NSM[k] = h.Vector()  # Membrane potential vector
        if k not in ["Glut"]:
            NSM[k].record(nmda_states[i])
        else:
            NSM[k].record(nmda_states[i])
    NSM["t"] = h.Vector()
    NSM["t"].record(h._ref_t)

    ASM = {}  # {k: nmda_states[i] for i, k in enumerate(nmda_statenames)}
    for i, k in enumerate(ampa_statenames):
        ASM[k] = h.Vector()  # Membrane potential vector
        if k not in ["Glut"]:
            ASM[k].record(ampa_states[i])
        else:
            ASM[k].record(ampa_states[i])
    ASM["t"] = h.Vector()
    ASM["t"].record(h._ref_t)

    h.finitialize()
    h.t = 0.0
    h.dt = 0.025
    while h.t < tstop:
        # force very high transmitter concentration after a delay
        if h.t < 1.0:
            term.XMTR[azone] = 1000.0
            term.XMTR[nzone] = 1000.0
        else:
            term.XMTR[azone] = 10000.0
            term.XMTR[nzone] = 10000.0
        postsec.v = 40.0
        h.fadvance()
        # print(term._ref_XMTR[azone], term._ref_XMTR[nzone])

    for s in NSM.keys():
        NSM[s] = np.array(NSM[s])  # no longer vector objects
    for s in ASM.keys():
        ASM[s] = np.array(ASM[s])  # no longer vector objects
    print(np.max(NSM["Open"]), npsd.MaxOpen)
    print(max(ASM["O1"] + ASM["O2"]), apsd.MaxOpen)
    if not savetype:
        sdata = NSM
    else:
        sdata = ASM
    sdata = {
        "ampa": ASM,
        "nmda": NSM,
        "machine": platform.node(),
        "datetime": datetime.datetime.now(),
    }
    with open(fn, "wb") as fh:
        pickle.dump(sdata, fh)

    for k in xmtrs.keys():
        xmtrs[k] = np.array(xmtrs[k])
    for i in range(len(xmtrs["ampa"])):
        print(f"{NSM['t'][i]:.3f}, {xmtrs['ampa'][i]:.6f}, {xmtrs['nmda'][i]:.6f}")
    # assert np.allclose(max(ASM['O1']+ASM['O2']), apsd.MaxOpen)
    # assert np.allclose(max(NSM['Open']), npsd.MaxOpen)

    return


def plot(fn=None, which="nmda"):
    with open(fn, "rb") as fh:
        sdata = pickle.load(fh)
    P = PH.regular_grid(
        4,
        3,
        margins={
            "bottommargin": 0.1,
            "leftmargin": 0.07,
            "rightmargin": 0.05,
            "topmargin": 0.1,
        },
    )
    axs = P.axarr.ravel()
    print(sdata.keys())
    t = np.array(sdata[which]["t"])
    for i, s in enumerate(list(sdata[which].keys())):
        if s == "t":
            continue
        axs[i].plot(t, np.array(sdata[which][s]))
        axs[i].set_title(s, fontsize=9)
    P.figure_handle.suptitle(
        f"{which:s}   {sdata['machine']:s}  {sdata['datetime'].strftime('%H:%M:%S %m-%d-%Y'):s}",
        fontsize=10,
    )
    mpl.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="test nmda_states in NMDA receptor model"
    )
    parser.add_argument(
        "-p", "--plot", action="store_true", dest="plot", help="enable plot",
    )
    parser.add_argument(
        "-a",
        "--ampa",
        type=str,
        default="nmda",
        dest="ampa",
        choices=["ampa", "nmda"],
        help="Set save to AMPA instead of NMDA",
    )
    parser.add_argument(
        "-r",
        "--retrieve",
        action="store_true",
        dest="retrieve",
        help="retrieve plot without running simulation",
    )
    parser.add_argument(
        "-f", "--filename", default="test_syn_mech.p", help="set filename"
    )
    args = parser.parse_args()
    if not args.retrieve:
        test_max_open_probability(args.filename, savetype=args.ampa)
    if args.plot or args.retrieve:
        import pylibrary.plotting.plothelpers as PH
        import matplotlib.pyplot as mpl

        plot(args.filename, args.ampa)
