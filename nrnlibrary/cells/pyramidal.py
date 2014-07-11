from neuron import h
import neuron as nrn
from ..util import nstomho

from .cell import Cell

__all__ = ['Pyramidal']


class Pyramidal(Cell):
    """
    DCN pyramidal cell
    Kanold and Manis, 1999, 2001, 2005
    """
    def __init__(self, scalefactor=1.0, debug=False):
        ndend = 1
        soma = h.Section()
        dend = h.Section()

        cm = 1.0
        scalefactor = 1.0 # This determines the relative size of the cell
        rinsf = 1.0           # input resistance adjustment (also current...)
        totcap = scalefactor * 12.0 # cap in pF for cell
        effcap = totcap # sometimes we change capacitance - that's effcap
        somaarea = totcap * 1E-6 / cm # pf -> uF, cm = 1uf/cm^2 nominal;
        lstd = 1E4 * ((somaarea / 3.14159) ** 0.5) # convert from cm to um

        gnab = nstomho(350, somaarea) * scalefactor
        gnap = 0.0
        gkb = nstomho(80, somaarea) * scalefactor # used to be 20?
        gkfb = nstomho(150, somaarea) * scalefactor
        gksb = nstomho(40, somaarea) * scalefactor
        glb = nstomho(2.8, somaarea) * scalefactor
        ghb = nstomho(3, somaarea) * scalefactor
        gkpksk = nstomho(0, somaarea) * scalefactor
        gkir = nstomho(0, somaarea) * scalefactor # incude KIR here, but set to 0
        stdrin = 300 * rinsf / scalefactor
        stdrmp = -60
        if debug:
            print ("in Mhos: gna: %9.3g  gkb: %9.3g  gkab: %9.3g  gksb: %9.3g"
                % (gnab, gkb, gkfb, gksb))
            print ("glb: %9.3g  ghb: %9.3g gkpksk: %9.3g" %
                (glb, ghb, gkpksk))

    # set up soma like a pyramidal cell
        soma.nseg = 1
        soma.diam = lstd
        soma.L = lstd # these are expressed in microns...
        soma.insert('pyr')
        soma.insert('kpksk')
        soma.insert('cadiff') # diffusion
        soma.insert('cap') # p-type calcium current
        #soma.insert('nacum') # sodium accumulation (yes!)
        soma.insert('nakpump') # and of course a pump to handle it.
        soma.insert('k_conc')
        soma.insert('na_conc')
        # soma.insert('kna')

        seg = soma()
        seg.kpksk.gbar = gkpksk
        seg.cap.pcabar = 0.00002
        seg.pyr.gbar = gnab
        seg.pyr.gbar = gnap
        seg.pyr.gbar = gkb
        seg.pyr.gbar = gkfb
        seg.pyr.gbar = gksb
        seg.pyr.gl = glb
        seg.pyr.gbar = ghb
    # seg.pyr.gbar = gkir
        seg.pyr.kif_ivh = -89.6
        if debug:
            print " "
            print "<< PYR: POK Pyramidal Cell created >>"
            print " "
            
        self.soma = soma
