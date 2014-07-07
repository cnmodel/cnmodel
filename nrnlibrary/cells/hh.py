from neuron import h
import neuron as nrn
from ..pynrnutilities import nstomho

from .cell import Cell

__all__ = ['HH']


class HH(Cell):
    """
    Standard Hodgkin-Huxley mechanisms from NEURON
    """
    def __init__(self, debug=False, message=None):
        soma = h.Section() # one compartment of about 29000 um2
        v_potassium = -80       # potassium reversal potential
        v_sodium = 50           # sodium reversal potential
        c_m = 1.0
        scalefactor = 1.0 # This determines the relative size of the cell
        rinsf = 1.0           # input resistance adjustment (also current...)
        totcap = 20.0 # scalefactor * 1.0 # cap in pF for cell
        effcap = totcap # sometimes we change capacitance - that's effcap
        somaarea = totcap * 1E-6 / c_m # pf -> uF, cm = 1uf/cm^2 nominal
        lstd = 1E4 * ((somaarea / 3.14159) ** 0.5) # convert from cm to um

        soma.nseg = 1
        soma.diam = lstd
        soma.L = lstd

        seg = soma
        seg.insert('hh')
        seg.insert('pas')
        if debug:
            if message is None:
                print "<< Standard HH model created >>"
            else:
                print message
        
        self.soma = soma
