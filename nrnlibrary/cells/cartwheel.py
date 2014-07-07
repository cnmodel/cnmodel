from neuron import h
import neuron as nrn
from ..pynrnutilities import nstomho

from .cell import Cell

__all__ = ['Cartwheel']


class Cartwheel(Cell):
    """
    DCN cartwheel cell model.
    
    """
    def __init__(self, debug=False):
        soma = h.Section() # one compartment of about 29000 um2
        soma.nseg = 1
        soma.diam = 96
        soma.L = 96
        cm = 1
        v_potassium = -80       # potassium reversal potential
        v_sodium = 50           # sodium reversal potential

        seg = soma
        seg.insert('naRsg')
        seg.insert('kpkj')
        seg.insert('kpkj2')
        seg.insert('kpkjslow')
        seg.insert('bkpkj')
        seg.insert('kpksk')
        seg.insert('cadiff')
        seg.insert('cap')
        seg.insert('lkpkj')
        seg.insert('hpkj')
        seg.ena = 60
        seg.ek = -80
        s = soma()
        s.kpksk.gbar = 0.002
        if debug:
            print "<< cartwheel: Raman Purkinje cell model (modified) created >>"
            
        self.soma = soma
