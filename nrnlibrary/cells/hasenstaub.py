from neuron import h
import neuron as nrn
from ..util import nstomho

from .cell import Cell

__all__ = ['Hasenstaub']


class Hasenstaub(Cell):
    def __init__(self, debug=False, ttx=False, message=None, pump=False):
        super(Hasenstaub, self).__init__()
        v_potassium = -70       # potassium reversal potential
        v_sodium = 50           # sodium reversal potential
        v_chloride = -20
        cm = 1.0
        Ra = 150
        scalefactor = 1.0 # This determines the relative size of the cell
        rinsf = 1.0           # input resistance adjustment (also current...)
        totcap = 1e5 # in pF

        effcap = totcap # sometimes we change capacitance - that's effcap
        somaarea = totcap * 1E-6 / cm # pf -> uF, cm = 1uf/cm^2 nominal
        lstd = 1E4 * ((somaarea / 3.14159) ** 0.5) # convert from cm to um

        soma = h.Section() # one compartment of about 29000 um2
        soma.nseg = 1
        soma.diam = lstd
        soma.L = lstd
        print 'soma: (um len, um dia, area in cm2)', soma.L, soma.diam, somaarea

        soma.insert('hstb')
        gkb = nstomho(7000, somaarea) * scalefactor
        gnab = nstomho(1500, somaarea) * scalefactor
        gl = nstomho(100, somaarea) * scalefactor
        soma().hstb.gbar = gnab
        soma().hstb.gbar = gkb
        soma().hstb.gl = gl
        print 'gkb, gnab: ', gkb, gnab
        if pump:
            soma.insert('k_conc')

            ki0_k_ion = 140
            soma().ki = ki0_k_ion
            soma().ki0_k_conc = ki0_k_ion
            soma().beta_k_conc = 0.075

            soma.insert('na_conc')
            nai0_na_ion = 5
            soma().nai = nai0_na_ion
            soma().nai0_na_conc = nai0_na_ion
            soma().beta_na_conc = 0.075

            soma.insert('nakpump')
            soma().nakpump.inakmax = 8
            soma().nao = 145
            soma().ko = 5
            soma().nakpump.Nai_inf = 5
            soma().nakpump.Ki_inf = 140
            soma().nakpump.ATPi = 5
        soma.ek = v_potassium
        soma().v = -60.
        if ttx:
            gbar = 0.0
            soma().na.gbar = gbar
            soma.ena = 50

        if debug:
            if message is None:
                print "<< hasenstaub created >>"
            else:
                print message
                
        self.add_section(soma, 'soma')
