from neuron import h
import neuron as nrn
from ..pynrnutilities import nstomho

from .cell import Cell

__all__ = ['DStellate', 'DStellateIF', 'DStellateEager']


class DStellate(Cell):
    """ 
    VCN D-stellate model.
    as a type 2-1 from RM03 
    """
    def __init__(self, debug=False, ttx=False, message=None):
        soma = h.Section() # one compartment of about 29000 um2
        v_potassium = -80       # potassium reversal potential
        v_sodium = 50           # sodium reversal potential

        cm = 1.0
        scalefactor = 1.0 # This determines the relative size of the cell
        rinsf = 1.0           # input resistance adjustment (also current...)
        totcap = scalefactor * 25.0
        #  Larger cap value of 25 best fit for mouse D-stellate cell
        effcap = totcap # sometimes we change capacitance - that's effcap
        somaarea = totcap * 1E-6 / cm # pf -> uF, cm = 1uf/cm^2 nominal
        lstd = 1E4 * ((somaarea / 3.14159) ** 0.5) # convert from cm to um

        soma.nseg = 1
        soma.diam = lstd
        soma.L = lstd

        seg = soma
        seg.insert('kht')
        seg.insert('klt')
        seg.insert('nav11')
        seg.insert('ka')
        seg.insert('ihvcn')
        #seg.insert('iH_std')
        seg.insert('leak')
        seg.ena = 20
        seg.ek = -77
        #seg.eh = -43 # Rodrigues and Oertel, 2006
        s = soma()
        if not ttx:
            s.nav11.gbar = nstomho(1000.0, somaarea) * scalefactor
        else:
            s.nav11.gbar = 0.0
        s.nav11.vsna = 8 # voltage shift
        s.kht.gbar = nstomho(250.0, somaarea) * scalefactor
        s.klt.gbar = nstomho(35.0, somaarea) * scalefactor
        s.ka.gbar = nstomho(0.0, somaarea) * scalefactor
        s.ihvcn.gbar = nstomho(3.5, somaarea) * scalefactor
        s.ihvcn.vshift = 0
        s.leak.gbar = nstomho(2, somaarea) * scalefactor
        vm0 = -64.1
        if debug:
            if message is None:
                print ("<< D-stellate: JSR Stellate",
                " Type I-II cell model created >>")
            else:
                print message
                
        self.soma = soma


class DStellateIF(Cell):
    """
    Integrate and fire version of the VCN D-stellate cell.
    """
    def __init__(self, debug=False, ttx=False, message=None):
        soma = h.Section() # one compartment of about 29000 um2
        v_potassium = -80       # potassium reversal potential
        v_sodium = 50           # sodium reversal potential

        cm = 1.0
        scalefactor = 1.0 # This determines the relative size of the cell
        rinsf = 1.0           # input resistance adjustment (also current...)
        totcap = scalefactor * 25.0
        #  Larger value of 25 best fit for mouse D-stellate cell
        effcap = totcap # sometimes we change capacitance - that's effcap
        somaarea = totcap * 1E-6 / cm # pf -> uF, cm = 1uf/cm^2 nominal;
        lstd = 1E4 * ((somaarea / 3.14159) ** 0.5) # convert from cm to um

        soma.nseg = 1
        soma.diam = lstd
        soma.L = lstd

        seg = soma
        seg.insert('kht')
        seg.insert('klt')
        seg.insert('nav11')
        seg.insert('ka')
        #seg.insert('ihvcn')
        seg.insert('iH_std')
        seg.insert('leak')
        seg.ena = 20
        seg.ek = -77
        seg.eh = -40 # Rodrigues and Oertel, 2006
        s = soma()
        if not ttx:
            s.nav11.gbar = nstomho(3500.0, somaarea) * scalefactor
        else:
            s.na.gbar = 0.0
        s.nav11.vsna = 8 # voltage shift
        s.kht.gbar = nstomho(150.0, somaarea) * scalefactor
        s.klt.gbar = nstomho(10.0, somaarea) * scalefactor
        s.ka.gbar = nstomho(0.0, somaarea) * scalefactor
        s.iH_std.gbar = nstomho(120.0, somaarea) * scalefactor
        s.iH_std.vshift = 16
        s.leak.gbar = nstomho(24, somaarea) * scalefactor
        vm0 = -64.1
        if debug:
            if message is None:
                print("<< D-stellate: JSR Stellate",
                " Type I-II cell model created >>")
            else:
                print message
                
        self.soma = soma


class DStellateEager(Cell):
    """
    This is a model of the VCN D-Stellate cells as proposed by
    Eager, M.A., Grayden, D.B., Burkitt, A.N., and Meffin, H.,
    "A neural circuit model of the ventral cochlear nucleus",
    Internet:
    http://citeseerx.ist.pus.edu/viewdoc/download?doi=10.1.79.9620.pdf&rep
    =rep&type=pdf
    also cited as:
    Proceedings of theh 10th AUstralian International Conference on
    Speech Science and Technology,
    pp. 539-544, 2004.
    it is based on the Rothman and Manis (2003c) model,
    with small modifications.
    Their model includes dendrites...
    """
    def __init__(self, debug=False, ttx=False, message=None):
        soma = h.Section() # one compartment of about 29000 um2
        v_potassium = -80       # potassium reversal potential
        v_sodium = 50           # sodium reversal potential

        cm = 0.9
        Ra = 150
        lstd = 25.0

        soma.nseg = 1
        soma.diam = lstd
        soma.L = lstd

        axon = h.Section()
        axon.nseg = 2
        axon.diam = 3.0
        axon.L = 70.0

        axon.connect(soma, 0, 0)

        nDendrite = 2
        dend = nDendrite * [None]

        for i in range(len(dend)):
            dend[i] = h.Section()
            dend[i].nseg = 5
            dend[i].L = 1100.0
            dend[i].diam = 3.5
            dend[i].Ra = 1500.0
            dend[i].cm = 0.9
            dend[i].connect(soma, 0, 1)

        activecompartments = [soma, axon, dend[0], dend[1]]
        #print activecompartments
        for seg in activecompartments:
            s = seg()
            seg.insert('kht')
            seg.insert('klt')
            seg.insert('nav11')
            seg.insert('ka')
            seg.insert('ihvcn')
            #seg.insert('iH_std')
            seg.insert('leak')
            seg.ena = 20
            seg.ek = -77
            # seg.eh = -40 # Rodrigues and Oertel, 2006
            seg.Ra = 150.0
            seg.cm = 0.9
            if not ttx:
                gna = 0.5
            else:
                gna = 0.0
            if seg not in dend:
                s.nav11.gbar = gna # nstomho(gna, somaarea) * scalefactor
                s.nav11.vsna = 8 # voltage shift
                s.kht.gbar = 0.01 # nstomho(500.0, somaarea) * scalefactor
                s.klt.gbar = 0.005 # nstomho(125.0, somaarea) * scalefactor
                s.ka.gbar = 0.0 # nstomho(0.0, somaarea) * scalefactor
                s.ihvcn.gbar = 0.0001 # nstomho(5.0, somaarea) * scalefactor
                s.ihvcn.vshift = 0 # 16
                s.leak.gbar = 0.00025 # nstomho(12.5, somaarea) * scalefactor
            else:
                s.nav11.gbar = 0 # nstomho(gna, somaarea) * scalefactor
                s.nav11.vsna = 0 # voltage shift
                s.kht.gbar = 0.00 # nstomho(500.0, somaarea) * scalefactor
                s.klt.gbar = 0.001 # nstomho(125.0, somaarea) * scalefactor
                s.ka.gbar = 0.0 # nstomho(0.0, somaarea) * scalefactor
                s.ihvcn.gbar = 0.0001 # nstomho(5.0, somaarea) * scalefactor
                s.ihvcn.vshift = 0 # 16
                s.leak.gbar = 0.00025 # nstomho(12.5, somaarea
            print 'nav: ', s.nav11.gbar
            print 'khe: ', s.kht.gbar
            print 'klt: ', s.klt.gbar
            print 'ih:  ', s.ihvcn.gbar
            print 'gleak:', s.leak.g
        vm0 = -64.1
        if debug:
            if message is None:
                print ("<< D-stellate: Eager et al. ",
                " Type I-II (D-stellate) cell model created >>")
            else:
                print message
                
        self.soma = soma
        self.dend = dend
