from neuron import h
import neuron as nrn
from ..pynrnutilities import nstomho

from .cell import Cell

__all__ = ['SGC']


class SGC(Cell):
    """
    Spiral ganglion cell model
    """
    def __init__(self, debug=False, ttx=False, message=None, nach='jsrna',
            species='mouse', chlist=None):
        if chlist is None:
            chlist = ['ih', 'klt', 'kht', 'na']
        v_potassium = -84       # potassium reversal potential
        v_sodium = 55           # sodium reversal potential
        v_chloride = -20
        v_eh = -41.3 # from R. Davis papers
        cm = 1.0
        scalefactor = 1.0 # This determines the relative size of the cell
        rinsf = 1.0           # input resistance adjustment (also current...)
        if species == 'guineapig':
            totcap = scalefactor * 12.0 # cap in pF for cell
        if species == 'mouse':
            totcap = scalefactor * 12.0 # cap in pF for cell from JS-S
        effcap = totcap # sometimes we change capacitance - that's effcap
        somaarea = totcap * 1E-6 / cm # pf -> uF, cm = 1uf/cm^2 nominal;
        lstd = 1E4 * ((somaarea / 3.14159) ** 0.5) # convert from cm to um

        soma = h.Section() # one compartment of about 29000 um2
        soma.nseg = 1
        soma.diam = lstd
        soma.L = lstd

        if 'klt' in chlist:
            soma.insert('klt')
        if 'kht' in chlist:
            soma.insert('kht')
        if 'na' in chlist:
            if nach == 'jsrna':
                soma.insert('jsrna')
            elif nach == 'nav11':
                soma.insert('nav11')
            else:
                soma.insert('na')
        #soma.insert('ka')
        if 'ih' in chlist:
            soma.insert('ihsgc')
        #soma.insert('hcno')
        soma.insert('leak')
        if 'klt' in chlist or 'klt' in chlist:
            soma.ek = v_potassium

        s = soma()
        gnabar = nstomho(1000.0, somaarea) * scalefactor
        if 'na' in chlist:
            if ttx:
                gnabar = 0.0
            if nach == 'jsrna':
                s.jsrna.gna = gnabar
                soma.ena = 50
            elif nach == 'nav11':
                s.nav11.gnatbar = gnabar * 0.5
                soma.ena = 50
                s.nav11.vsna = 4.3
                print "sgc using inva11"
            else:
                s.na.gbar = gnabar
                soma.ena = 50
        if species == 'mouse':
            # use conductance levels from Cao et al.,  J. Neurophys., 2007.
            #if debug:
            print 'Mouse sgc cell'
            #s.nav11.gbar = gnabar*0.3
            if 'kht' in chlist:
                s.kht.gbar = nstomho(58.0, somaarea) * scalefactor
            if 'klt' in chlist:
                s.klt.gbar = nstomho(80.0, somaarea) * scalefactor
                # nstomho(200.0, somaarea) * scalefactor
            #s.ka.gkabar = nstomho(0.0, somaarea) * scalefactor
            if 'ih' in chlist:
                s.ihsgc.gbar = nstomho(10.0, somaarea) * scalefactor
                s.ihsgc.eh = v_eh
            #s.hcno.gbar = 0.0
            s.leak.gbar = nstomho(0.5, somaarea) * scalefactor
            vm0 = -63.6
        if species == 'guineapig-sgc-II':
            # guinea pig data from Rothman and Manis, 2003, type II
            if debug:
                print 'Guinea pig sgc cell type II'
            if 'kht' in chlist:
                s.kht.gbar = nstomho(150.0, somaarea) * scalefactor
            if 'klt' in chlist:
                s.klt.gbar = nstomho(200.0, somaarea) * scalefactor
                # nstomho(200.0, somaarea) * scalefactor
            #s.ka.gbar = nstomho(0.0, somaarea) * scalefactor
            if 'ih' in chlist:
                s.ihsgc.gbar = nstomho(1.0, somaarea) * scalefactor
            #s.hcno.gbar = 0.0
            s.leak.gbar = nstomho(2.0, somaarea) * scalefactor
            vm0 = -63.6

        if debug:
            if message is None:
                print "<< sgc: SGC cell model created >>"
            else:
                print message
                
        self.soma = soma
