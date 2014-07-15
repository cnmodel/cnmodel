from neuron import h
import neuron as nrn
from ..pynrnutilities import nstomho
import numpy as np
import scipy.optimize

from .cell import Cell

__all__ = ['Bushy']


class Bushy(Cell):
    """
    VCN bushy cell model.
    Rothman and Manis, 2003abc (Type II)    
    """

    def __init__(self, nach='nacn', ttx=False, debug=False):
        """
        initialize the bushy cell, using the default parameters for guinea pig from
        R&M2003, as a type II cell.
        Modifications to the cell can be made by calling methods below.
        """
        super(Bushy, self).__init__()

        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': 'guineapig-bushy-II', 'ttx': ttx}
        self.e_k = -70  # potassium reversal potential, mV
        self.e_na = 55
        self.e_h = -43
        self.c_m = 0.9  # specific membrane capacitance,  uf/cm^2
        self.R_a = 150  # axial resistivity of cytoplasm/axoplasm, ohm.cm
        self.vm0 = -63.6467358   # nominal for type II
        self.i_test_range=(-0.5, 0.5, 0.05)
        self.e_leak = -65

        soma = h.Section()  # one compartment of about 29000 um2
        soma.nseg = 1
        soma.insert('klt')
        soma.insert('kht')
        soma.insert('ihvcn')
        soma.insert('leak')
        soma.ek = self.e_k
        soma().ihvcn.eh = self.e_h
        soma().leak.erev = self.e_leak

        # insert the requested sodium channel
        if nach == 'jsrna':
            soma.insert('jsrna')
        elif nach == 'nav11':
            soma.insert('nav11')
        elif nach in ['na', 'nacn']:
            soma.insert('na')
        else:
            raise ValueError('Sodium channel %s not available for Bushy cells' % nach)

        soma.ena = self.e_na
        self.soma = soma
        self.species_scaling()  # set the default type II cell parameters
        self.all_sections['soma'].append(soma)
        self.add_section(soma)
        if debug:
            print "<< bushy: JSR bushy cell model created >>"

    def find_i0(self, vrange=[-70., -55.]):
        v0 = scipy.optimize.brentq(self.i_bushy, vrange[0], vrange[1])
        #print 'found V0 = %f' % v0
        return v0

    def i_bushy(self, V):
        """
        For the steady-state case, return the total current at voltage V
        Used to find the zero current point
        vrange brackets the interval

        """
        ix = {}

        # klt
        zss = 0.5   # <0,1>   : steady state inactivation of glt
        q10g = 2.0
        qg = q10g**((h.celsius-22.)/10.)
        winf = (1.0 / (1.0 + np.exp(-(V + 48.) / 6 )))**0.25
        zinf = zss + ((1.0-zss) / (1.0 + np.exp((V + 71.) / 10 )))
        gklt = qg*self.soma().klt.gbar*(winf**4.)*zinf
        ix['klt'] = gklt*(V - self.e_k)

        # ihvcn
        q10g = 2.0
        qg = q10g**((h.celsius-22.)/10.)
        rinf = 1. / (1+np.exp((V + 76.) / 7.))
        gh = qg*self.soma().ihvcn.gbar*rinf
        ix['ih'] = gh*(V - self.e_h)

        # kht
        q10g = 2.0
        qg = q10g**((h.celsius-22)/10)
        nf = 0.85  # proportion of n vs p kinetics
        ninf =   (1 + np.exp(-(V + 15.) / 5.))**(-0.5)
        pinf =  1 / (1 + np.exp(-(V + 23.) / 6.))
        gkht = qg*self.soma().kht.gbar*(nf*(ninf**2) + (1.-nf)*pinf)
        ix['kht'] = gkht*(V - self.e_k)

        # nacn
        q10g = 2.0
        qg = q10g**((h.celsius-22.)/10.)
        minf = 1 / (1+np.exp(-(V + 38.) / 7.))
        hinf = 1 / (1+np.exp((V + 65.) / 6.))
        gna = qg*self.soma().na.gbar*(minf**3.0)*hinf
        ix['na'] = gna*(V - self.e_na)

        # leak
        ix['leak'] = self.soma().leak.gbar*(V - self.e_leak)
#        print ix
        return np.sum([ix[i] for i in ix])

    def set_soma_size_from_Cm(self, cap):
        self.totcap = cap
        self.somaarea = self.totcap * 1E-6 / self.c_m  # pf -> uF, cm = 1uf/cm^2 nominal
        lstd = 1E4 * ((self.somaarea / np.pi) ** 0.5)  # convert from cm to um
        self.soma.diam = lstd
        self.soma.L = lstd

    def species_scaling(self, species='guineapig-bushy-II'):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        """
        soma = self.soma
        if species == 'mouse':
            # use conductance levels from Cao et al.,  J. Neurophys., 2007.
            #print 'Mouse bushy cell'
            self.set_soma_size_from_Cm(26.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(58.0, self.somaarea)
            soma().klt.gbar = nstomho(80.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(30.0, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            self.vm0 = -63.6
            self.axonsf = 0.57
        elif species == 'guineapig-bushy-II':
            self.set_soma_size_from_Cm(12.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().klt.gbar = nstomho(200.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(20.0, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            self.vm0 = self.find_i0() # -63.6467358  # "-63.6"
            self.axonsf = 0.57
            print 'set cell as: ', species
        elif species == 'guineapig-bushy-II-I':
            # guinea pig data from Rothman and Manis, 2003, type II=I
            self.set_soma_size_from_Cm(12.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().klt.gbar = nstomho(35.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(3.5, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            self.vm0 = -63.648175
            self.axonsf = 0.57
        elif species == 'cat':  # a cat is a big guinea pig
            self.set_soma_size_from_Cm(35.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().klt.gbar = nstomho(200.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(20.0, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            self.vm0 = -63.6
            self.axonsf = 1.0
        else:
            raise ValueError('Species %s is not recognized for Bushy cells', species)

        self.status['species'] = species

    def adjust_na_chans(self, soma, debug=False):
        """
        adjust the sodium channel conductance
        :param soma: a soma object whose sodium channel complement will have it's 
        conductances adjusted depending on the channel type
        :return nothing:
        """
        if self.status['ttx']:
            gnabar = 0.0
        else:
            gnabar = nstomho(1000.0, self.somaarea)
        nach = self.status['na']
        if nach == 'jsrna':
            soma().jsrna.gbar = gnabar
            soma.ena = self.e_na
            if debug:
                print 'jsrna gbar: ', soma().jsrna.gbar
        elif nach == 'nav11':
            soma().nav11.gbar = gnabar * 0.5
            soma.ena = self.e_na
            soma().nav11.vsna = 4.3
            if debug:
                print "bushy using inva11"
            print 'nav11 gbar: ', soma().nav11.gbar
        elif nach in ['nacn', 'na']:
            soma().na.gbar = gnabar
            soma.ena = self.e_na
            if debug:
                print 'nacn gbar: ', soma().na.gbar
        else:
            raise ValueError('Sodium channel %s is not recognized for Bushy cells', nach)

    def add_axon(self):
        Cell.add_axon(self, self.c_m, self.R_a, self.axonsf)

    def add_pumps(self):
        """
        Insert mechanisms for potassium ion, sodium ion, and a
        sodium-potassium pump at the soma.
        """
        soma = self.soma
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
        self.status['pumps'] = True

    def add_dendrites(self, debug=False):
        """
        Add a simple dendrite to the bushy cell.
        """
        if debug:
            print 'Adding dendrite to Bushy model'
        section = h.Section
        maindend = section(cell=self.soma)
        maindend.connect(self.soma)
        maindend.nseg = 10
        maindend.L = 100.0
        maindend.diam = 2.5
        maindend.insert('klt')
        maindend.insert('ihvcn')
        maindend().klt.gbar = self.soma().klt.gbar / 2.0
        maindend().ihvcn.gbar = self.soma().ihvcn.gbar / 2.0

        maindend.cm = self.c_m
        maindend.Ra = self.R_a
        nsecd = range(0, 5)
        secdend = []
        for ibd in nsecd:
            secdend.append(section(cell=self.soma))
        for ibd in nsecd:
            secdend[ibd].connect(maindend)
            secdend[ibd].diam = 1.0
            secdend[ibd].L = 15.0
            secdend[ibd].cm = self.c_m
            secdend[ibd].Ra = self.R_a
        self.maindend = maindend
        self.secdend = secdend
        self.status['dendrite'] = True
        if debug:
            print 'Bushy: added dendrite'
            h.topology()
        self.all_sections['maindend'].extend(maindend)
        self.all_sections['secdend'].extend(secdend)

