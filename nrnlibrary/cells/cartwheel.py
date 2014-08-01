from neuron import h
import neuron as nrn
from ..util import nstomho
import numpy as np
from .cell import Cell

__all__ = ['Cartwheel', 'CartwheelDefault']

class Cartwheel(Cell):

    @classmethod
    def create(cls, model='RM03', **kwds):
        if model == 'RM03':
            return CartwheelDefault(**kwds)
        else:
            raise ValueError ('DStellate type %s is unknown', type)

class CartwheelDefault(Cartwheel, Cell):
    """
    DCN cartwheel cell model.
    
    """
    def __init__(self, debug=False, ttx=False, nach='naRsg', species='rat', type=None):
        super(CartwheelDefault, self).__init__()
        if type == None:
            type = 'I'
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'type': type, 'ttx': ttx, 'name': 'Cartwheel'}

        self.i_test_range=(-0.2, 0.2, 0.02)

        soma = h.Section(name="Cartwheel_Soma_%x" % id(self)) # one compartment of about 29000 um2
        cm = 1
        v_potassium = -80       # potassium reversal potential
        v_sodium = 50           # sodium reversal potential

        self.mechanisms = ['naRsg', 'bkpkj', 'hpkj', 'kpkj', 'kpkj2',
                           'kpkjslow', 'kpksk', 'lkpkj', 'cap']
        for mech in self.mechanisms:
            soma.insert(mech)
        soma.insert('cadiff')

       # soma().kpksk.gbar = 0.002
       # soma().lkpkj.gbar = 3e-4
        self.add_section(soma, 'soma')
        self.species_scaling(silent=True, species=species, type=type)  # set the default type I-c  cell parameters
        self.get_mechs(soma)
        if debug:
            print "<< cartwheel: Raman Purkinje cell model (modified) created >>"

    def species_scaling(self, silent=True, species='rat', type='I'):
        soma = self.soma
        dia = 18.
        self.set_soma_size_from_Diam(dia)# if species == 'rat' and type == 'I':
        #self.print_mechs(self.soma)
        #     self.set_soma_size_from_Cm(12.0)
        self.soma().bkpkj.gbar = nstomho(2., self.somaarea) # 2030
        self.soma().hpkj.gbar = nstomho(5, self.somaarea) # 29
        self.soma().kpkj.gbar = nstomho(100, self.somaarea) # 1160
        self.soma().kpkj2.gbar = nstomho(50, self.somaarea) #579
        self.soma().kpkjslow.gbar = nstomho(150, self.somaarea) # 1160
        self.soma().kpksk.gbar = nstomho(25, self.somaarea) # 2900
        self.soma().lkpkj.gbar = nstomho(5, self.somaarea) # 14.5
        self.soma().naRsg.gbar = nstomho(500, self.somaarea)  # 4340
        self.soma().cap.pcabar = 0.00015 # * (dia/2.)**2/(21./2.)**2 # 0.00005
        self.soma().ena = 50
        self.soma().ek = -80
        self.soma().lkpkj.e = -65
        self.soma().hpkj.eh = -43
        self.soma().eca = 50
        # else:
        #     raise ValueError('Species %s or species-type %s is not recognized for T-stellate cells' % (species, type))

        self.status['species'] = species
        self.status['type'] = type
        self.cell_initialize(showinfo=False)
        if not silent:
            print 'set cell as: ', species
            print ' with Vm rest = %f' % self.vm0

    def i_currents(self, V):
        """
        For the steady-state case, return the total current at voltage V
        Used to find the zero current point
        vrange brackets the interval
        Overrides i_currents in cells.py
        """
        for part in self.all_sections.keys():
            for sec in self.all_sections[part]:
                sec.v = V
        h.finitialize()
        self.ix = {}

        if 'naRsg' in self.mechanisms:
             self.ix['naRsg'] = self.soma().naRsg.gna*(V - self.soma().ena)
        if 'cap' in self.mechanisms:
            a = self.soma().cap.pcabar*self.soma().cap.minf
            self.ix['cap'] = a * self.ghk(V, self.soma().cao, self.soma().cai, 2)
        if 'kpkj' in self.mechanisms:
             self.ix['kpkj'] = self.soma().kpkj.gk*(V - self.soma().ek)
        if 'kpkj2' in self.mechanisms:
             self.ix['kpkj2'] = self.soma().kpkj2.gk*(V - self.soma().ek)
        if 'kpkjslow' in self.mechanisms:
             self.ix['kpkjslow'] = self.soma().kpkjslow.gk*(V - self.soma().ek)
        if 'kpksk' in self.mechanisms:
             self.ix['kdpyr'] = self.soma().kpksk.gk*(V - self.soma().ek)
        if 'bkpkj' in self.mechanisms:
             self.ix['bkpkj'] = self.soma().bkpkj.gbkpkj*(V - self.soma().ek)
        if 'hpkj' in self.mechanisms:
             self.ix['hpkj'] = self.soma().hpkj.gh*(V - self.soma().hpkj.eh)
        # leak
        if 'lkpkj' in self.mechanisms:
            self.ix['lkpkj'] = self.soma().lkpkj.gbar*(V - self.soma().lkpkj.e)
        return np.sum([self.ix[i] for i in self.ix])

    def ghk(self, v, ci, co, z):
        F = 9.6485e4  # (coul)
        R = 8.3145 # (joule/degC)
        T = h.celsius + 273.19  # Kelvin
        E = (1e-3) * v
        Ci = ci + (self.soma().cap.monovalPerm) * (self.soma().cap.monovalConc)  #       : Monovalent permeability
        if (np.fabs(1-np.exp(-z*(F*E)/(R*T))) < 1e-6):  #denominator is small -> Taylor series
            ghk = (1e-6) * z * F * (Ci-co*np.exp(-z*(F*E)/(R*T)))*(1-(z*(F*E)/(R*T)))
        else:
            ghk = (1e-6) * z**2.*(E*F**2.)/(R*T)*(Ci-co*np.exp(-z*(F*E)/(R*T)))/(1-np.exp(-z*(F*E)/(R*T)))
        return ghk
