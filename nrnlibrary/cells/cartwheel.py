from neuron import h
import neuron as nrn
from ..util import nstomho

from .cell import Cell

__all__ = ['Cartwheel']


class Cartwheel(Cell):
    """
    DCN cartwheel cell model.
    
    """
    def __init__(self, debug=False, ttx=False, nach='naRsg', species='rat', type=None):
        super(Cartwheel, self).__init__()
        if type == None:
            type = 'I'
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'type': type, 'ttx': ttx, 'name': 'Cartwheel'}

        self.i_test_range=(-0.15, 0.15, 0.01)

        soma = h.Section(name="Cartwheel_Soma_%x" % id(self)) # one compartment of about 29000 um2
        cm = 1
        v_potassium = -80       # potassium reversal potential
        v_sodium = 50           # sodium reversal potential

        soma.insert('naRsg')
        soma.insert('kpkj')
        soma.insert('kpkj2')
        soma.insert('kpkjslow')
        soma.insert('bkpkj')
        soma.insert('kpksk')
        soma.insert('cadiff')
        soma.insert('cap')
        soma.insert('lkpkj')
        soma.insert('hpkj')
        soma.ena = v_sodium
        soma.ek = v_potassium
        soma().kpksk.gbar = 0.002
        self.mechanisms = ['naRsg', 'bkpkj', 'hpkj', 'kpkj', 'kpkj2',
                           'kpkjslow', 'kpksk', 'pkjlk',]
        self.add_section(soma, 'soma')
        self.species_scaling(silent=True, species=species, type=type)  # set the default type I-c  cell parameters
        self.get_mechs(soma)
        if debug:
            print "<< cartwheel: Raman Purkinje cell model (modified) created >>"

    def species_scaling(self, silent=True, species='rat', type='I'):
        soma = self.soma
        self.set_soma_size_from_Diam(50.0)# if species == 'rat' and type == 'I':
        #     self.set_soma_size_from_Cm(12.0)
        #     soma().napyr.gbar = nstomho(350, self.somaarea)
        #     #soma().pyr.gnapbar = 0.0
        #     soma().kdpyr.gbar = nstomho(80, self.somaarea) # used to be 20?
        #     soma().kif.gbar = nstomho(150, self.somaarea)
        #     soma().kis.gbar = nstomho(40, self.somaarea)
        #     soma().ihpyr.gbar = nstomho(3, self.somaarea)
        #     soma().leak.gbar = nstomho(2.8, self.somaarea)
        #     soma().leak.erev = -57.7  # override default values in cell.py
        #     soma().ena = 50.0
        #     soma().ek = -81.5
        #     soma().eh = -43

        # else:
        #     raise ValueError('Species %s or species-type %s is not recognized for T-stellate cells' % (species, type))

        self.status['species'] = species
        self.status['type'] = type
        self.cell_initialize(showinfo=True)
        if not silent:
            print 'set cell as: ', species
            print ' with Vm rest = %f' % self.vm0
