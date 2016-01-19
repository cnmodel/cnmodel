from neuron import h
import neuron as nrn
import numpy as np
import scipy.optimize

from .cell import Cell
from .. import synapses
from ..util import nstomho
from .. import data

__all__ = ['Bushy', 'BushyRothman', 'BushyHoc']


class Bushy(Cell):
    
    type = 'bushy'

    @classmethod
    def create(cls, model='RM03', **kwds):
        if model == 'RM03':
            return BushyRothman(**kwds)
        if model in ['hoc']:
            return BushyHoc(**kwds)
        else:
            raise ValueError ('Bushy type %s is unknown', type)

    def make_psd(self, terminal, **kwds):
        from .. import cells
        
        pre_sec = terminal.section
        pre_cell = terminal.cell
        post_sec = self.soma
        
        if isinstance(pre_cell, cells.SGC):
            # Max conductances for the glu mechanisms are calibrated by 
            # running `synapses/tests/test_psd.py`. The test should fail
            # if these values are incorrect:
            AMPA_gmax = 3.314707700918133*1e3  # factor of 1e3 scales to pS (.mod mechanisms) from nS.
            NMDA_gmax = 0.4531929783503451*1e3
            
            # Get AMPAR kinetic constants from database 
            params = data.get('sgc_ampa_kinetics', species='mouse', post_type='bushy',
                              field=['Ro1', 'Ro2', 'Rc1', 'Rc2', 'PA'])
            for key, value in kwds.iteritems():
                        print "make_psd: keywords  %s == %s" %(key,value)
            return synapses.GluPSD(post_sec, terminal,
                                   ampa_gmax=AMPA_gmax,
                                   nmda_gmax=NMDA_gmax,
                                   ampa_params=dict(
                                        Ro1=params['Ro1'],
                                        Ro2=params['Ro2'],
                                        Rc1=params['Rc1'],
                                        Rc2=params['Rc2'],),
                                   **kwds)
        elif isinstance(pre_cell, cells.DStellate):
            # Get GLY kinetic constants from database 
            params = data.get('gly_kinetics', species='mouse', post_type='bushy',
                              field=['KU', 'KV', 'XMax'])
            return synapses.GlyPSD(post_sec, terminal, params=params,
                                   psdType='glyslow', **kwds)
        else:
            raise TypeError("Cannot make PSD for %s => %s" % 
                            (pre_cell.__class__.__name__, 
                             self.__class__.__name__))

class BushyHoc(Bushy):
    """
    VCN bushy cell model - with dendritic structure from a hoc file
    The model dendritic structure will have been read in HocReader
    The model structures should be passed to the first parameter, hoc
    Ion channels are likely alread decorated, so this routine only inserts the
    model sections into the scaffolding of the cell class, so we can use that class.
    
    Based on Rothman and Manis, 2003abc (Type II, Type II-I)
    """

    def __init__(self, hoc=None, nach='na', ttx=False, debug=False, species='guineapig', modeltype=None, decorator=None):
        """
        initialize the bushy cell, using the default parameters for guinea pig from
        R&M2003, as a type II cell.
        Modifications to the cell can be made by calling methods below.
        """
        super(BushyHoc, self).__init__()
        print "<< bushy: Creating Bushy model from HOC file >>"
        if modeltype == None:
            modeltype = 'II'
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'type': type, 'ttx': ttx, 'name': 'Bushy'}
        self.i_test_range=(-0.5, 0.5, 0.05)
        self.spike_threshold = -40
        self.vrange = [-70., -57.]  # set a default vrange for searching for rmp

        self.mechanisms = ['klt', 'kht', 'ihvcn', 'leak', nach]
        nach = self.status['na']
        if decorator is not None:
            self.cd = decorator(hoc, celltype='Bushy', modeltype=modeltype,
                             parMap=None)
        hf = self.cd.hf  # get the hoc information
        self.cd.channelValidate(hf, verify=False)
        hf._read_section_info()  # this should populate the section channel densities as well from an neurovis read file
            # these were not instantiated when the file was read, but when the decorator was run.
        for s in hf.sec_groups.keys():
            for sec in hf.sec_groups[s]:
                section = hf.get_section(sec)
                mechs = hf.get_mechanisms(sec)
                self.add_section(section, s) # add the section to the cell.
               # print '\nmechanisms for section: %s', section
               # self.print_mechs(section)
        soma = hf.get_section('sections[0]')
        self.set_soma_size_from_Section(soma)  # this is used for reporting and setting g values...
        self.cell_initialize(vrange=self.vrange)
#        if debug:
        print "<< bushy: Bushy model created from HOC file >>"
        #print 'Cell created: ', self.status


class BushyRothman(Bushy):
    """
    VCN bushy cell model.
    Rothman and Manis, 2003abc (Type II, Type II-I)
    """

    def __init__(self, nach='na', ttx=False, debug=False, species='guineapig', type=None):
        """
        initialize the bushy cell, using the default parameters for guinea pig from
        R&M2003, as a type II cell.
        Modifications to the cell can be made by calling methods below.
        """
        super(BushyRothman, self).__init__()
        print "<< Bushy model: Creating point cell using JSR parameters >>"

        if type == None:
            type = 'II'
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'type': type, 'ttx': ttx, 'name': 'Bushy'}
        self.i_test_range=(-0.5, 0.5, 0.05)
        self.spike_threshold = -40
        self.vrange = [-70., -57.]  # set a default vrange for searching for rmp

        soma = h.Section(name="Bushy_Soma_%x" % id(self))  # one compartment of about 29000 um2
        soma.nseg = 1

        self.mechanisms = ['klt', 'kht', 'ihvcn', 'leak', nach]
        for mech in self.mechanisms:
            soma.insert(mech)
        soma.ena = self.e_na
        soma.ek = self.e_k
        soma().ihvcn.eh = self.e_h
        soma().leak.erev = self.e_leak
        self.add_section(soma, 'soma')
        self.species_scaling(silent=True, species=species, type=type)  # set the default type II cell parameters
        self.get_mechs(soma)
        self.cell_initialize(vrange=self.vrange)
        if debug:
            pass
        print "<< Bushy model created, point cell using JSR parameters >>"
        #print 'Cell created: ', self.status

    def species_scaling(self, species='guineapig', type='II', silent=True):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        """
        #print '\nSpecies scaling: %s   %s' % (species, type)
        soma = self.soma
        if species == 'mouse' and type == 'II':
            # use conductance levels from Cao et al.,  J. Neurophys., 2007.
           # print 'Mouse bushy cell'
            self.set_soma_size_from_Cm(26.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(58.0, self.somaarea)
            soma().klt.gbar = nstomho(80.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(30.0, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            self.vrange = [-70., -55.]  # need to specify non-default range for convergence
            self.axonsf = 0.57
        elif species == 'guineapig' and type =='II':
            self.set_soma_size_from_Cm(12.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().klt.gbar = nstomho(200.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(20.0, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            self.axonsf = 0.57
        elif species == 'guineapig' and type =='II-I':
            # guinea pig data from Rothman and Manis, 2003, type II=I
            self.i_test_range=(-0.4, 0.4, 0.02)
            self.set_soma_size_from_Cm(12.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().klt.gbar = nstomho(35.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(3.5, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            self.axonsf = 0.57
        elif species == 'cat' and type == 'II':  # a cat is a big guinea pig
            self.set_soma_size_from_Cm(35.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().klt.gbar = nstomho(200.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(20.0, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            self.axonsf = 1.0
        else:
            raise ValueError('Species "%s" or species-type "%s" is not recognized for Bushy cells' %  (species, type))
        self.status['species'] = species
        self.status['type'] = type
#        self.cell_initialize(showinfo=False)
#        if not silent:
#            print 'set cell as: ', species
#            print ' with Vm rest = %6.3f' % self.vm0


       # print 'Rescaled, status: ', self.status

    def adjust_na_chans(self, soma, gbar=1000., debug=False):
        """
        adjust the sodium channel conductance
        :param soma: a soma object whose sodium channel complement will have it's 
        conductances adjusted depending on the channel type
        :return nothing:
        """
        if self.status['ttx']:
            gnabar = 0.0
        else:
            gnabar = nstomho(gbar, self.somaarea)
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
        elif nach in ['na', 'nacn']:
            soma().na.gbar = gnabar
            soma.ena = self.e_na
            if debug:
                print 'na gbar: ', soma().na.gbar
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
        self.add_section(maindend, 'maindend')
        self.add_section(secdend, 'secdend')

