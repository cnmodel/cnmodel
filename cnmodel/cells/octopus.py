from neuron import h
from ..util import nstomho
import numpy as np
"""
Original hoc code from RMmodel.hoc
// including the "Octopus" cell:
proc set_Type2o() {
    gbar_na = nstomho(1000)
    gbar_kht = nstomho(150)
    gbar_klt = nstomho(600)
    gbar_ka = nstomho(0)
    gbar_ih = nstomho(0)
    gbar_hcno = nstomho(40)
    gbar_leak = nstomho(2)
    model = 6
    modelname = "Type IIo (Octopus)"
    vm0 = -66.67
}

"""

from .cell import Cell

__all__ = ['Octopus', 'OctopusRothman']

class Octopus(Cell):

    @classmethod
    def create(cls, model='RM03', **kwds):
        if model == 'RM03':
            print 'making RM03'
            return OctopusRothman(**kwds)
        else:
            raise ValueError ('DStellate type %s is unknown', type)

class OctopusRothman(Octopus, Cell):
    """
    VCN octopus cell model (point cell).
    Rothman and Manis, 2003abc (Type II, with high gklt and hcno - octopus cell h current).
    """

    def __init__(self, morphology=None, decorator=None, nach='jsrna', ttx=False,
                debug=False, species='guineapig', modelType=None):
        """
        initialize the octopus cell, using the default parameters for guinea pig from
        R&M2003, as a type II cell with modified conductances.
        Modifications to the cell can be made by calling methods below.
        
        Parameters
        ----------
        morphology : string (default: None)
            a file name to read the cell morphology from. If a valid file is found, a cell is constructed
            as a cable model from the hoc file.
            If None (default), the only a point model is made, exactly according to RM03.
            
        decorator : Python function (default: None)
            decorator is a function that "decorates" the morphology with ion channels according
            to a set of rules.
            If None, a default set of channels aer inserted into the first soma section, and the
            rest of the structure is "bare".
        
        nach : string (default: 'na')
            nach selects the type of sodium channel that will be used in the model. A channel mechanims
            by that name must exist. 
        
        ttx : Boolean (default: False)
            If ttx is True, then the sodium channel conductance is set to 0 everywhere in the cell.
            Currently, this is not implemented.
        
        species: string (default 'guineapig')
            species defines the channel density that will be inserted for different models. Note that
            if a decorator function is specified, this argument is ignored.
            
        modelType: string (default: None)
            modelType specifies the type of the model that will be used (e.g., "II", "II-I", etc).
            modelType is passed to the decorator, or to species_scaling to adjust point models.
            
        debug: boolean (default: False)
            debug is a boolean flag. When set, there will be multiple printouts of progress and parameters.
            
        Returns
        -------
            Nothing
        """
        
        super(OctopusRothman, self).__init__()
        if modelType == None:
            modelType = 'II-o'
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'type': modelType, 'ttx': ttx, 'name': 'Octopus',
                        'morphology': morphology, 'decorator': decorator}
        self.i_test_range=(-4.0, 4.0, 0.2)
        self.spike_threshold = -50
        if morphology is None:
            """
            instantiate a basic soma-only ("point") model
            """
            soma = h.Section(name="Octopus_Soma_%x" % id(self))  # one compartment of about 29000 um2
            soma.nseg = 1
            self.add_section(soma, 'soma')
        else:
            """
            instantiate a structured model with the morphology as specified by 
            the morphology file
            """
            self.set_morphology(morphology_file=morphology)

        # decorate the morphology with ion channels
        if decorator is None:   # basic model, only on the soma
            self.e_leak = -62
            self.e_h = -38
            self.R_a = 100
            self.mechanisms = ['klt', 'kht', 'hcno', 'leak', nach]
            for mech in self.mechanisms:
                self.soma.insert(mech)
            self.soma.ek = self.e_k
            self.soma.ena = self.e_na
            self.soma().hcno.eh = self.e_h
            self.soma().leak.erev = self.e_leak
            self.soma.Ra = self.R_a
            self.species_scaling(silent=True, species=species, modelType=modelType)  # set the default type II cell parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorate()
#        print 'Mechanisms inserted: ', self.mechanisms
        self.get_mechs(self.soma)
        self.cell_initialize()
        
        if debug:
            print "<< octopus: octopus cell model created >>"
        #print 'Cell created: ', self.status

    def species_scaling(self, species='guineapig', modelType='II-o', silent=True):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        """
        #print '\nSpecies scaling: %s   %s' % (species, type)
        soma = self.soma
        # if species == 'mouse' and type == 'II-o':
        #     # use conductance levels from Cao et al.,  J. Neurophys., 2007.
        #     #print 'Mouse octopus cell'
        #     self.set_soma_size_from_Cm(26.0)
        #     self.adjust_na_chans(soma)
        #     soma().kht.gbar = nstomho(58.0, self.somaarea)
        #     soma().klt.gbar = nstomho(80.0, self.somaarea)
        #     soma().hcno.gbar = nstomho(30.0, self.somaarea)
        #     soma().leak.gbar = nstomho(2.0, self.somaarea)
        #     self.vm0 = self.find_i0()
        #     self.axonsf = 0.57
        if species == 'guineapig' and modelType =='II-o':
            self.set_soma_size_from_Cm(25.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = 0.0061  # nstomho(150.0, self.somaarea)  # 6.1 mmho/cm2
            soma().klt.gbar = 0.0407  # nstomho(3196.0, self.somaarea)  #  40.7 mmho/cm2
            soma().hcno.gbar = 0.0076  #nstomho(40.0, self.somaarea)  # 7.6 mmho/cm2, cf. Bal and Oertel, Spencer et al. 25 u dia cell
            soma().leak.gbar = 0.0005  # nstomho(2.0, self.somaarea)
            self.axonsf = 1.0
        # elif species == 'guineapig' and type =='II-I':
        #     # guinea pig data from Rothman and Manis, 2003, type II=I
        #     self.i_test_range=(-0.4, 0.4, 0.02)
        #     self.set_soma_size_from_Cm(12.0)
        #     self.adjust_na_chans(soma)
        #     soma().kht.gbar = nstomho(150.0, self.somaarea)
        #     soma().klt.gbar = nstomho(35.0, self.somaarea)
        #     soma().hcno.gbar = nstomho(3.5, self.somaarea)
        #     soma().leak.gbar = nstomho(2.0, self.somaarea)
        #     self.axonsf = 0.57
        # elif species == 'cat' and type == 'II':  # a cat is a big guinea pig
        #     self.set_soma_size_from_Cm(35.0)
        #     self.adjust_na_chans(soma)
        #     soma().kht.gbar = nstomho(150.0, self.somaarea)
        #     soma().klt.gbar = nstomho(200.0, self.somaarea)
        #     soma().hcno.gbar = nstomho(20.0, self.somaarea)
        #     soma().leak.gbar = nstomho(2.0, self.somaarea)
        #     self.axonsf = 1.0
        else:
            raise ValueError('Species "%s" or species-type "%s" is not recognized for octopus cells' %  (species, type))
        self.status['species'] = species
        self.status['modelType'] = modelType
        self.cell_initialize(showinfo=False)
        if not silent:
            print 'set cell as: ', species
            print ' with Vm rest = %6.3f' % self.vm0

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
            gnabar = 0.04244 # 4.2441  # nstomho(1000.0, self.somaarea)  # mmho/cm2 - 4244.1 moh - 4.2441
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
                print "octopus using inva11"
            print 'nav11 gbar: ', soma().nav11.gbar
        elif nach in ['na', 'nacn']:
            soma().na.gbar = gnabar
            soma.ena = self.e_na
            if debug:
                print 'na gbar: ', soma().na.gbar
        else:
            raise ValueError('Sodium channel %s is not recognized for octopus cells', nach)

    def i_currents(self, V):
        """
        For the steady-state case, return the total current at voltage V
        Used to find the zero current point
        vrange brackets the interval
        Implemented here are the basic RM03 mechanisms
        This function should be replaced for specific cell types.
        """
        for part in self.all_sections.keys():
            for sec in self.all_sections[part]:
                sec.v = V

        h.t = 0.
        h.finitialize()
        self.ix = {}
        if 'na' in self.mechanisms:
            #print dir(self.soma().na)
            self.ix['na'] = self.soma().na.gna*(V - self.soma().ena)
        if 'jsrna' in self.mechanisms:
            #print dir(self.soma().na)
            self.ix['jsrna'] = self.soma().jsrna.gna*(V - self.soma().ena)
        if 'klt' in self.mechanisms:
            self.ix['klt'] = self.soma().klt.gklt*(V - self.soma().ek)
        if 'kht' in self.mechanisms:
            self.ix['kht'] = self.soma().kht.gkht*(V - self.soma().ek)
        if 'hcno' in self.mechanisms:
            self.ix['hcno'] = self.soma().hcno.gh*(V - self.soma().hcno.eh)
        if 'leak' in self.mechanisms:
            self.ix['leak'] = self.soma().leak.gbar*(V - self.soma().leak.erev)
#        print self.status['name'], self.status['type'], V, self.ix
        return np.sum([self.ix[i] for i in self.ix])


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
        Add a simple dendrite to the octopus cell.
        """
        if debug:
            print 'Adding dendrite to octopus model'
        section = h.Section
        maindend = section(cell=self.soma)
        maindend.connect(self.soma)
        maindend.nseg = 10
        maindend.L = 100.0
        maindend.diam = 2.5
        maindend.insert('klt')
        maindend.insert('ihcno')
        maindend().klt.gbar = self.soma().klt.gbar / 2.0
        maindend().hcno.gbar = self.soma().hcno.gbar / 2.0

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
            print 'octopus: added dendrite'
            h.topology()
        self.add_section(maindend, 'maindend')
        self.add_section(secdend, 'secdend')

