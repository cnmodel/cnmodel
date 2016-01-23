from neuron import h

from .cell import Cell
from .. import synapses
from ..util import nstomho
from .. import data

__all__ = ['Bushy', 'BushyRothman']


class Bushy(Cell):
    
    type = 'bushy'

    @classmethod
    def create(cls, model='RM03', **kwds):
        if model == 'RM03':
            return BushyRothman(**kwds)
        else:
            raise ValueError ('Bushy model %s is unknown', model)

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


class BushyRothman(Bushy):
    """
    VCN bushy cell model.
    Rothman and Manis, 2003abc (Type II, Type II-I)
    """

    def __init__(self, morphology=None, decorator=None, nach='na',
            ttx=False, species='guineapig', modelType=None, debug=False):
        """
        Initialize the bushy cell, using the default parameters for guinea pig from
        R&M2003, as a type II cell.
        Additional modifications to the cell can be made by calling methods below.
        
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
        super(BushyRothman, self).__init__()
        print "<< Bushy model: Creating point cell using JSR parameters >>"

        if modelType == None:
            modelType = 'II'
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'ttx': ttx, 'name': 'Bushy',
                       'morphology': morphology, 'decorator': decorator}
        self.i_test_range=(-0.5, 0.5, 0.05)
        self.spike_threshold = -40
        self.vrange = [-70., -57.]  # set a default vrange for searching for rmp
        
        if morphology is None:
            """
            instantiate a basic soma-only ("point") model
            """
            soma = h.Section(name="Bushy_Soma_%x" % id(self))  # one compartment of about 29000 um2
            soma.nseg = 1
        else:
            """
            instantiate a structured model with the morphology as specified by 
            the morphology file
            """
            soma = self.morphology_from_hoc(morphology=morphology, somasection='sections[0]')

        # decorate the morphology with ion channels
        if decorator is None:   # basic model, only on the soma
            self.mechanisms = ['klt', 'kht', 'ihvcn', 'leak', nach]
            for mech in self.mechanisms:
                soma.insert(mech)
            soma.ena = self.e_na
            soma.ek = self.e_k
            soma().ihvcn.eh = self.e_h
            soma().leak.erev = self.e_leak
            self.add_section(soma, 'soma')
            self.species_scaling(silent=True, species=species, modelType=modelType)  # set the default type II cell parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorated = decorator(self.hr, cellType='Bushy', modelType=modelType,
                                 parMap=None)
            self.decorated.channelValidate(self.hr, verify=False)
            self.mechanisms = self.decorated.hf.mechanisms  # copy out all of the mechanisms that were inserted
#        print 'Mechanisms inserted: ', self.mechanisms
        self.get_mechs(soma)
        self.cell_initialize(vrange=self.vrange)
        if debug:
            print "<< Bushy model created, point cell using JSR parameters >>"

    def species_scaling(self, species='guineapig', modelType='II', silent=True):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        Used only for point models.
        
        Parameters
        ----------
        species : string (default: 'guineapig')
            name of the species to use for scaling the conductances in the base point model
            Must be one of mouse, cat, guineapig
        
        modelType: string (default: 'II')
            definition of model type from RM03 models, type II or type II-I
        
        silent : boolean (default: True)
            run silently (True) or verbosely (False)
        """
        #print '\nSpecies scaling: %s   %s' % (species, type)
        knownspecies = ['mouse', 'guineapig', 'cat']
        soma = self.soma
        if species == 'mouse' and modelType == 'II':
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
        elif species == 'guineapig' and modelType =='II':
            self.set_soma_size_from_Cm(12.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().klt.gbar = nstomho(200.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(20.0, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            self.axonsf = 0.57
        elif species == 'guineapig' and modelType =='II-I':
            # guinea pig data from Rothman and Manis, 2003, type II=I
            self.i_test_range=(-0.4, 0.4, 0.02)
            self.set_soma_size_from_Cm(12.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().klt.gbar = nstomho(35.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(3.5, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            self.axonsf = 0.57
        elif species == 'cat' and modelType == 'II':  # a cat is a big guinea pig
            self.set_soma_size_from_Cm(35.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().klt.gbar = nstomho(200.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(20.0, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            self.axonsf = 1.0
        else:
            errmsg = 'Species "%s" or species-type "%s" is not recognized for Bushy cells' %  (species, modelType)
            errmsg += 'Valid species are: \n'
            for s in knownspecies:
                errmsg += '   %s\n' % s
            errmsg += '-'*40
            raise ValueError(errmsg)
        self.status['species'] = species
        self.status['modelType'] = modelType
        if not silent:
           print ' set cell as: ', species
           print ' with Vm rest = %6.3f' % self.vm0


       # print 'Rescaled, status: ', self.status

    def adjust_na_chans(self, soma, gbar=1000., debug=False):
        """
        adjust the sodium channel conductance
        Parameters
        ----------
        soma : neuron section object
            a soma object whose sodium channel complement will have it's 
            conductances adjusted depending on the channel type
        
        gbar : float (default: 1000.)
            the maximal conductance for the sodium channel
        
        debug : boolean (false):
            verbose printing
            
        Returns
        -------
            nothing
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

