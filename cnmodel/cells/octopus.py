from neuron import h
from ..util import nstomho
from ..util import Params
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

__all__ = ['Octopus', 'OctopusRothman', 'OctopusSpencer']

class Octopus(Cell):

    type = 'octopus'
    
    @classmethod
    def create(cls, modelType='RM03', **kwds):
        print 'modelType: ', modelType
        if modelType in ['RM03', 'II-o']:
            print 'making RM03'
            return OctopusRothman(**kwds)
        elif modelType == 'Spencer':
            return OctopusSpencer(**kwds)
        else:
            raise ValueError ('Octopus cell type %s is unknown', type)

    def make_psd(self, terminal, psd_type, **kwds):
        """
        Connect a presynaptic terminal to one post section at the specified location, with the fraction
        of the "standard" conductance determined by gbar.
        The default condition is to try to pass the default unit test (loc=0.5)
        
        Parameters
        ----------
        terminal : Presynaptic terminal (NEURON object)
        
        psd_type : either simple or multisite PSD for bushy cell
        
        kwds: dict of options. Two are currently handled:
        postsize : expect a list consisting of [sectionno, location (float)]
        AMPAScale : float to scale the ampa currents
        
        """
        if 'postsite' in kwds:  # use a defined location instead of the default (soma(0.5)
            postsite = kwds['postsite']
            loc = postsite[1]  # where on the section?
            uname = 'sections[%d]' % postsite[0]  # make a name to look up the neuron section object
            post_sec = self.hr.get_section(uname)  # Tell us where to put the synapse.
        else:
            loc = 0.5
            post_sec = self.soma
        
        if psd_type == 'simple':
            return self.make_exp2_psd(post_sec, terminal, loc=loc)
        elif psd_type == 'multisite':
            if terminal.cell.type == 'sgc':
                # Max conductances for the glu mechanisms are calibrated by 
                # running `synapses/tests/test_psd.py`. The test should fail
                # if these values are incorrect:
                AMPA_gmax = 3.314707700918133*1e3  # factor of 1e3 scales to pS (.mod mechanisms) from nS.
                NMDA_gmax = 0.4531929783503451*1e3
                if 'AMPAScale' in kwds:
                    AMPA_gmax = AMPA_gmax * kwds['AMPAScale']  # allow scaling of AMPA conductances
                if 'NMDAScale' in kwds:
                    NMDA_gmax = NMDA_gmax*kwds['NMDAScale']
                return self.make_glu_psd(post_sec, terminal, AMPA_gmax, NMDA_gmax, loc=loc)
            elif terminal.cell.type == 'dstellate':
                return self.make_gly_psd(post_sec, terminal, type='glyslow', loc=loc)
            else:
                raise TypeError("Cannot make PSD for %s => %s" % 
                            (terminal.cell.type, self.type))
        else:
            raise ValueError("Unsupported psd type %s" % psd_type)


class OctopusRothman(Octopus, Cell):
    """
    VCN octopus cell model (point cell).
    Rothman and Manis, 2003abc (Type II, with high gklt and hcno - octopus cell h current).
    """

    def __init__(self, morphology=None, decorator=None, nach='jsrna', ttx=False,
                species='guineapig', modelType=None, debug=False):
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
                       'na': nach, 'species': species, 'modelype': modelType, 'ttx': ttx, 'name': 'Octopus',
                        'morphology': morphology, 'decorator': decorator}
        self.i_test_range=(-4.0, 4.0, 0.2)
        self.spike_threshold = -50
        self.vrange = [-70., -57.]  # set a default vrange for searching for rmp
        
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
            self.e_leak = -73. # from McGinley et al., 2016
            self.e_h = -38. # from McGinley et al. 
            self.R_a = 195  # McGinley et al. 
            self.mechanisms = ['klt', 'kht', 'hcnobo', 'leak', nach]
            for mech in self.mechanisms:
                self.soma.insert(mech)
            self.soma.ek = self.e_k
            self.soma.ena = self.e_na
            self.soma().hcnobo.eh = self.e_h
            self.soma().leak.erev = self.e_leak
            self.soma.Ra = self.R_a
            self.species_scaling(silent=True, species=species, modelType=modelType)  # set the default type II cell parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorate()
#        print 'Mechanisms inserted: ', self.mechanisms
        self.get_mechs(self.soma)
        self.cell_initialize(vrange=self.vrange)
        print 'soma area: ', self.somaarea
        print 'soma cap: ', self.totcap
        
        if debug:
            print "<< octopus: octopus cell model created >>"
        #print 'Cell created: ', self.status

    def species_scaling(self, species='guineapig', modelType='II-o', silent=True):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        Used ONLY for point models.
        
        Parameters
        ----------
        species : string (default: 'guineapig')
            name of the species to use for scaling the conductances in the base point model
            Must be guineapig
        
        modelType: string (default: 'II-o')
            definition of model type from RM03 models, currently limited to type II-o
        
        silent : boolean (default: True)
            run silently (True) or verbosely (False)
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
            self.print_soma_info()
            
            # print 'soma area: ', self.somaarea
            # print 'soma cap: ', self.totcap
            # print 'soma L', self.soma.L
            # print 'diam: ', self.soma.diam
            # print 'cm: ', self.c_m
            self.adjust_na_chans(soma)
            soma().kht.gbar = 0.0061  # nstomho(150.0, self.somaarea)  # 6.1 mmho/cm2
            soma().klt.gbar = 0.0407  # nstomho(3196.0, self.somaarea)  #  40.7 mmho/cm2
            soma().hcnobo.gbar = 0.0076  #nstomho(40.0, self.somaarea)  # 7.6 mmho/cm2, cf. Bal and Oertel, Spencer et al. 25 u dia cell
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
        self.cell_initialize(showinfo=True)
        if not silent:
            print 'set cell as: ', species
            print ' with Vm rest = %6.3f' % self.vm0

    def adjust_na_chans(self, soma, debug=False):
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
        Nothing
        """
        if self.status['ttx']:
            gnabar = 0.0
        else:
            gnabar = 0.04244 # 4.2441  # nstomho(1000.0, self.somaarea)  # mmho/cm2 - 4244.1 moh - 4.2441
        nach = self.status['na']
        if nach == 'jsrna':
            soma().jsrna.gbar = gnabar * 0.2
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


class OctopusSpencer(Octopus, Cell):
    """
    VCN octopus cell model (with dendrites).
    Based on Spencer et al Front. Comput. Neurosci., 22 October 2012 | https://doi.org/10.3389/fncom.2012.00083
    """

    def __init__(self, morphology=None, decorator=None, nach='jsrna', ttx=False,
                species='guineapig', modelType=None, debug=False):
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
        
        super(OctopusSpencer, self).__init__()
        if modelType == None:
            modelType = 'Spencer'
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'ttx': ttx, 'name': 'Octopus',
                        'morphology': morphology, 'decorator': decorator}
        self.i_test_range=(-4.0, 4.0, 0.2)
        self.spike_threshold = -50
        self.vrange = [-70., -57.]  # set a default vrange for searching for rmp
        
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
            self.e_leak = -62. # from Spencer et al., 2012
            self.e_h = -38. ## from Spencer et al., 2012
            self.R_a = 100.  # from Spencer et al., 2012 
            self.mechanisms = ['klt', 'kht', 'hcnobo', 'leak', nach]
            for mech in self.mechanisms:
                self.soma.insert(mech)
            self.soma.ek = -70. # self.e_k
            self.soma.ena = 55.0 # self.e_na
            self.soma().hcnobo.eh = self.e_h
            self.soma().leak.erev = self.e_leak
            self.soma.Ra = self.R_a
            #self.species_scaling(silent=True, species=species, modelType=modelType)  # set the default type II cell parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorate()
            self.decorated.channelValidate(self, verify=True)
#        print 'Mechanisms inserted: ', self.mechanisms
        self.get_mechs(self.soma)
        self.cell_initialize(vrange=self.vrange)
        
        if debug:
            print "<< octopus: octopus cell model created >>"
        #print 'Cell created: ', self.status

    def channel_manager(self, modelType='Spencer'):
        """
        This routine defines channel density maps and distance map patterns
        for each type of compartment in the cell. The maps
        are used by the ChannelDecorator class (specifically, it's private
        \_biophys function) to decorate the cell membrane.
        
        Parameters
        ----------
        modelType : string (default: 'Spencer')
            A string that defines the type of the model. Currently, 1 type is implemented:
            Spencer : Spencer et al Front. Comput. Neurosci. 2012
        
        Returns
        -------
        Nothing
        
        Notes
        -----
        This routine defines the following variables for the class:
        # conductances (gBar)
        # a channelMap (dictonary of channel densities in defined anatomical compartments)
        # a current injection range for IV's (when testing)
        # a distance map, which defines how selected conductances in selected compartments
        will change with distance. This includes both linear and exponential gradients,
        the minimum conductance at the end of the gradient, and the space constant or
        slope for the gradient.
        
        """
        
        #
        # Create a model based on the Spencer model
        # 6
        if modelType == 'Spencer':
#            print self.c_m
            self.c_m = 0.9
            self.set_soma_size_from_Section(self.soma)
            totcap = self.totcap
            refarea = self.somaarea # totcap / self.c_m  # see above for units
            self.print_soma_info()
            
            self.gBar = Params(nabar=0., #0.0407,  # S/cm2
                               nabar_ais=0., # 0.42441,
                               khtbar=0.0061,
                               khtbar_hillock=0.0061,
                               khtbar_dend=0.0061,
                               kltbar=0.0407,
                               kltbar_dend=0.0027,
                               kltbar_hillock=0.0407,
                               ihbar=0.0076,
                               ihbar_dend=0.0006,
                               leakbar=0.0002,
            )
            
            self.channelMap = {
                # 'axon': {'nacn': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar, 'ihvcn': 0.,
                #          'leak': self.gBar.leakbar / 2.},
                'hillock': {'jsrna': 0., 'klt': self.gBar.kltbar_hillock, 'kht': self.gBar.khtbar_hillock, 'ihvcn': 0.,
                            'leak': self.gBar.leakbar, },
                # initial segment:
                'unmyelinatedaxon': {'jsrna': self.gBar.nabar_ais, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar,
                            'ihvcn': self.gBar.ihbar, 'leak': self.gBar.leakbar, },
                'soma': {'jsrna': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar,
                         'ihvcn': self.gBar.ihbar, 'leak': self.gBar.leakbar, },
                'primarydendrite': {'jsrna': self.gBar.nabar*0., 'klt': self.gBar.kltbar_dend, 'kht': self.gBar.khtbar * 0.,
                         'ihvcn': self.gBar.ihbar_dend, 'leak': self.gBar.leakbar, },
            }
            self.distMap = {'primardendrite': {'klt': {'gradient': 'flat', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'flat', 'gminf': 0., 'lambda': 100.},
                                     'ihvcn': {'gradient': 'flat', 'gminf': 0., 'lambda': 100.}}, # all flat with distance
                            }

        else:
            raise ValueError('model type %s is not implemented' % modelType)
