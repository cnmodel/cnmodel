from neuron import h
import numpy as np
import neuron as nrn

from .cell import Cell
#from .. import synapses
from ..util import nstomho
from ..util import Params
from .. import data

__all__ = ['TStellate', 'TStellateRothman', 'TStellateNav11'] 


class TStellate(Cell):
    
    type = 'tstellate'

    @classmethod
    def create(cls, model='RM03', **kwds):
        if model == 'RM03':  # original Rothman-Manis 2003, 22C, point cell, extendable
            return TStellateRothman(**kwds)
        elif model == 'Nav11':   # Xie-Manis, 2013, 37C, pointe cell, extendable
            return TStellateNav11(**kwds)
        else:
            raise ValueError ('TStellate type %s is unknown', type)
        

    def make_psd(self, terminal, psd_type, **kwds):
        """
        Connect a presynaptic terminal to one post section at the specified location, with the fraction
        of the "standard" conductance determined by gbar.
        The default condition is to try to pass the default unit test (loc=0.5)
        
        Scaling is corrected by initial release probability now.
        
        Parameters
        ----------
        terminal : Presynaptic terminal (NEURON object)
        
        psd_type : either simple or multisite PSD for bushy cell
        
        kwds: dict of options.
            Available options:
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
                # if these values are incorrect
                self.AMPAR_gmax = data.get('sgc_synapse', species=self.species,
                        post_type=self.type, field='AMPAR_gmax')*1e3
                self.NMDAR_gmax = data.get('sgc_synapse', species=self.species,
                        post_type=self.type, field='NMDAR_gmax')*1e3
                self.Pr = data.get('sgc_synapse', species=self.species,
                        post_type=self.type, field='Pr')
                # adjust gmax to correct for initial Pr
                self.AMPAR_gmax = self.AMPAR_gmax/self.Pr
                self.NMDAR_gmax = self.NMDAR_gmax/self.Pr
                # old values:
                # AMPA_gmax = 0.22479596944138733*1e3  # factor of 1e3 scales to pS (.mod mechanisms) from nS.
                # NMDA_gmax = 0.12281291946623739*1e3
                if 'AMPAScale' in kwds:
                    self.AMPAR_gmax = self.AMPAR_gmax * kwds['AMPAScale']  # allow scaling of AMPA conductances
                if 'NMDAScale' in kwds:
                    self.NMDAR_gmax = self.NMDAR_gmax*kwds['NMDAScale']
                return self.make_glu_psd(post_sec, terminal, self.AMPAR_gmax, self.NMDAR_gmax, loc=loc)
            elif terminal.cell.type == 'dstellate':
                return self.make_gly_psd(post_sec, terminal, type='glyfast', loc=loc)
            else:
                raise TypeError("Cannot make PSD for %s => %s" % 
                            (terminal.cell.type, self.type))
        else:
            raise ValueError("Unsupported psd type %s" % psd_type)


class TStellateRothman(TStellate):
    """
    VCN T-stellate base model.
    Rothman and Manis, 2003abc (Type I-c, Type I-t)
    """
    def __init__(self, morphology=None, decorator=None, nach=None,
                ttx=False,
                species='guineapig', modelType=None, debug=False):
        """
        Initialize a planar stellate (T-stellate) cell, using the default parameters for guinea pig from
        R&M2003, as a type I cell.
        Modifications to the cell can be made by calling methods below. These include:
        Converting to a type IA model (add transient K current) (species: guineapig-TypeIA).
        Changing "species" to mouse or cat (scales conductances)
        
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
    
        nach : string (default: None)
            nach selects the type of sodium channel that will be used in the model. A channel mechanism
            by that name must exist. The default is 'nacn', from R&M2003.
    
        ttx : Boolean (default: False)
            If ttx is True, then the sodium channel conductance is set to 0 everywhere in the cell.
            Currently, this is not implemented.
    
        species: string (default 'guineapig')
            species defines the channel density that will be inserted for different models. Note that
            if a decorator function is specified, this argument is ignored.
        
        modelType: string (default: None)
            modelType specifies the type of the model that will be used (e.g., "I-c", "I-t").
            modelType is passed to the decorator, or to species_scaling to adjust point models.
        
        debug: boolean (default: False)
            debug is a boolean flag. When set, there will be multiple printouts of progress and parameters.
        
        Returns
        -------
        Nothing
        """
        
        super(TStellateRothman, self).__init__()
        self.i_test_range={'pulse': (-0.15, 0.15, 0.01)}
        if modelType == None:
            modelType = 'I-c'
        if nach == None and species == 'guineapig':
            nach = 'nacn'
        if nach == None and species == 'mouse':
            nach = 'nav11'
            self.i_test_range={'pulse': (-1.0, 1.0, 0.05)}
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'ttx': ttx, 'name': 'TStellate',
                       'morphology': morphology, 'decorator': decorator, 'temperature': None}
        self.vrange = [-70., -55.]
        if morphology is None:
            """
            instantiate a basic soma-only ("point") model
            """
            print( "<< TStellate model: Creating point cell, type={:s} >>".format(modelType))
            soma = h.Section(name="TStellate_Soma_%x" % id(self))  # one compartment of about 29000 um2
            soma.nseg = 1
            self.add_section(soma, 'soma')
        else:
            """
            instantiate a structured model with the morphology as specified by 
            the morphology file
            """
            print "<< TStellate: Creating cell with morphology = %s>>" % morphology
            self.set_morphology(morphology_file=morphology)

        # decorate the morphology with ion channels
        if decorator is None:   # basic model, only on the soma
            self.mechanisms = ['kht', 'ka', 'ihvcn', 'leak', nach]
            for mech in self.mechanisms:
                self.soma.insert(mech)
            self.soma.ek = self.e_k
            self.soma.ena = self.e_na
            self.soma().ihvcn.eh = self.e_h
            self.soma().leak.erev = self.e_leak
            self.species_scaling(silent=True, species=species, modelType=modelType)  # set the default type II cell parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorate()
#        print 'Mechanisms inserted: ', self.mechanisms
        
        self.get_mechs(self.soma)
#        self.cell_initialize(vrange=self.vrange)
        if debug:
                print "<< T-stellate: JSR Stellate Type 1 cell model created >>"

    def species_scaling(self, species='guineapig', modelType='I-c', silent=True):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        Used ONLY for point models.
        
        This scaling routine also sets the temperature for the model to a default value. Some models
        can be run at multiple temperatures, and so a default from one of the temperatures is used.
        The calling cell.set_temperature(newtemp) will change the conductances and reinitialize
        the cell to the new temperature settings.
        
        Parameters
        ----------
        species : string (default: 'guineapig')
            name of the species to use for scaling the conductances in the base point model
            Must be one of mouse, cat, guineapig
        
        modelType: string (default: 'I-c')
            definition of model type from RM03 models, type I-c or type I-t
        
        silent : boolean (default: True)
            run silently (True) or verbosely (False)
        """
        soma = self.soma
        if species == 'mouse' and modelType == 'I-c':
            # use conductance levels from Cao et al.,  J. Neurophys., 2007.
            # model description in Xie and Manis 2013. Note that
            # conductances were not scaled for temperature (rates were)
            # so here we reset the default Q10's for conductance (g) to 1.0
            print '  Setting Conductances for mouse I-c Tstellate cell, Xie and Manis, 2013'
            self.c_m = 0.9  # default in units of F/cm^2
            self.vrange = [-75., -55.]
            self.set_soma_size_from_Cm(25.0)
            self._valid_temperatures = (34.,)
            if self.status['temperature'] is None:
                self.set_temperature(34.)
            #self.adjust_na_chans(soma, gbar=800.)
            self.e_k = -84.
            self.e_na = 50.
            soma().nav11.gbar = nstomho(1800., self.somaarea)
            soma().nav11.vsna = 4.3
            soma.ena = self.e_na
            soma.ek = self.e_k
            soma().kht.gbar = nstomho(250.0, self.somaarea)
            soma().ka.gbar = nstomho(0.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(18.0, self.somaarea)
            soma().ihvcn.eh = -43 # Rodrigues and Oertel, 2006
            soma().leak.gbar = nstomho(8.0, self.somaarea)
            soma().leak.erev = -65.0
            self.axonsf = 0.5
            
        elif species == 'guineapig' and modelType == 'I-c':  # values from R&M 2003, Type I
            print '  Setting Conductances for Guinea Pig I-c, Rothman and Manis, 2003'
            self.c_m = 0.9  # default in units of F/cm^2
            self.vrange = [-75., -55.]
            self.set_soma_size_from_Cm(12.0)
            self._valid_temperatures = (22., 38.)
            if self.status['temperature'] is None:
                self.set_temperature(22.)
            sf = 1.0
            if self.status['temperature'] == 38.:  # adjust for 2003 model conductance levels at 38
                sf = 3.03  # Q10 of 2, 22->38C. (p3106, R&M2003c)
                # note that kinetics are scaled in the mod file.
            self.adjust_na_chans(soma, sf=sf)
            soma().kht.gbar = nstomho(150.0, self.somaarea)*sf
            soma().ka.gbar = nstomho(0.0, self.somaarea)*sf
            soma().ihvcn.gbar = nstomho(0.5, self.somaarea)*sf
            soma().leak.gbar = nstomho(2.0, self.somaarea)*sf
            soma().leak.erev = -65.0
            self.axonsf = 0.5
            
        elif species == 'guineapig' and modelType =='I-t':
            print '  Setting Conductances for Guinea Pig, I-t, Rothman and Manis, 2003'
            # guinea pig data from Rothman and Manis, 2003, type It
            self.c_m = 0.9  # default in units of F/cm^2
            self.set_soma_size_from_Cm(12.0)
            self._valid_temperatures = (22., 38.)
            if self.status['temperature'] is None:
                self.set_temperature(22.)
            sf = 1.0
            if self.status['temperature'] == 38.:  # adjust for 2003 model conductance levels at 38
                sf = 3.03  # Q10 of 2, 22->38C. (p3106, R&M2003c)
                # note that kinetics are scaled in the mod file.
            self.adjust_na_chans(soma, sf)
            soma().kht.gbar = nstomho(80.0, self.somaarea)*sf
            soma().ka.gbar = nstomho(65.0, self.somaarea)*sf
            soma().ihvcn.gbar = nstomho(0.5, self.somaarea)*sf
            soma().leak.gbar = nstomho(2.0, self.somaarea)*sf
            soma().leak.erev = -65.0
            self.axonsf = 0.5

        else:
            raise ValueError('Species %s or species-type %s is not recognized for T-stellate cells' % (species, type))

        self.status['species'] = species
        self.status['modelType'] = modelType
        self.check_temperature()

    def channel_manager(self, modelType='RM03'):
        """
        This routine defines channel density maps and distance map patterns
        for each type of compartment in the cell. The maps
        are used by the ChannelDecorator class (specifically, called from it's private
        _biophys function) to decorate the cell membrane with channels.
        
        Parameters
        ----------
        modelType : string (default: 'RM03')
            A string that defines the type of the model. Currently, 3 types are implemented:
            RM03: Rothman and Manis, 2003 somatic densities for guinea pig
            XM13: Xie and Manis, 2013, somatic densities for mouse
            XM13PasDend: XM13, but with only passive dendrites, no channels.
        
        Returns
        -------
        Nothing
        
        Notes
        -----
        
        This routine defines the following variables for the class:
            
            - conductances (gBar)
            - a channelMap (dictonary of channel densities in defined anatomical compartments)
            - a current injection range for IV's (when testing)
            - a distance map, which defines how selected conductances in selected compartments
                will change with distance. This includes both linear and exponential gradients,
                the minimum conductance at the end of the gradient, and the space constant or
                slope for the gradient.
        
        """
        if modelType == 'RM03':
            totcap = 12.0E-12  # TStellate cell (type I) from Rothman and Manis, 2003, as base model
            refarea = totcap / self.c_m  # see above for units
            # Type I stellate Rothman and Manis, 2003c
            self._valid_temperatures = (22., 38.)
            if self.status['temperature'] is None:
                self.set_temperature(22.)
            sf = 1.0
            if self.status['temperature'] == 38.:  # adjust for 2003 model conductance levels at 38
                sf = 3.03  # Q10 of 2, 22->38C. (p3106, R&M2003c)
            self.gBar = Params(nabar=sf*1000.0E-9/refarea,
                               khtbar=sf*150.0E-9/refarea,
                               kltbar=sf*0.0E-9/refarea,
                               ihbar=sf*0.5E-9/refarea,
                               leakbar=sf*2.0E-9/refarea,
            )
            
            self.channelMap = {
                'axon': {'nacn': 0.0, 'klt': 0., 'kht': self.gBar.khtbar,
                         'ihvcn': 0., 'leak': self.gBar.leakbar / 4.},
                'hillock': {'nacn': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar,
                             'ihvcn': 0., 'leak': self.gBar.leakbar, },
                'initseg': {'nacn': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar,
                            'ihvcn': self.gBar.ihbar / 2.,
                            'leak': self.gBar.leakbar, },
                'soma': {'nacn': self.gBar.nabar, 'klt': self.gBar.kltbar,
                         'kht': self.gBar.khtbar, 'ihvcn': self.gBar.ihbar,
                         'leak': self.gBar.leakbar, },
                'dend': {'nacn': self.gBar.nabar / 2.0, 'klt': 0., 'kht': self.gBar.khtbar * 0.5,
                         'ihvcn': self.gBar.ihbar / 3., 'leak': self.gBar.leakbar * 0.5, },
                'apic': {'nacn': 0.0, 'klt': 0., 'kht': self.gBar.khtbar * 0.2,
                         'ihvcn': self.gBar.ihbar / 4.,
                         'leak': self.gBar.leakbar * 0.2, },
            }
            self.irange = np.linspace(-0.1, 0.1, 7)
            self.distMap = {'dend': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.}}, # linear with distance, gminf (factor) is multiplied by gbar
                            'apic': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.}}, # gradients are: flat, linear, exponential
                            }

        elif modelType  == 'XM13':
            totcap = 25.0E-12  # Base model from Xie and Manis, 2013 for type I stellate cell
            refarea = totcap / self.c_m  # see above for units
            self._valid_temperatures = (34.,)
            if self.status['temperature'] is None:
                self.set_temperature(34.)
            self.gBar = Params(nabar=1800.0E-9/refarea,
                               khtbar=250.0E-9/refarea,
                               kltbar=0.0E-9/refarea,
                               ihbar=18.0E-9/refarea,
                               leakbar=8.0E-9/refarea,
            )
            self.channelMap = {
                'axon': {'nav11': 0.0, 'klt': 0., 'kht': self.gBar.khtbar,
                         'ihvcn': 0., 'leak': self.gBar.leakbar / 4.},
                'hillock': {'nav11': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar,
                            'ihvcn': 0.,
                            'leak': self.gBar.leakbar, },
                'initseg': {'nav11': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar,
                            'ihvcn': self.gBar.ihbar / 2.,
                            'leak': self.gBar.leakbar, },
                'soma': {'nav11': self.gBar.nabar, 'klt': self.gBar.kltbar,
                         'kht': self.gBar.khtbar, 'ihvcn': self.gBar.ihbar,
                         'leak': self.gBar.leakbar, },
                'dend': {'nav11': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar * 0.5,
                         'ihvcn': self.gBar.ihbar / 3., 'leak': self.gBar.leakbar * 0.5, },
                'apic': {'nav11': 0.0, 'klt': 0., 'kht': self.gBar.khtbar * 0.2,
                         'ihvcn': self.gBar.ihbar / 4.,
                         'leak': self.gBar.leakbar * 0.2, },
            }
            self.irange = np.linspace(-0.5, 0.5, 9)
            self.distMap = {'dend': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'nav11': {'gradient': 'exp', 'gminf': 0., 'lambda': 100.}}, # linear with distance, gminf (factor) is multiplied by gbar
                            'apic': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'nav11': {'gradient': 'exp', 'gminf': 0., 'lambda': 100.}}, # gradients are: flat, linear, exponential
                            }

        elif modelType == 'XM13PasDend':
            # bushy form Xie and Manis, 2013, based on Cao and Oertel mouse conductances
            # passive dendritestotcap = 26.0E-12 # uF/cm2 
            totcap = 26.0E-12 # uF/cm2 
            refarea = totcap  / self.c_m  # see above for units
            self._valid_temperatures = (34.,)
            if self.status['temperature'] is None:
                self.set_temperature(34.)
            self.gBar = Params(nabar=1000.0E-9/refarea,
                               khtbar=150.0E-9/refarea,
                               kltbar=0.0E-9/refarea,
                               ihbar=0.5E-9/refarea,
                               leakbar=2.0E-9/refarea,
            )
            self.channelMap = {
                'axon': {'nav11': self.gBar.nabar*0, 'klt': self.gBar.kltbar * 0.25, 'kht': self.gBar.khtbar, 'ihvcn': 0.,
                         'leak': self.gBar.leakbar * 0.25},
                'hillock': {'nav11': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar, 'ihvcn': 0.,
                            'leak': self.gBar.leakbar, },
                'initseg': {'nav11': self.gBar.nabar*3, 'klt': self.gBar.kltbar*2, 'kht': self.gBar.khtbar*2,
                            'ihvcn': self.gBar.ihbar * 0.5, 'leak': self.gBar.leakbar, },
                'soma': {'nav11': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar,
                         'ihvcn': self.gBar.ihbar, 'leak': self.gBar.leakbar, },
                'dend': {'nav11': self.gBar.nabar * 0.0, 'klt': self.gBar.kltbar*0 , 'kht': self.gBar.khtbar*0,
                         'ihvcn': self.gBar.ihbar*0, 'leak': self.gBar.leakbar*0.5, },
                'apic': {'nav11': self.gBar.nabar * 0.0, 'klt': self.gBar.kltbar * 0, 'kht': self.gBar.khtbar * 0.,
                         'ihvcn': self.gBar.ihbar *0., 'leak': self.gBar.leakbar * 0.25, },
            }
            self.irange = np.linspace(-1, 1, 21)
            self.distMap = {'dend': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 200.},
                                     'kht': {'gradient': 'llinear', 'gminf': 0., 'lambda': 200.},
                                     'nav11': {'gradient': 'linear', 'gminf': 0., 'lambda': 200.}}, # linear with distance, gminf (factor) is multiplied by gbar
                            'apic': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'nav11': {'gradient': 'exp', 'gminf': 0., 'lambda': 200.}}, # gradients are: flat, linear, exponential
                            }
        else:
            raise ValueError('model type %s is not implemented' % modelType)
        self.check_temperature()
        
    def adjust_na_chans(self, soma, sf=1.0, gbar=1000., debug=False):
        """
        Adjust the sodium channel conductance, depending on the type of conductance
        
        Parameters
        ----------
        soma : NEURON section object (required)
            This identifies the soma object whose sodium channel complement will have it's
            conductances adjusted depending on the sodium channel type
        gbar : float (default: 1000.)
            The "maximal" conductance to be set in the model.
        debug : boolean (default: False)
            A flag the prints out messages to confirm the operations applied.
            
        Returns
        -------
        Nothing
        """
        if self.status['ttx']:
            gnabar = 0.0
        else:
            gnabar = nstomho(gbar, self.somaarea)*sf
        nach = self.status['na']
        if nach == 'jsrna':
            soma().jsrna.gbar = gnabar*sf
            soma.ena = self.e_na
            if debug:
                print 'jsrna gbar: ', soma().jsrna.gbar
        elif nach == 'nav11':
            soma().nav11.gbar = gnabar
            soma.ena = self.e_na
            soma().nav11.vsna = 4.3
            if debug:
                print "tstellate using inva11"
            print 'nav11 gbar: ', soma().nav11.gbar
            print 'nav11 vsna: ', soma().nav11.vsna
        elif nach == 'na':
            soma().na.gbar = gnabar
            soma.ena = self.e_na
            if debug:
                print 'na gbar: ', soma().na.gbar
        elif  nach == 'nacn':
            soma().nacn.gbar = gnabar
            soma.ena = self.e_na
            if debug:
                print 'nacn gbar: ', soma().nacn.gbar
        else:
            raise ValueError("tstellate setting Na channels: channel %s not known" % nach)

    def add_axon(self):
        Cell.add_axon(self, self.soma, self.somaarea, self.c_m, self.R_a, self.axonsf)

    def add_dendrites(self):
        """
        Add simple unbranched dendrites to basic Rothman Type I models.
        The dendrites have some kht and ih current
        """
        cs = False  # not implemented outside here - internal Cesium.
        nDend = range(4) # these will be simple, unbranced, N=4 dendrites
        dendrites=[]
        for i in nDend:
            dendrites.append(h.Section(cell=self.soma))
        for i in nDend:
            dendrites[i].connect(self.soma)
            dendrites[i].L = 200 # length of the dendrite (not tapered)
            dendrites[i].diam = 1.5 # dendrite diameter
            dendrites[i].nseg = 21 # # segments in dendrites
            dendrites[i].Ra = 150 # ohm.cm
            dendrites[i].insert('kht')
            if cs is False:
                dendrites[i]().kht.gbar = 0.005 # a little Ht
            else:
                dendrites[i]().kht.gbar = 0.0
            dendrites[i].insert('leak') # leak
            dendrites[i]().leak.gbar = 0.0001
            dendrites[i].insert('ihvcn') # some H current
            dendrites[i]().ihvcn.gbar = 0.# 0.001
            dendrites[i]().ihvcn.eh = -43.0
        self.maindend = dendrites
        self.status['dendrites'] = True
        self.add_section(self.maindend, 'maindend')


class TStellateNav11(TStellate):
    """
    VCN T-stellate cell (Mouse) from Xie and Manis, 2013. 
    Using nav11 sodium channel model.
    
    """
    def __init__(self, morphology=None, decorator=None, nach='nav11', ttx=False,
                species='mouse', modelType=None, debug=False):
        """
        Initialize a planar stellate (T-stellate) cell as a point model, using the default parameters for
        mouse from Xie and Manis, 2013. 
        Modifications to the cell can be made by calling methods below.
        Changing "species": This routine only supports "mouse"
        Note: in the original model, the temperature scaling applied only to the rate constants, and not
                to the conductance. Therefore, the conductances here need to be adjusted to compensate for the
                way the mechanisms are currently implemented (so that they scale correctly to the values
                used in Xie and Manis, 2013). This is done by setting q10g (the q10 for conductances) to 1
                before setting up the rest of the model parameters. For those conducantances in which a Q10 for 
                conductance is implemented, the value is typically 2.
                
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

        nach : string (default: 'nav11')
            nach selects the type of sodium channel that will be used in the model. A channel mechanims
            by that name must exist. 

        ttx : Boolean (default: False)
            If ttx is True, then the sodium channel conductance is set to 0 everywhere in the cell.
            Currently, this is not implemented.

        species: string (default 'mouse')
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
    
        super(TStellateNav11, self).__init__()
        if modelType == None:
            modelType = 'XM13'
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'ttx': ttx, 'name': 'TStellate',
                       'morphology': morphology, 'decorator': decorator}

        self.i_test_range=(-1, 1., 0.05)
    
        if morphology is None:
            """
            instantiate a basic soma-only ("point") model
            """
            print( "<< TStellate Xie&Manis 2013 model: Creating point cell, type={:s} >>".format(modelType))
            soma = h.Section(name="TStellate_Soma_%x" % id(self))  # one compartment of about 29000 um2
            soma.nseg = 1
            self.add_section(soma, 'soma')
        else:
            """
            instantiate a structured model with the morphology as specified by 
            the morphology file
            """
            print "<< TStellate Xie&Manis 2013 model: Creating structured cell >>"
            self.set_morphology(morphology_file=morphology)

        # decorate the morphology with ion channels
        if decorator is None:   # basic model, only on the soma
            self.mechanisms = ['kht', 'ka', 'ihvcn', 'leak', nach]
            for mech in self.mechanisms:
                self.soma.insert(mech)
            self.soma.ek = self.e_k
            self.soma.ena = self.e_na
            self.soma().ihvcn.eh = self.e_h
            self.soma().leak.erev = self.e_leak
            self.soma().cm = 1.0
            self.species_scaling(silent=True, species=species, modelType=modelType)  # set the default type II cell parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorate()
    
        self.get_mechs(self.soma)
        self.cell_initialize(vrange=self.vrange)
#        self.print_mechs(self.soma)
        if debug:
                print "<< T-stellate: Xie&Manis 2013 cell model created >>"

    def species_scaling(self, species='mouse', modelType='I-c', silent=True):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        Used ONLY for point models.
                
        Parameters
        ----------
        species : string (default: 'guineapig')
            name of the species to use for scaling the conductances in the base point model
            Must be one of mouse, cat, guineapig
        
        modelType: string (default: 'I-c')
            definition of model type from RM03 models, type I-c or type I-t
        
        silent : boolean (default: True)
            run silently (True) or verbosely (False)
        """
        soma = self.soma
        if species == 'mouse' and modelType == 'XM13':
            # use conductance levels from Cao et al.,  J. Neurophys., 2007.
            # original temp for model: 32 C
            print 'Mouse Tstellate cell, Xie and Manis, 2013'
            self.set_soma_size_from_Cm(25.0)
            self.adjust_na_chans(soma, gbar=800.)  # inav11 does not scale conductance
            self.e_k = -84.
            self.e_na = 50.
            soma.ek = self.e_k
            soma.ena = self.e_na
            soma().kht.gbar = nstomho(250.0, self.somaarea)
            soma().ka.gbar = nstomho(0.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(18.0, self.somaarea)
            soma().ihvcn.eh = -43 # Rodrigues and Oertel, 2006
            soma().leak.gbar = nstomho(8.0, self.somaarea)
            soma().leak.erev = -65.0

        else:
            raise ValueError('Species %s or species-type %s is not recognized for T-stellate XM13 cells' % (species, type))

        self.status['species'] = species
        self.status['modelType'] = modelType
        # self.cell_initialize(showinfo=False)
        # if not silent:
        #     print 'set cell as: ', species
        #     print ' with Vm rest = %f' % self.vm0

    def channel_manager(self, modelType='XM13'):
        """
        This routine defines channel density maps and distance map patterns
        for each type of compartment in the cell. The maps
        are used by the ChannelDecorator class (specifically, called from it's private
        _biophys function) to decorate the cell membrane with channels.
        
        Parameters
        ----------
        modelType : string (default: 'XM13')
            A string that defines the type of the model. Currently, 3 types are implemented:
            XM13: Xie and Manis, 2013, somatic densities for mouse
            XM13PasDend: XM13, but with only passive dendrites, no channels.
        
        Returns
        -------
        Nothing
        
        Notes
        -----
        
        This routine defines the following variables for the class:
            
            - conductances (gBar)
            - a channelMap (dictonary of channel densities in defined anatomical compartments)
            - a current injection range for IV's (when testing)
            - a distance map, which defines how selected conductances in selected compartments
                will change with distance. This includes both linear and exponential gradients,
                the minimum conductance at the end of the gradient, and the space constant or
                slope for the gradient.
        
        """
        if modelType  == 'XM13':
            totcap = 25.0E-12  # Base model from Xie and Manis, 2013 for type I stellate cell
            refarea = totcap / self.c_m  # see above for units
            self.gBar = Params(nabar=800.0E-9/refarea,
                               khtbar=250.0E-9/refarea,
                               kltbar=0.0E-9/refarea,
                               ihbar=18.0E-9/refarea,
                               leakbar=8.0E-9/refarea,
            )
            self.channelMap = {
                'axon': {'nav11': 0.0, 'klt': 0., 'kht': self.gBar.khtbar,
                         'ihvcn': 0., 'leak': self.gBar.leakbar / 4.},
                'hillock': {'nav11': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar,
                            'ihvcn': 0.,
                            'leak': self.gBar.leakbar, },
                'initseg': {'nav11': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar,
                            'ihvcn': self.gBar.ihbar / 2.,
                            'leak': self.gBar.leakbar, },
                'soma': {'nav11': self.gBar.nabar, 'klt': self.gBar.kltbar,
                         'kht': self.gBar.khtbar, 'ihvcn': self.gBar.ihbar,
                         'leak': self.gBar.leakbar, },
                'dend': {'nav11': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar * 0.5,
                         'ihvcn': self.gBar.ihbar / 3., 'leak': self.gBar.leakbar * 0.5, },
            }
            self.irange = np.linspace(-1.0, 1.0, 21)
            self.distMap = {'dend': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'nav11': {'gradient': 'exp', 'gminf': 0., 'lambda': 100.}}, # linear with distance, gminf (factor) is multiplied by gbar
                            }

        elif modelType == 'XM13PasDend':
            # bushy form Xie and Manis, 2013, based on Cao and Oertel mouse conductances
            # passive dendrites
            totcap = 26.0E-12 # uF/cm2 
            refarea = totcap  / self.c_m  # see above for units
            self.gBar = Params(nabar=1000.0E-9/refarea,
                               khtbar=150.0E-9/refarea,
                               kltbar=0.0E-9/refarea,
                               ihbar=0.5E-9/refarea,
                               leakbar=2.0E-9/refarea,
            )
            self.channelMap = {
                'axon': {'nav11': self.gBar.nabar*0, 'klt': self.gBar.kltbar * 0.25, 'kht': self.gBar.khtbar, 'ihvcn': 0.,
                         'leak': self.gBar.leakbar * 0.25},
                'hillock': {'nav11': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar, 'ihvcn': 0.,
                            'leak': self.gBar.leakbar, },
                'initseg': {'nav11': self.gBar.nabar*3, 'klt': self.gBar.kltbar*2, 'kht': self.gBar.khtbar*2,
                            'ihvcn': self.gBar.ihbar * 0.5, 'leak': self.gBar.leakbar, },
                'soma': {'nav11': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar,
                         'ihvcn': self.gBar.ihbar, 'leak': self.gBar.leakbar, },
                'dend': {'nav11': self.gBar.nabar * 0.0, 'klt': self.gBar.kltbar*0 , 'kht': self.gBar.khtbar*0,
                         'ihvcn': self.gBar.ihbar*0, 'leak': self.gBar.leakbar*0.5, },
            }
            self.irange = np.linspace(-1, 1, 21)
            self.distMap = {'dend': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 200.},
                                     'kht': {'gradient': 'llinear', 'gminf': 0., 'lambda': 200.},
                                     'nav11': {'gradient': 'linear', 'gminf': 0., 'lambda': 200.}}, # linear with distance, gminf (factor) is multiplied by gbar
                            }
        else:
            raise ValueError('model type %s is not implemented' % modelType)

    def adjust_na_chans(self, soma, gbar=800., debug=False):
        """
        Adjust the sodium channel conductance, depending on the type of conductance
        
        Parameters
        ----------
        soma : NEURON section object (required)
            This identifies the soma object whose sodium channel complement will have it's
            conductances adjusted depending on the sodium channel type
        gbar : float (default: 800.)
            The "maximal" conductance to be set in the model.
        debug : boolean (default: False)
            A flag the prints out messages to confirm the operations applied.
            
        Returns
        -------
        Nothing
        """
        if self.status['ttx']:
            gnabar = 0.0
        else:
            gnabar = nstomho(gbar, self.somaarea)
        nach = self.status['na']
        if nach == 'nav11':
            soma().nav11.gbar = gnabar
            soma.ena = self.e_na
            soma().nav11.vsna = 4.3
            if debug:
                print "tstellate using inva11"
        else:
            raise ValueError("tstellate setting Na channels only supporting nav11: channel %s not known" % nach)

