from neuron import h
import numpy as np
import neuron as nrn

from .cell import Cell
#from .. import synapses
from ..util import nstomho
from ..util import Params
#from .. import data

__all__ = ['TStellate', 'TStellateRothman', 'TStellateNav11', 'TStellateFast'] 


class TStellate(Cell):
    
    type = 'tstellate'

    @classmethod
    def create(cls, model='RM03', **kwds):
        if model == 'RM03':
            return TStellateRothman(**kwds)
        # elif model == 'Nav11':   # these models are not supported.
        #     return TStellateNav11(**kwds)
        # elif model == 'fast':
        #     return TStellateFast(**kwds)
        else:
            raise ValueError ('TStellate type %s is unknown', type)

    def make_psd(self, terminal, psd_type, **kwds):
        pre_sec = terminal.section
        pre_cell = terminal.cell
        post_sec = self.soma

        if psd_type == 'simple':
            return self.make_exp2_psd(post_sec, terminal)
        
        elif psd_type == 'multisite':
            if pre_cell.type == 'sgc':
                # Max conductances for the glu mechanisms are calibrated by 
                # running `synapses/tests/test_psd.py`. The test should fail
                # if these values are incorrect:
                AMPA_gmax = 0.22479596944138733*1e3  # factor of 1e3 scales to pS (.mod mechanisms) from nS.
                NMDA_gmax = 0.12281291946623739*1e3
                return self.make_glu_psd(post_sec, terminal, AMPA_gmax, NMDA_gmax)
            elif pre_cell.type == 'dstellate':
                # Get GLY kinetic constants from database 
                return self.make_gly_psd(post_sec, terminal, type='glyfast')
            else:
                raise TypeError("Cannot make PSD for %s => %s" % 
                                (pre_cell.type, self.type))
        else:
            raise ValueError("Unsupported psd type %s" % psd_type)


class TStellateRothman(TStellate):
    """
    VCN T-stellate base model.
    Rothman and Manis, 2003abc (Type I-c, Type I-t)
    """
    def __init__(self, morphology=None, decorator=None, nach='na', ttx=False,
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
        
        super(TStellateRothman, self).__init__()
        if modelType == None:
            modelType = 'I-c'
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'ttx': ttx, 'name': 'TStellate',
                       'morphology': morphology, 'decorator': decorator}

        self.i_test_range=(-0.15, 0.15, 0.01)
        
        if morphology is None:
            """
            instantiate a basic soma-only ("point") model
            """
            print "<< TStellate model: Creating point cell using JSR parameters >>"
            soma = h.Section(name="TStellate_Soma_%x" % id(self))  # one compartment of about 29000 um2
            soma.nseg = 1
            self.add_section(soma, 'soma')
        else:
            """
            instantiate a structured model with the morphology as specified by 
            the morphology file
            """
            print "<< TStellate model: Creating structured cell using JSR parameters >>"
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
        self.cell_initialize()
        if debug:
                print "<< T-stellate: JSR Stellate Type 1 cell model created >>"

    def species_scaling(self, species='guineapig', modelType='I-c', silent=True):
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
        if species == 'mouse' and modelType == 'I-c':
            # use conductance levels from Cao et al.,  J. Neurophys., 2007.
            #print 'Mouse Tstellate cell'
            self.set_soma_size_from_Cm(25.0)
            self.adjust_na_chans(soma, gbar=800.)
            soma().kht.gbar = nstomho(250.0, self.somaarea)
            soma().ka.gbar = nstomho(0.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(18.0, self.somaarea)
            soma().ihvcn.eh = -43 # Rodrigues and Oertel, 2006
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            soma().leak.erev = -65.0
            self.axonsf = 0.5
        elif species == 'guineapig' and modelType == 'I-c':  # values from R&M 2003, Type I
            self.set_soma_size_from_Cm(12.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().ka.gbar = nstomho(0.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(0.5, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            soma().leak.erev = -65.0
            self.axonsf = 0.5
        elif species == 'guineapig' and modelType =='I-t':
            # guinea pig data from Rothman and Manis, 2003, type It
            self.set_soma_size_from_Cm(12.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(80.0, self.somaarea)
            soma().ka.gbar = nstomho(65.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(0.5, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            soma().leak.erev = -65.0
            self.axonsf = 0.5
        elif species == 'cat' and modelType == 'I-c':  # a cat is a big guinea pig Type I
            self.set_soma_size_from_Cm(30.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().ka.gbar = nstomho(0.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(0.5, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            soma().leak.erev = -65.0
            self.axonsf = 1.0
        else:
            raise ValueError('Species %s or species-type %s is not recognized for T-stellate cells' % (species, type))

        self.status['species'] = species
        self.status['modelType'] = modelType
        # self.cell_initialize(showinfo=False)
        # if not silent:
        #     print 'set cell as: ', species
        #     print ' with Vm rest = %f' % self.vm0

    def channel_manager(self, modelType='RM03'):
        """
        This routine defines channel density maps and distance map patterns
        for each type of compartment in the cell. The maps
        are used by the ChannelDecorator class to(specifically, it's private
        _biophys function) to decorate the cell membrane.
        
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
            self.gBar = Params(nabar=1000.0E-9/refarea,
                               khtbar=150.0E-9/refarea,
                               kltbar=0.0E-9/refarea,
                               ihbar=0.5E-9/refarea,
                               leakbar=2.0E-9/refarea,
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
        
    def adjust_na_chans(self, soma, gbar=1000., debug=False):
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
                print "tstellate using inva11"
            print 'nav11 gbar: ', soma().nav11.gbar
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
    VCN T-stellate cell setup from Rothman and Manis, 2003, 
    using nav11 sodium channel model
    
    ttx: if True, turns off sodium channels
    cs: if True, turns off K channels (e.g., cesium in pipette). 
    dend: if True, adds dendrites to the model, based roughly on White et al.,
    1994)
    NOTE: This has been modified from it's original from
    for use in simulating MOUSE stellate cells.
    """
    def __init__(self, debug=False, ttx=False, cs = False, message=None, dend=False):
        super(TStellateNav11, self).__init__()
        print ("T-STELLATE NAV11",
            "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        soma = h.Section() # one compartment of about 29000 um2
        v_potassium = -80       # potassium reversal potential
        v_sodium = 50           # sodium reversal potential

        cm = 1.0
        scalefactor = 1.0 # This determines the relative size of the cell
        rinsf = 1.0           # input resistance adjustment (also current...)
        totcap = scalefactor * 25.0 # cap in pF for cell
        effcap = totcap # sometimes we change capacitance - that's effcap
        somaarea = totcap * 1E-6 / cm # pf -> uF, cm = 1uf/cm^2 nominal;
        lstd = 1E4 * ((somaarea / 3.14159) ** 0.5) # convert from cm to um

        soma.nseg = 1
        soma.diam = lstd
        soma.L = lstd

        dendrites=[]
        if dend is True:
            nDend = range(4) # these will be simple, unbranced, N=4 dendrites

            
    #    print nnodes
            for i in nDend:
                dendrites.append(h.Section(cell=soma))
            for i in nDend:
                dendrites[i].connect(soma)
                dendrites[i].L = 200 # length of the dendrite (not tapered)
                dendrites[i].diam = 1.5 # dendritic diameter
                dendrites[i].nseg = 21 # # segments in dendrites
                dendrites[i].Ra = 150 # ohm.cm
                d = dendrites[i]
                ds = d()
                d.insert('kht')
                if cs is False:
                    ds.kht.gbar = 0.005 # a little Ht
                else:
                    ds.kht.gbar = 0.0
                d.insert('leak') # leak
                ds.leak.gbar = 0.0001
                d.insert('ihvcn') # some H current
                ds.ihvcn.gbar = 0.# 0.001
                ds.ihvcn.eh = -43.0
        seg = soma
        seg.insert('kht')
        seg.insert('nav11')
        seg.insert('ka')
        seg.insert('ihvcn')
        seg.insert('leak')
        seg.ena = 10
        seg.ek = -84
        s = soma()
        if ttx is False:
            s.nav11.gbar = nstomho(1800.0, somaarea) * scalefactor
        else:
            s.nav11.gbar = 0.0 # print s.nav11.gnat
        s.nav11.vsna = 4.3 # was 8

        if cs is False:
            s.kht.gbar = nstomho(200.0, somaarea) * scalefactor
            s.ka.gbar = nstomho(0.0, somaarea) * scalefactor
        else:
            s.kht.gbar = 0.
            s.ka.gbar = 0.
        s.ihvcn.gbar = nstomho(18.0, somaarea) * scalefactor # was 10
        s.ihvcn.eh = -43 # Rodrigues and Oertel, 2006
    #    print 'ih vcn vh: %f ' % (s.ihvcn.vh)
        s.leak.gbar = nstomho(2.0, somaarea) * scalefactor
        vm0 = -63.9
        if debug:
            if message is None:
                print ("<< T-stellate: JSR Stellate Type 1 cell model created",
                " - modified for mouse >>")
            else:
                print message
    # print dendrites
        self.add_section(soma, 'soma')
        self.add_section(dendrites, 'dendrite')

class TStellateFast(TStellate):
    """ 
    VCN t-stellate model based on Rothman and Manis 2003, but with fast sodium
    channel 
    """
    def __init__(self, debug=False, ttx=False, message=None, dend=False):
        super(TStellateFast, self).__init__()
        soma = h.Section() # one compartment of about 29000 um2
        v_potassium = -80       # potassium reversal potential
        v_sodium = 50           # sodium reversal potential

        cm = 1.0
        scalefactor = 1.0 # This determines the relative size of the cell
        rinsf = 1.0           # input resistance adjustment (also current...)
        totcap = scalefactor * 20.0 # cap in pF for cell
        effcap = totcap # sometimes we change capacitance - that's effcap
        somaarea = totcap * 1E-6 / cm # pf -> uF, cm = 1uf/cm^2 nominal;
        lstd = 1E4 * ((somaarea / 3.14159) ** 0.5) # convert from cm to um

        soma.nseg = 1
        soma.diam = lstd
        soma.L = lstd

        seg = soma
        seg.insert('kht')
        seg.insert('nav11')
        seg.insert('ka')
        # 'it' is not part of canonical model; 
        # just trying it to reproduce some data.
        #seg.insert('it') # low-voltage activated ca channel
        seg.insert('ihvcn')
        #seg.insert('iH_std')
        seg.insert('leak')
        seg.ena = 10
        seg.ek = -80
        #seg.eh = -40 # Rodrigues and Oertel, 2006

        s = soma()
        if ttx is False:
            s.nav11.gbar = nstomho(1500.0, somaarea) * scalefactor
        else:
            s.nav11.gbar = 0.0 # print s.nav11.gnat
        s.nav11.vsna = 4.3 # was 8
        s.kht.gbar = nstomho(380.0 * 2.0, somaarea) * scalefactor
        s.ka.gbar = nstomho(90.0, somaarea) * scalefactor # was 280
        
        #s.it.gbar = nstomho(14.0 * 4.0, somaarea) * scalefactor
        # this creates a better rebound! with it
        #s.it.vshift = -16
        
        #  was 16.5to allow the tau shift to be about right so it is not so fast.
        #s.iH_std.gbar = nstomho(100.0, somaarea) * scalefactor
        #s.iH_std.vshift = 1.8
        s.ihvcn.gbar = nstomho(220.0, somaarea) * scalefactor
        s.ihvcn.vshift = 16.0
        s.ihvcn.eh = -43
        s.leak.gbar = nstomho(18.0, somaarea) * scalefactor
        s.leak.e = -61
        vm0 = -60.0
        if debug:
            if message is None:
                print "<< T-stellate: JSR Stellate Type 1 cell model created >>"
            else:
                print message
        
        self.add_section(soma, 'soma')

