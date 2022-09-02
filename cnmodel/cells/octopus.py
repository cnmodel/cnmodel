from __future__ import print_function
from neuron import h
from ..util import nstomho
from ..util import Params
import numpy as np
from .cell import Cell
from .. import data
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

__all__ = ['Octopus', 'OctopusRothman', 'OctopusSpencer']

class Octopus(Cell):

    celltype = 'octopus'
    scaled = False
    
    @classmethod
    def create(cls, modelType='RM03', **kwds):
        if modelType in ['RM03', 'II-o']:
            return OctopusRothman(**kwds)
        elif modelType == 'Spencer':
            return OctopusSpencer(**kwds)
        else:
            raise ValueError ('Octopus cell type %s is unknown' % modelType)

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
            if terminal.cell.celltype in ['sgc']:
                weight = data.get('%s_synapse' % terminal.cell.celltype, species=self.species,
                        post_type=self.celltype, field='weight')
                tau1 = data.get('%s_synapse' % terminal.cell.celltype, species=self.species,
                        post_type=self.celltype, field='tau1')
                tau2 = data.get('%s_synapse' % terminal.cell.celltype, species=self.species,
                        post_type=self.celltype, field='tau2')
                erev = data.get('%s_synapse' % terminal.cell.celltype, species=self.species,
                        post_type=self.celltype, field='erev')
                return self.make_exp2_psd(post_sec, terminal, weight=weight, loc=loc,
                        tau1=tau1, tau2=tau2, erev=erev)
            else:
                raise TypeError("Cannot make simple PSD for %s => %s" % 
                            (terminal.cell.celltype, self.celltype))

            
        elif psd_type == 'multisite':
            if terminal.cell.celltype == 'sgc':
                # Max conductances for the glu mechanisms are calibrated by 
                # running `synapses/tests/test_psd.py`. The test should fail
                # if these values are incorrect
                self.AMPAR_gmax = data.get('sgc_synapse', species=self.species,
                        post_type=self.celltype, field='AMPAR_gmax')*1e3
                self.NMDAR_gmax = data.get('sgc_synapse', species=self.species,
                        post_type=self.celltype, field='NMDAR_gmax')*1e3
                self.Pr = data.get('sgc_synapse', species=self.species,
                        post_type=self.celltype, field='Pr')
                # adjust gmax to correct for initial Pr
                self.AMPAR_gmax = self.AMPAR_gmax/self.Pr
                self.NMDAR_gmax = self.NMDAR_gmax/self.Pr
                # AMPA_gmax = 3.314707700918133*1e3  # factor of 1e3 scales to pS (.mod mechanisms) from nS.
                # NMDA_gmax = 0.4531929783503451*1e3
                if 'AMPAScale' in kwds:
                    self.AMPAR_gmax = self.AMPAR_gmax * kwds['AMPAScale']  # allow scaling of AMPA conductances
                if 'NMDAScale' in kwds:
                    self.NMDAR_gmax = self.NMDAR_gmax*kwds['NMDAScale']
                return self.make_glu_psd(post_sec, terminal, self.AMPAR_gmax, self.NMDAR_gmax, loc=loc)
            elif terminal.cell.celltype == 'dstellate':
                return self.make_gly_psd(post_sec, terminal, psdtype='glyslow', loc=loc)
            else:
                raise TypeError("Cannot make PSD for %s => %s" % 
                            (terminal.cell.celltype, self.celltype))
        else:
            raise ValueError("Unsupported psd type %s" % psd_type)


class OctopusRothman(Octopus, Cell):
    """
    VCN octopus cell model (point cell).
    Rothman and Manis, 2003abc (Type II, with high gklt and hcno - octopus cell h current).
    """

    def __init__(self, morphology=None, decorator=None, nach=None,
                ttx=False, temperature=None,
                species='guineapig', modelType=None, modelName=None,
                debug=False):
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
        
        nach : string (default: None)
            nach selects the type of sodium channel that will be used in the model. A channel mechanism
            by that name must exist. None implies the default channel (jsrna for this model).
        
        ttx : Boolean (default: False)
            If ttx is True, then the sodium channel conductance is set to 0 everywhere in the cell.
            Currently, this is not implemented.
        
        species: string (default 'guineapig')
            species defines the channel density that will be inserted for different models. Note that
            if a decorator function is specified, this argument is ignored.

        modelName: string (default: None)
            modelName specifies the source conductance pattern (RM03, XM13, etc).
            modelName is passed to the decorator, or to species_scaling to adjust point (single cylinder) models.
            
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
        if modelType is None:
            modelType = 'II-o'
        if species == 'guineapig':
            modelName = 'RM03'
            dataset = 'RM03_channels'
            temp = 22.
            nach = 'jsrna'
        elif species == 'mouse':
            if modelName is None:
                modelName = 'XM13'
            if modelName == 'XM13':
                dataset = 'XM13_channels'
                temp = 34.0
            elif modelName  == 'XM13nacncoop':
                dataset = 'XM13_channels_nacncoop'
            else:
                raise ValueError(f"ModelName {self.status['modelName']:s} not recognized for mouse T-stellate cells")
        else:
            raise ValueError(f"Species {species:s} not recognized for {self.celltype:s} cells")

        self.debug = debug
        self.status = {'species': species, 'cellClass': self.celltype, 'modelType': modelType, 'modelName': modelName,
                       self.somaname: True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'ttx': ttx, 'name': self.celltype,
                       'morphology': morphology, 'decorator': decorator, 'temperature': None}
        self.spike_threshold = -50.
        
        soma = self.do_morphology(morphology)

        self.pars = self.get_cellpars(dataset, species=species, modelType=modelType)
        self.status['na'] = self.pars.natype

        # decorate the morphology with ion channels
        if decorator is None:   # basic model, only on the soma
            self.e_leak = -73. # from McGinley et al., 2016
            self.e_h = -38. # from McGinley et al. 
            self.R_a = 195  # McGinley et al. 
            if self.status['species'] == 'mouse':
                # self.mechanisms = ['klt', 'kht', 'hcnobo', 'leak', self.pars.natype]
                self.mechanisms = ['klt', 'kht', 'ihvcn', 'leak', self.pars.natype]
            elif self.status['species'] == 'guineapig':
                self.mechanisms = ['klt', 'kht', 'ihvcn', 'leak', self.pars.natype]
            for mech in self.mechanisms:
                self.soma.insert(mech)
            self.soma.ek = self.e_k
            self.soma.ena = self.e_na
            # if self.status['species'] == 'mouse':
            #     self.soma().hcnobo.eh = self.e_h
            # else:
            self.soma().ihvcn.eh = self.e_h
            self.soma().leak.erev = self.e_leak
            self.soma.Ra = self.R_a
            self.species_scaling(silent=True)  # set the default type II cell parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorate()
        self.save_all_mechs()  # save all mechanisms inserted, location and gbar values...
        self.get_mechs(self.soma)
        
        if self.debug:
                print("<< T-stellate: JSR Stellate Type 1 cell model created >>")

    def get_cellpars(self, dataset, species='guineapig', modelType='I-c'):
        """
        Retrieve parameters for the specifed model type and species from the data tables
        
        dataset : str (no default)
            name of the data table to use
        
        species : str (default: 'guineapig')
            Species table to use
        
        modelType : str (default: 'I-c')
            Model type to get parameters from the table.
        
        """
        cellcap = data.get(dataset, species=species, model_type=modelType,
            field='soma_Cap')
        chtype = data.get(dataset, species=species, model_type=modelType,
            field='na_type')
        pars = Params(cap=cellcap, natype=chtype)

        if self.status['modelName'] == 'RM03':
            for g in ['%s_gbar' % pars.natype, 'kht_gbar', 'ka_gbar', 'ih_gbar', 'leak_gbar', 'leak_erev', 'ih_eh', 'e_k', 'e_na']:
                pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
                    field=g))
        elif self.status['modelName'] == 'XM13':
            for g in ['%s_gbar' % pars.natype,'klt_gbar', 'kht_gbar', 'ka_gbar', 'ihvcn_gbar', 'leak_gbar', 'leak_erev', 'ih_eh', 'e_k', 'e_na']:
                pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
                    field=g))
        # elif self.status['modelName'] == 'mGBC':
        #     for g in ['%s_gbar' % pars.natype, 'kht_gbar', 'ka_gbar', 'ihvcn_gbar', 'leak_gbar', 'leak_erev', 'ih_eh', 'e_k', 'e_na']:
        #         pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
        #             field=g))
        else:
            raise ValueError(f"get_cellpars: Model name {self.status['modelName']} is not yet implemented for cell type {self.celltype.title():s}")
            
        if self.debug:
            pars.show()
        return pars
    
    def species_scaling(self, silent=True):
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
        assert self.scaled is False  # block double scaling!
        self.scaled = True

        soma = self.soma

        if self.status['species'] == 'guineapig' and self.status['modelType'] =='II-o':
            self.c_m = 0.9
            self.i_test_range = {'pulse': (-4.0, 4.0, 0.2)}
            self.vrange = [-70., -57.]  # set a default vrange for searching for rmp
            self.set_soma_size_from_Cm(25.0)
            self._valid_temperatures = (22., 38.)
            if self.status['temperature'] is None:
                self.set_temperature(22.)
            sf = 1.0
            if self.status['temperature'] == 38.:  # adjust for 2003 model conductance levels at 38
                sf = 3.03  # Q10 of 2, 22->38C. (p3106, R&M2003c)
                # note that kinetics are scaled in the mod file.
            self.adjust_na_chans(soma, sf=sf)
            soma().kht.gbar = sf*nstomho(150.0, self.somaarea)  # 6.1 mmho/cm2
            soma().klt.gbar = sf*nstomho(1000.0, self.somaarea)  #  40.7 mmho/cm2  3195?
            soma().ihvcn.gbar = sf*nstomho(30.0, self.somaarea)  # 7.6 mmho/cm2, cf. Bal and Oertel, Spencer et al. 25 u dia cell 40ns?
            soma().leak.gbar = sf*nstomho(2.0, self.somaarea)
            self.axonsf = 1.0
        elif self.status['species'] == 'mouse' and self.status['modelType'] =='II-o':
            self.i_test_range = {'pulse': (-4.0, 4.0, 0.2)}
            self.vrange = [-70., -57.]  # set a default vrange for searching for rmp
            self.set_soma_size_from_Cm(self.pars.cap)
            self._valid_temperatures = (34., )
            if self.status['temperature'] is None:
                self.set_temperature(34.)
            self.adjust_na_chans(soma, sf=1.0)
            soma().kht.gbar = nstomho(self.pars.kht_gbar, self.somaarea)  # 6.1 mmho/cm2
            soma().klt.gbar = nstomho(self.pars.klt_gbar, self.somaarea)  #  40.7 mmho/cm2
            # soma().hcnobo.gbar = nstomho(40.0, self.somaarea)  # 7.6 mmho/cm2, cf. Bal and Oertel, Spencer et al. 25 u dia cell
            soma().ihvcn.gbar = nstomho(self.pars.ihvcn_gbar, self.somaarea)  # 7.6 mmho/cm2, cf. Bal and Oertel, Spencer et al. 25 u dia cell 40ns?
            soma().leak.gbar = nstomho(self.pars.leak_gbar, self.somaarea)
            self.axonsf = 1.0
        else:
            raise ValueError('Species "%s" or species-type "%s" is not recognized for octopus cells' %  (species, type))
        self.check_temperature()


"""
*****************************************************************************************************************************

****************The following is experimental code and has not been fully tested or integrated ******************************

In other words, use at your own risk.
*****************************************************************************************************************************
"""

class OctopusSpencer(Octopus, Cell):
    """
    VCN octopus cell model (with dendrites).
    Based on Spencer et al Front. Comput. Neurosci., 22 October 2012
    https://doi.org/10.3389/fncom.2012.00083
    """

    def __init__(self, morphology=None, decorator=None, nach='jsrna', ttx=False,
                species='guineapig', modelType=None, debug=False):
        """
        initialize the octopus cell, using the parameters Spencer et al. 2012
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
            nach selects the type of sodium channel that will be used in the model. A channel mechanism
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
        self.status = {self.somaname: True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'ttx': ttx, 'name': 'Octopus',
                        'morphology': morphology, 'decorator': decorator, 'temperature': None}
        self.i_test_range=(-4.0, 6.0, 0.25)
        self.spike_threshold = -50
        self.vrange = [-75., -63.]  # set a default vrange for searching for rmp
        
        if morphology is None:
            """
            instantiate a basic soma-only ("point") model
            """
            soma = h.Section(name="Octopus_Soma_%x" % id(self))  # one compartment of about 29000 um2
            soma.nseg = 1
            self.add_section(soma, self.somaname)
            self.set_soma_size_from_Section(self.soma)
            
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
            self.species_scaling(silent=True, species=species, modelType=modelType)  # set the default type II cell parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorate()
            self.decorated.channelValidate(self, verify=True)
#        print 'Mechanisms inserted: ', self.mechanisms
        self.get_mechs(self.soma)
#        self.cell_initialize(vrange=self.vrange)
        
        if debug:
            print("<< octopus: octopus cell model created >>")
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
        # Channel decoration and stick model from Figure 2
        # densities from Tables 2 and 3
        if modelType == 'Spencer':
#            print self.c_m
            self.c_m = 0.9
#            self.set_soma_size_from_Section(self.soma)
            totcap = self.totcap
            refarea = self.somaarea # totcap / self.c_m  # see above for units
#            self.print_soma_info()
            self._valid_temperatures = (34., )  # 34 for consistency with other mouse models, but
                                                # Spencer data used "33". This affects very slightly
                                                # the HCN channel conductance.
            if self.status['temperature'] is None:
                self.set_temperature(34.)            
            self.gBar = Params(nabar=0., #0.0407,  # S/cm2
                               nabar_ais=0.42441,
                               kltbar_ais=0.,
                               khtbar_ais=0.,
                               ihbar_ais=0.,
                               kltbar_soma=0.0407,
                               khtbar_soma=0.0061,
                               ihbar_soma=0.0076,
                               kltbar_dend=0.0027,
                               khtbar_dend=0.0,
                               ihbar_dend=0.0006,
                               khtbar_hillock=0.0,
                               kltbar_hillock=0.0,
                               ihbar_hillock=0.,
                               leakbar=0.0020,
            )
            
            self.channelMap = {
                self.somaname: {'jsrna': self.gBar.nabar, 'klt': self.gBar.kltbar_soma, 'kht': self.gBar.khtbar_soma,
                         'hcnobo': self.gBar.ihbar_soma, 'leak': self.gBar.leakbar, },
                'hillock': {'jsrna': 0., 'klt': self.gBar.kltbar_hillock, 'kht': self.gBar.khtbar_hillock,
                            'hcnobo': self.gBar.ihbar_hillock,
                            'leak': self.gBar.leakbar, },
                # axon initial segment:
                'unmyelinatedaxon': {'jsrna': self.gBar.nabar_ais, 'klt': self.gBar.kltbar_ais, 'kht': self.gBar.khtbar_ais,
                            'hcnobo': self.gBar.ihbar_ais, 'leak': self.gBar.leakbar, },
                'primarydendrite': {'jsrna': 0., 'klt': self.gBar.kltbar_dend, 'kht': self.gBar.khtbar_dend,
                         'hcnobo': self.gBar.ihbar_dend, 'leak': self.gBar.leakbar, },
            }
            
            self.distMap = {'primarydendrite': {'klt': {'gradient': 'flat', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'flat', 'gminf': 0., 'lambda': 100.},
                                     'hcnobo': {'gradient': 'flat', 'gminf': 0., 'lambda': 100.}}, # all flat with distance
                            }
            # reversal potential map
            self.channelErevMap = {
                self.somaname: {'jsrna': 55., 'klt': -70, 'kht': -70,
                         'hcnobo': -38, 'leak': -62., },
                'hillock': {'jsrna': 55., 'klt': -70, 'kht': -70,
                         'hcnobo': -38, 'leak': -62., },
                'unmyelinatedaxon': {'jsrna': 55., 'klt': -70, 'kht': -70,
                         'hcnobo': -38, 'leak': -62., },
                'primarydendrite': {'jsrna': 55., 'klt': -70, 'kht': -70,
                         'hcnobo': -38, 'leak': -62., },
            }

        else:
            raise ValueError('model type %s is not implemented' % modelType)
        self.check_temperature()

    def species_scaling(self, species='mouse', modelType='Spencer', silent=True):
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
        soma = self.soma

        if species == 'mouse' and modelType =='Spencer':
            print('Octopus: Mouse, Spencer point model - not a valid model')
            self.set_soma_size_from_Cm(25.0)
            self._valid_temperatures = (34., )
            if self.status['temperature'] is None:
                self.set_temperature(34.)
            self.print_soma_info()
#            self.adjust_na_chans(soma)
            # soma().kht.gbar = 0.0061  # nstomho(150.0, self.somaarea)  # 6.1 mmho/cm2
            # soma().klt.gbar = 0.0407  # nstomho(3196.0, self.somaarea)  #  40.7 mmho/cm2
            # soma().hcnobo.gbar = 0.0076  #nstomho(40.0, self.somaarea)  # 7.6 mmho/cm2, cf. Bal and Oertel, Spencer et al. 25 u dia cell
            # soma().leak.gbar = 0.0005  # nstomho(2.0, self.somaarea)
            self.axonsf = 1.0
        else:
            raise ValueError('Species "%s" or species-type "%s" is not recognized for octopus cells' %  (species, type))
        self.status['species'] = species
        self.status['modelType'] = modelType
        self.cell_initialize(showinfo=True)
        if not silent:
            print('set cell as: ', species)
            print(' with Vm rest = %6.3f' % self.vm0)

