from __future__ import print_function
from neuron import h
from ..util import nstomho
import numpy as np
from .cell import Cell
from ..util import Params
from .. import data

__all__ = ['Pyramidal', 'PyramidalKanold', 'PyramidalCeballos']

class Pyramidal(Cell):

    celltype = 'pyramidal'
    scaled = False
    
    @classmethod
    def create(cls, model='POK', **kwds):
        if model == 'POK':
            return PyramidalKanold(**kwds)
        if model == 'Ceballos':
            return PyramidalCeballos(**kwds)
        else:
            raise ValueError ('Pyramidal model %s is unknown', model)

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
            if terminal.cell.celltype in ['sgc', 'dstellate', 'tuberculoventral', 'cartwheel']:
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
                if 'AMPAScale' in kwds:
                    self.AMPA_gmax = self.AMPA_gmax * kwds['AMPAScale']  # allow scaling of AMPA conductances
                if 'NMDAScale' in kwds:
                    self.NMDA_gmax = self.NMDA_gmax*kwds['NMDAScale']
                return self.make_glu_psd(post_sec, terminal, self.AMPAR_gmax, self.NMDAR_gmax, loc=loc)
            elif terminal.cell.celltype == 'dstellate':  # WBI input -Voigt, Nelken, Young
                return self.make_gly_psd(post_sec, terminal, psdtype='glyfast', loc=loc)
            elif terminal.cell.celltype == 'tuberculoventral':  # TV cells talk to each other-Kuo et al.
                return self.make_gly_psd(post_sec, terminal, psdtype='glyfast', loc=loc)
            else:
                raise TypeError("Cannot make PSD for %s => %s" % 
                            (terminal.cell.celltype, self.celltype))
        else:
            raise ValueError("Unsupported psd type %s" % psd_type)

class PyramidalKanold(Pyramidal, Cell):
    """
    DCN pyramidal cell
    Kanold and Manis, 1999, 2001, 2005
    """
    def __init__(self,  morphology=None, decorator=None, nach=None,
                 ttx=False, species='rat', modelType=None, modelName=None,
                 debug=False, temperature=None):
        """
        initialize a pyramidal cell, based on the Kanold-Manis (2001) pyramidal cell model.
        Modifications to the cell can be made by calling methods below. These include
        converting to a model with modified size and conductances (experimental).
        
        Parameters
        ----------
        morphology : string (default: None)
            a file name to read the cell morphology from. If a valid file is found, a cell is constructed
            as a cable model from the hoc file.
            If None (default), the only a point model is made, exactly according to RM03.
            
        decorator : Python function (default: None)
            decorator is a function that "decorates" the morphology with ion channels according
            to a set of rules.
            If None, a default set of channels is inserted into the first soma section, and the
            rest of the structure is "bare".
        
        nach : string (default: None)
            nach selects the type of sodium channel that will be used in the model. A channel mechanim
            by that name must exist. None implies the default channel, 'napyr'.
        
        ttx : Boolean (default: False)
            If ttx is True, then the sodium channel conductance is set to 0 everywhere in the cell.
            Currently, this is not implemented.
        
        species: string (default 'rat')
            species defines the channel density that will be inserted for different models. Note that
            if a decorator function is specified, this argument is ignored (overridden by decorator).

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
        super(PyramidalKanold, self).__init__()
        if modelType == None or modelType == 'I':
            modelName = 'POK'
            modelType = 'pyramidal'
            dataset = 'POK_channels'
            temp = 32.  # as defined in the ms. and in the hoc files on modelDB

        else:
            raise ValueError(f"Species {species:s} and modeltype {modelType:s} not recognized for {self.celltype:s} cells")

        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'modelName': modelName, 'ttx': ttx, 'name': 'Pyramidal',
                       'morphology': morphology, 'decorator': decorator, 'temperature': None,
                   }
        self.debug=debug
        self._valid_temperatures = (temp, )
        if self.status['temperature'] == None:
            self.status['temperature'] = temp

        soma = self.do_morphology(morphology)

        self.pars = self.get_cellpars(dataset, species=species, modelType=modelType)
        self.status['na'] = self.pars.natype

        # decorate the morphology with ion channels
        if decorator is None:   # basic model, only on the soma
            self.mechanisms = ['napyr', 'kdpyr', 'kif', 'kis', 'ihpyr', 'leak']
            for mech in self.mechanisms:
                try:
                    self.soma.insert(mech)
                except ValueError:
                    print('WARNING: Mechanism %s not found' % mech)
            self.soma().kif.kif_ivh = -89.6
            self.species_scaling(silent=True)  # set the default type I-c  cell parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorate()
        self.save_all_mechs()  # save all mechanisms inserted, location and gbar values...
        self.get_mechs(self.soma)
        if debug:
            print("<< PYR: POK Pyramidal Cell created >>")


    def get_cellpars(self, dataset, species='guineapig', modelType='II'):
        cellcap = data.get(dataset, species=species, model_type=modelType,
            field='soma_Cap')
        chtype = data.get(dataset, species=species, model_type=modelType,
            field='soma_natype')
        pars = Params(cap=cellcap, natype=chtype)
        for g in ['soma_napyr_gbar', 'soma_kdpyr_gbar', 'soma_kif_gbar', 'soma_kis_gbar',
                  'soma_ihpyr_gbar', 'soma_leak_gbar',
                  'soma_e_h','soma_leak_erev', 'soma_e_k', 'soma_e_na']:
            pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
            field=g))
        if self.debug:
            pars.show()
        return pars

    def species_scaling(self, silent=True):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        Used ONLY for point models.
        
        Parameters
        ----------
        species : string (default: 'rat')
            name of the species to use for scaling the conductances in the base point model
            Must be 'rat'
        
        modelType: string (default: 'I')
            definition of model type from Kanold and Manis, 2001
            choices are 'I' or 'POK' (canonical model) or
            'II', a modified model with more physiological surface area and KCNQ channels
        
        silent : boolean (default: True)
            run silently (True) or verbosely (False)
        """
        assert self.scaled is False  # block double scaling!
        self.scaled = True

        soma = self.soma
        if self.status['species'] in ['rat', 'mouse']:
            if self.status['modelType'] not in ['pyramidal']:  # canonical K&M2001 model cell
                raise ValueError(f"\nModel type {self.status['modelType']:s} is not implemented for mouse {self.celltype.title():s} cells")
            
            self._valid_temperatures = (32.,)
            if self.status['temperature'] is None:
              self.set_temperature(32.)
            self.i_test_range = {'pulse': (-0.3, 0.401, 0.02)}
            self.vrange = [-75., -60.]
            self.set_soma_size_from_Cm(self.pars.cap)
            soma().napyr.gbar = nstomho(self.pars.soma_napyr_gbar, self.somaarea)
            # soma().nap.gbar = nstomho(self.pars.soma_nap_gbar, self.somaarea) # does not exist in canonical model
            soma().kdpyr.gbar = nstomho(self.pars.soma_kdpyr_gbar, self.somaarea)
            # soma().kcnq.gbar = nstomho(self.pars.soma_kcnq_gbar, self.somaarea) # does not exist in canonical model.
            # soma().kpksk.gbar = nstomho(self.pars.soma_kpksk_gbar, self.somaarea) # does not exist in canonical model.
            # soma().kir.gbar = nstomho(self.pars.soma_kir_gbar, self.somaarea)
            soma().kif.gbar = nstomho(self.pars.soma_kif_gbar, self.somaarea)
            soma().kis.gbar = nstomho(self.pars.soma_kis_gbar, self.somaarea)
            soma().ihpyr.gbar = nstomho(self.pars.soma_ihpyr_gbar, self.somaarea)
#            soma().ihpyr_adj.q10 = 3.0  # no temp scaling to sta
            soma().leak.gbar = nstomho(self.pars.soma_leak_gbar, self.somaarea)
            soma().leak.erev = self.pars.soma_leak_erev
            soma().ena = self.pars.soma_e_na
            soma().ek = self.pars.soma_e_k
            soma().ihpyr.eh = self.pars.soma_e_h

        # elif species in 'rat' and modelType == 'II':
        #     """
        #     Modified canonical K&M2001 model cell
        #     In this model version, the specific membrane capacitance is modified
        #     so that the overall membrane time constant is consistent with experimental
        #     measures in slices. However, this is not a physiological value. Attempts
        #     to use the normal 1 uF/cm2 value were unsuccessful in establishing the expected
        #     ~12 msec time constant.
        #     This model also adds a KCNQ channel, as described by Li et al., 2012.
        #     """
        #     self.c_m = 6.0
        #     self.set_soma_size_from_Diam(30.0)
        #     # self.set_soma_size_from_Cm(80.0)
        #     # print 'diameter: %7.1f' % self.soma.diam
        #     self._valid_temperatures = (34.,)
        #     if self.status['temperature'] is None:
        #         self.set_temperature(34.)
        #     self.refarea = self.somaarea
        #     soma().napyr.gbar = nstomho(550, self.refarea)
        #     soma().nap.gbar = nstomho(60.0, self.refarea)
        #     soma().kcnq.gbar = nstomho(2, self.refarea)  # pyramidal cells have kcnq: Li et al, 2011 (Thanos)
        #     soma().kdpyr.gbar = nstomho(180, self.refarea) # Normally 80.
        #     soma().kif.gbar = nstomho(150, self.refarea) # normally 150
        #     soma().kis.gbar = nstomho(40, self.refarea) # 40
        #     soma().ihpyr.gbar = nstomho(2.8, self.refarea)
        #     soma().leak.gbar = nstomho(0.5, self.refarea)
        #     soma().leak.erev = -62.  # override default values in cell.py
        #     soma().ena = 50.0
        #     soma().ek = -81.5
        #     soma().ihpyr.eh = -43
        #     if not self.status['dendrites']:
        #         self.add_dendrites()

        else:
            raise ValueError(f"Species {self.status['species']:s} or species-type {self.status['modelType']:s} is not recognized for T-stellate cells")

#        self.cell_initialize(showinfo=True)
        self.check_temperature()
        if not silent:
            print(f"Set cell as: {self.status['species']:s}, {self.status['modelType']:s}")
            print(self.status)
            for m in self.mechanisms:
                print('%s.gbar = %f' % (m, eval('soma().%s.gbar' % m)))

    def _ghk(self, V, ci, co, Z, mc):
        """
        GHK equation - duplicate what is in .mod file
        """
        F = 9.6485e4 #  (coul)
        R = 8.3145 # (joule/degC)
        T = h.celsius + 273.19  # Kelvin
        E = (1e-3) * V
        Ci = ci + (mc.monovalPerm) * (mc.monovalConc)        ##Monovalent permeability
        if (np.fabs(1.0-np.exp(-Z*(F*E)/(R*T))) < 1e-6): # denominator is small -> Taylor series
            ghk = (1e-6) * Z * F * (Ci-co*np.exp(-z*(F*E)/(R*T)))*(1-(z*(F*E)/(R*T)))
        else:
            ghk = (1e-6) * Z**2*(E*F**2)/(R*T)*(Ci-co*np.exp(-Z*(F*E)/(R*T)))/(1-np.exp(-Z*(F*E)/(R*T)))
        return ghk
    
    def i_currents(self, V):
        """
        For the steady-state case, return the total current at voltage V
        Used to find the zero current point
        vrange brackets the interval
        Overrides i_currents in cells.py because we have a different set of currents
        to compute.
        """
        for part in self.all_sections.keys():
            for sec in self.all_sections[part]:
                sec.v = V
        h.celsius = self.status['temperature']
        h.finitialize()
        self.ix = {}

        if 'napyr' in self.mechanisms:
             self.ix['napyr'] = self.soma().napyr.gna*(V - self.soma().ena)
        # if 'nap' in self.mechanisms:
        #      self.ix['nap'] = self.soma().nap.gnap*(V - self.soma().ena)
        # if 'cap' in self.mechanisms:
        #      mc = self.soma().cap
        #      self.ix['cap'] =mc.pcabar*mc.m*self._ghk(V, self.soma().cai, self.soma().cao, 2, mc) # (V - self.soma().ena)
        if 'kdpyr' in self.mechanisms:
             self.ix['kdpyr'] = self.soma().kdpyr.gk*(V - self.soma().ek)
        if 'kif' in self.mechanisms:
             self.ix['kif'] = self.soma().kif.gkif*(V - self.soma().ek)
        if 'kis' in self.mechanisms:
             self.ix['kis'] = self.soma().kis.gkis*(V - self.soma().ek)
        # if 'kcnq' in self.mechanisms:
        #      self.ix['kcnq'] = self.soma().kcnq.gk*(V - self.soma().ek)
        # if 'kpksk' in self.mechanisms:
        #      self.ix['kpksk'] = self.soma().kpksk.gk*(V - self.soma().ek)
        # if 'kir' in self.mechanisms:
        #      self.ix['kir'] = self.soma().kir.gk*(V - self.soma().ek)
        if 'ihpyr' in self.mechanisms:
             self.ix['ihpyr'] = self.soma().ihpyr.gh*(V - self.soma().ihpyr.eh)
        # if 'ihpyr_adj' in self.mechanisms:
        #      self.ix['ihpyr_adj'] = self.soma().ihpyr_adj.gh*(V - self.soma().ihpyr_adj.eh)
        # leak
        if 'leak' in self.mechanisms:
            self.ix['leak'] = self.soma().leak.gbar*(V - self.soma().leak.erev)
        return np.sum([self.ix[i] for i in self.ix])

    def add_dendrites(self):
        """
        Add simple unbranched dendrite.
        The dendrites have some kd, kif and ih current
        """
        nDend = range(2) # these will be simple, unbranced, N=4 dendrites
        dendrites=[]
        for i in nDend:
            dendrites.append(h.Section(cell=self.soma))
        for i in nDend:
            dendrites[i].connect(self.soma)
            dendrites[i].L = 250 # length of the dendrite (not tapered)
            dendrites[i].diam = 1
            dendrites[i].cm = self.c_m
            #h('dendrites[i].diam(0:1) = 2:1') # dendrite diameter, with tapering
            dendrites[i].nseg = 21 # # segments in dendrites
            dendrites[i].Ra = 150 # ohm.cm
            dendrites[i].insert('napyr')
            dendrites[i]().napyr.gbar = 0.00
            dendrites[i].insert('kdpyr')
            dendrites[i]().kdpyr.gbar = 0.002 # a little Ht
            dendrites[i].insert('kif')
            dendrites[i]().kif.gbar = 0.0001 # a little Ht
            dendrites[i].insert('leak') # leak
            dendrites[i]().leak.gbar = 0.00001
            dendrites[i].insert('ihpyr') # some H current
            dendrites[i]().ihvcn.gbar =  0. # 0.00002
            dendrites[i]().ihvcn.eh = -43.0
        self.maindend = dendrites
        self.status['dendrites'] = True
        self.add_section(self.maindend, 'maindend')


class PyramidalCeballos(Pyramidal, Cell):
    """
    DCN pyramidal cell
    Ceballos et al., Front Cell Neurosci. 2016
    (based on POK cell, above)
    """
    def __init__(self,  morphology=None, decorator=None, nach=None,
                 ttx=False, species='mouse', modelType=None, modelName=None,
                 debug=False, temperature=None):
        """
        initialize a pyramidal cell, based on the Kanold-Manis (2001) pyramidal cell model.
        Modifications to the cell can be made by calling methods below. These include
        converting to a model with modified size and conductances (experimental).
        
        Parameters
        ----------
        morphology : string (default: None)
            a file name to read the cell morphology from. If a valid file is found, a cell is constructed
            as a cable model from the hoc file.
            If None (default), the only a point model is made, exactly according to RM03.
            
        decorator : Python function (default: None)
            decorator is a function that "decorates" the morphology with ion channels according
            to a set of rules.
            If None, a default set of channels is inserted into the first soma section, and the
            rest of the structure is "bare".
        
        nach : string (default: None)
            nach selects the type of sodium channel that will be used in the model. A channel mechanim
            by that name must exist. None implies the default channel, 'napyr'.
        
        ttx : Boolean (default: False)
            If ttx is True, then the sodium channel conductance is set to 0 everywhere in the cell.
            Currently, this is not implemented.
        
        species: string (default 'rat')
            species defines the channel density that will be inserted for different models. Note that
            if a decorator function is specified, this argument is ignored (overridden by decorator).

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
        super(PyramidalCeballos, self).__init__()
        if modelType == None or modelType == 'I':
            modelName = 'Ceballos'
            modelType = 'pyramidal'
            dataset = 'Ceballos_channels'
            temp = 34.

        else:
            raise ValueError(f"Species {species:s} and modeltype {modelType:s} not recognized for {self.celltype:s} cells")

        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'modelName': modelName, 'ttx': ttx, 'name': 'Pyramidal',
                       'morphology': morphology, 'decorator': decorator, 'temperature': None,
                   }
        self.debug=debug
        self._valid_temperatures = (temp, )
        if self.status['temperature'] == None:
            self.status['temperature'] = temp

        soma = self.do_morphology(morphology)

        self.pars = self.get_cellpars(dataset, species=species, modelType=modelType)
        self.status['na'] = self.pars.natype

        # decorate the morphology with ion channels
        if decorator is None:   # basic model, only on the soma
            self.mechanisms = ['napyr', 'kdpyr', 'kif', 'kis', 'ihpyr', 'leak', 'kcnq', 'kir', 'nap', 'kpksk', 'cap']
            for mech in self.mechanisms:
                try:
                    self.soma.insert(mech)
                except ValueError:
                    print('WARNING: Mechanism %s not found' % mech)
            self.soma().kif.kif_ivh = -89.6
            self.species_scaling(silent=True)  # set the default type I-c  cell parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorate()
        self.save_all_mechs()  # save all mechanisms inserted, location and gbar values...
        self.get_mechs(self.soma)
        if debug:
            print("<< PYR: Ceballos Pyramidal Cell created >>")


    def get_cellpars(self, dataset, species='mouse', modelType='I'):
        cellcap = data.get(dataset, species=species, model_type=modelType,
            field='soma_Cap')
        chtype = data.get(dataset, species=species, model_type=modelType,
            field='soma_natype')
        pars = Params(cap=cellcap, natype=chtype)
        for g in ['soma_napyr_gbar', 'soma_kdpyr_gbar', 'soma_kif_gbar', 'soma_kis_gbar',
                  'soma_kcnq_gbar', 'soma_kir_gbar', 'soma_kpksk_gbar', 'soma_cap_pcabar',
                  'soma_nap_gbar', 'soma_ihpyr_gbar', 'soma_leak_gbar',
                  'soma_e_h','soma_leak_erev', 'soma_e_k', 'soma_e_na']:
            pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
            field=g))
        if self.debug:
            pars.show()
        return pars

    def species_scaling(self, silent=True):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        Used ONLY for point models.
        
        Parameters
        ----------
        species : string (default: 'rat')
            name of the species to use for scaling the conductances in the base point model
            Must be 'rat'
        
        modelType: string (default: 'I')
            definition of model type from Kanold and Manis, 2001
            choices are 'I' or 'POK' (canonical model) or
            'II', a modified model with more physiological surface area and KCNQ channels
        
        silent : boolean (default: True)
            run silently (True) or verbosely (False)
        """
        assert self.scaled is False  # block double scaling!
        self.scaled = True

        soma = self.soma
        if self.status['species'] in ['rat', 'mouse']:
            if self.status['modelType'] not in ['pyramidal']:  # canonical K&M2001 model cell
                raise ValueError(f"\nModel type {self.status['modelType']:s} is not implemented for mouse {self.celltype.title():s} cells")
            
            self._valid_temperatures = (34.,)
            if self.status['temperature'] is None:
              self.set_temperature(34.)
            self.i_test_range = {'pulse': (-0.3, 0.401, 0.02)}
            self.vrange = [-75., -60.]
            self.set_soma_size_from_Cm(self.pars.cap)
            soma().napyr.gbar = nstomho(self.pars.soma_napyr_gbar, self.somaarea)
            soma().nap.gbar = nstomho(self.pars.soma_nap_gbar, self.somaarea) # does not exist in canonical model
            soma().kdpyr.gbar = nstomho(self.pars.soma_kdpyr_gbar, self.somaarea)
            soma().kcnq.gbar = nstomho(self.pars.soma_kcnq_gbar, self.somaarea) # does not exist in canonical model.
            soma().kpksk.gbar = nstomho(self.pars.soma_kpksk_gbar, self.somaarea) # does not exist in canonical model.
            soma().kir.gbar = nstomho(self.pars.soma_kir_gbar, self.somaarea)
            soma().kif.gbar = nstomho(self.pars.soma_kif_gbar, self.somaarea)
            soma().kis.gbar = nstomho(self.pars.soma_kis_gbar, self.somaarea)
            soma().ihpyr.gbar = nstomho(self.pars.soma_ihpyr_gbar, self.somaarea)
#            soma().ihpyr_adj.q10 = 3.0  # no temp scaling to sta
            soma().leak.gbar = nstomho(self.pars.soma_leak_gbar, self.somaarea)
            soma().leak.erev = self.pars.soma_leak_erev
            soma().ena = self.pars.soma_e_na
            soma().ek = self.pars.soma_e_k
            soma().ihpyr.eh = self.pars.soma_e_h

        # elif species in 'rat' and modelType == 'II':
        #     """
        #     Modified canonical K&M2001 model cell
        #     In this model version, the specific membrane capacitance is modified
        #     so that the overall membrane time constant is consistent with experimental
        #     measures in slices. However, this is not a physiological value. Attempts
        #     to use the normal 1 uF/cm2 value were unsuccessful in establishing the expected
        #     ~12 msec time constant.
        #     This model also adds a KCNQ channel, as described by Li et al., 2012.
        #     """
        #     self.c_m = 6.0
        #     self.set_soma_size_from_Diam(30.0)
        #     # self.set_soma_size_from_Cm(80.0)
        #     # print 'diameter: %7.1f' % self.soma.diam
        #     self._valid_temperatures = (34.,)
        #     if self.status['temperature'] is None:
        #         self.set_temperature(34.)
        #     self.refarea = self.somaarea
        #     soma().napyr.gbar = nstomho(550, self.refarea)
        #     soma().nap.gbar = nstomho(60.0, self.refarea)
        #     soma().kcnq.gbar = nstomho(2, self.refarea)  # pyramidal cells have kcnq: Li et al, 2011 (Thanos)
        #     soma().kdpyr.gbar = nstomho(180, self.refarea) # Normally 80.
        #     soma().kif.gbar = nstomho(150, self.refarea) # normally 150
        #     soma().kis.gbar = nstomho(40, self.refarea) # 40
        #     soma().ihpyr.gbar = nstomho(2.8, self.refarea)
        #     soma().leak.gbar = nstomho(0.5, self.refarea)
        #     soma().leak.erev = -62.  # override default values in cell.py
        #     soma().ena = 50.0
        #     soma().ek = -81.5
        #     soma().ihpyr.eh = -43
        #     if not self.status['dendrites']:
        #         self.add_dendrites()

        else:
            raise ValueError(f"Species {self.status['species']:s} or species-type {self.status['modelType']:s} is not recognized for T-stellate cells")

#        self.cell_initialize(showinfo=True)
        self.check_temperature()
        if not silent:
            print(f"Set cell as: {self.status['species']:s}, {self.status['modelType']:s}")
            print(self.status)
            for m in self.mechanisms:
                print('%s.gbar = %f' % (m, eval('soma().%s.gbar' % m)))

    def _ghk(self, V, ci, co, Z, mc):
        """
        GHK equation - duplicate what is in .mod file
        """
        F = 9.6485e4 #  (coul)
        R = 8.3145 # (joule/degC)
        T = h.celsius + 273.19  # Kelvin
        E = (1e-3) * V
        Ci = ci + (mc.monovalPerm) * (mc.monovalConc)        ##Monovalent permeability
        if (np.fabs(1.0-np.exp(-Z*(F*E)/(R*T))) < 1e-6): # denominator is small -> Taylor series
            ghk = (1e-6) * Z * F * (Ci-co*np.exp(-z*(F*E)/(R*T)))*(1-(z*(F*E)/(R*T)))
        else:
            ghk = (1e-6) * Z**2*(E*F**2)/(R*T)*(Ci-co*np.exp(-Z*(F*E)/(R*T)))/(1-np.exp(-Z*(F*E)/(R*T)))
        return ghk
    
    def i_currents(self, V):
        """
        For the steady-state case, return the total current at voltage V
        Used to find the zero current point
        vrange brackets the interval
        Overrides i_currents in cells.py because we have a different set of currents
        to compute.
        """
        for part in self.all_sections.keys():
            for sec in self.all_sections[part]:
                sec.v = V
        h.celsius = self.status['temperature']
        h.finitialize()
        self.ix = {}

        if 'napyr' in self.mechanisms:
             self.ix['napyr'] = self.soma().napyr.gna*(V - self.soma().ena)
        if 'nap' in self.mechanisms:
             self.ix['nap'] = self.soma().nap.gnap*(V - self.soma().ena)
        if 'cap' in self.mechanisms:
             mc = self.soma().cap
             self.ix['cap'] =mc.pcabar*mc.m*self._ghk(V, self.soma().cai, self.soma().cao, 2, mc) # (V - self.soma().ena)
        if 'kdpyr' in self.mechanisms:
             self.ix['kdpyr'] = self.soma().kdpyr.gk*(V - self.soma().ek)
        if 'kif' in self.mechanisms:
             self.ix['kif'] = self.soma().kif.gkif*(V - self.soma().ek)
        if 'kis' in self.mechanisms:
             self.ix['kis'] = self.soma().kis.gkis*(V - self.soma().ek)
        if 'kcnq' in self.mechanisms:
             self.ix['kcnq'] = self.soma().kcnq.gk*(V - self.soma().ek)
        if 'kpksk' in self.mechanisms:
             self.ix['kpksk'] = self.soma().kpksk.gk*(V - self.soma().ek)
        if 'kir' in self.mechanisms:
             self.ix['kir'] = self.soma().kir.gk*(V - self.soma().ek)
        if 'ihpyr' in self.mechanisms:
             self.ix['ihpyr'] = self.soma().ihpyr.gh*(V - self.soma().ihpyr.eh)
        if 'ihpyr_adj' in self.mechanisms:
             self.ix['ihpyr_adj'] = self.soma().ihpyr_adj.gh*(V - self.soma().ihpyr_adj.eh)
        # leak
        if 'leak' in self.mechanisms:
            self.ix['leak'] = self.soma().leak.gbar*(V - self.soma().leak.erev)
        return np.sum([self.ix[i] for i in self.ix])

    def add_dendrites(self):
        """
        Add simple unbranched dendrite.
        The dendrites have some kd, kif and ih current
        """
        nDend = range(2) # these will be simple, unbranced, N=4 dendrites
        dendrites=[]
        for i in nDend:
            dendrites.append(h.Section(cell=self.soma))
        for i in nDend:
            dendrites[i].connect(self.soma)
            dendrites[i].L = 250 # length of the dendrite (not tapered)
            dendrites[i].diam = 1
            dendrites[i].cm = self.c_m
            #h('dendrites[i].diam(0:1) = 2:1') # dendrite diameter, with tapering
            dendrites[i].nseg = 21 # # segments in dendrites
            dendrites[i].Ra = 150 # ohm.cm
            dendrites[i].insert('napyr')
            dendrites[i]().napyr.gbar = 0.00
            dendrites[i].insert('kdpyr')
            dendrites[i]().kdpyr.gbar = 0.002 # a little Ht
            dendrites[i].insert('kif')
            dendrites[i]().kif.gbar = 0.0001 # a little Ht
            dendrites[i].insert('leak') # leak
            dendrites[i]().leak.gbar = 0.00001
            dendrites[i].insert('ihpyr_adj') # some H current
            dendrites[i]().ihvcn.gbar =  0. # 0.00002
            dendrites[i]().ihvcn.eh = -43.0
        self.maindend = dendrites
        self.status['dendrites'] = True
        self.add_section(self.maindend, 'maindend')
    