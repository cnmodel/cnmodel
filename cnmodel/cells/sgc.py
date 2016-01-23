from neuron import h
import neuron as nrn
from ..util import nstomho
import numpy as np
from .cell import Cell
from .. import synapses
from .. import an_model
from .. import data

__all__ = ['SGC', 'SGC_TypeI', 'DummySGC']


class SGC(Cell):
    type = 'sgc'

    @classmethod
    def create(cls, model='I', species='mouse', **kwds):
        if model == 'dummy':
            return DummySGC(**kwds)
        elif model == 'I':
            return SGC_TypeI(species=species, **kwds)
        else:
            raise ValueError ('SGC model %s is unknown', model)
        
    def __init__(self, cf=None, sr=None):
        Cell.__init__(self)
        self._cf = cf
        self._sr = sr
        self.spike_source = None  # used by DummySGC to connect VecStim to terminal

    @property
    def cf(self):
        """ Center frequency
        """
        return self._cf
    
    @property
    def sr(self):
        """ Spontaneous rate group. 1=low, 2=mid, 3=high
        """
        return self._sr
        
    def make_terminal(self, post_cell, term_type, **kwds):
        """Create a StochasticTerminal and configure it according to the 
        postsynaptic cell type.
        """
        pre_sec = self.soma
        
        # Return a simple terminal unless a stochastic terminal was requested.
        if term_type == 'simple':
            return synapses.SimpleTerminal(pre_sec, post_cell, 
                                           spike_source=self.spike_source, **kwds)
        elif term_type == 'multisite':
            n_rsites = data.get('sgc_synapse', species='mouse', post_type=post_cell.type,
                            field='n_rsites')
            
            opts = {'nzones': n_rsites, 'delay': 0, 'dep_flag' : 1}
            opts.update(kwds)
            # when created, depflag is set True (1) so that we compute the DKR D*F to get release
            # this can be modified prior to the run by setting the terminal(s) so that dep_flag is 0
            # (no DKR: constant release probability)
            term = synapses.StochasticTerminal(pre_sec, post_cell,
                    spike_source=self.spike_source, **opts)
            
            kinetics = data.get('sgc_ampa_kinetics', species='mouse', post_type=post_cell.type,
                            field=['tau_g', 'amp_g'])
            term.set_params(**kinetics)
            dynamics = data.get('sgc_release_dynamics', species='mouse', post_type=post_cell.type,
                                field=['F', 'k0', 'kmax', 'kd', 'kf', 'taud', 'tauf', 'dD', 'dF'])
            term.set_params(**dynamics)
            return term
        else:
            raise ValueError("Unsupported terminal type %s" % term_type)
    

class DummySGC(SGC):
    """ SGC class with no cell body; this cell only replays a predetermined
    spike train.
    """
    def __init__(self, cf=None, sr=None):
        SGC.__init__(self, cf, sr)
        self.vecstim = h.VecStim()
        
        # this causes the terminal to receive events from the VecStim:
        self.spike_source = self.vecstim
        
        # just an empty section for holding the terminal
        self.add_section(h.Section(), 'soma')
        
    def set_spiketrain(self, times):
        """ Set the times of spikes to be replayed by the cell.
        """
        self._spiketrain = times
        self._stvec = h.Vector(times)
        self.vecstim.play(self._stvec)

    def set_sound_stim(self, stim, seed):
        """ Set the sound stimulus used to generate this cell's spike train.
        """
        self._sound_stim = stim
        spikes = an_model.get_spiketrain(cf=self.cf, sr=self.sr, seed=seed, stim=stim)
        self.set_spiketrain(spikes * 1000)


class SGC_TypeI(SGC):
    """
    Spiral ganglion cell model
    """
    def __init__(self, morphology=None, decorator=None, nach='jsrna', ttx=False,
                 debug=False, species='guineapig', 
                 modelType='bm', cf=None, sr=None):
        super(SGC_TypeI, self).__init__(cf=cf, sr=sr)

        if modelType == None:
            modelType = 'bm'  # modelTypes are: a (apical), bm (basal middle)
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'ttx': ttx, 'name': 'SGC',
                        'morphology': morphology, 'decorator': decorator}

        self.i_test_range=[(-0.3, 0.3, 0.02), (-0.03, 0., 0.005)]

        if morphology is None:
            """
            instantiate a basic soma-only ("point") model
            """
            soma = h.Section(name="SGC_Soma_%x" % id(self)) # one compartment of about 29000 um2
            soma.nseg = 1
        else:
            """
            instantiate a structured model with the morphology as specified by 
            the morphology file
            """
            soma = self.morphology_from_hoc(morphology=morphology, somasection='sections[0]')

        # decorate the morphology with ion channels
        if decorator is None:   # basic model, only on the soma
            self.mechanisms = [nach, 'klt', 'kht', 'leak']
            if modelType == 'a':
                self.mechanisms.append('ihsgcApical')
            elif modelType == 'bm':
                self.mechanisms.append('ihsgcBasalMiddle')
            else:
                raise ValueError ('Type %s not know for SGC model' % modelType)
            for mech in self.mechanisms:
                soma.insert(mech)
            soma.ek = self.e_k
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
        self.cell_initialize()
        if debug:
            print "<< SGC: Spiral Ganglion Cell created >>"

    def species_scaling(self, silent=True, species='guineapig', modelType='a'):
        soma = self.soma
        if species == 'mouse':
            self.set_soma_size_from_Cm(12.0)
            self.adjust_na_chans(soma, gbar=350.)
            soma().kht.gbar = nstomho(58.0, self.somaarea)
            soma().klt.gbar = nstomho(80.0, self.somaarea)
                # nstomho(200.0, somaarea) * scalefactor
            if modelType == 'a':
                soma().ihsgcApical.gbar = nstomho(3.0, self.somaarea)
                soma().ihsgcApical.eh = -41
            elif modelType == 'bm':
                soma().ihsgcBasalMiddle.gbar = nstomho(3.0, self.somaarea)
                soma().ihsgcBasalMiddle.eh = -41
            else:
                raise ValueError('Ihsgc modelType %s not recognized for species %s' % (species, modelType))
            soma().leak.gbar = nstomho(2.0, self.somaarea)

        elif species == 'guineapig':
            # guinea pig data from Rothman and Manis, 2003, modelType II
            self.set_soma_size_from_Cm(12.0)
            self.adjust_na_chans(soma, gbar=1000.)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().klt.gbar = nstomho(200.0, self.somaarea)
                # nstomho(200.0, somaarea) * scalefactor
            if modelType == 'a':
                soma().ihsgcApical.gbar = nstomho(3.0, self.somaarea)
                soma().ihsgcApical.eh = -41
            elif modelType == 'bm':
                soma().ihsgcBasalMiddle.gbar = nstomho(3.0, self.somaarea)
                soma().ihsgcBasalMiddle.eh = -41
            else:
                raise ValueError('Ihsgc modelType %s not recognized for species %s' % (species, modelType))
            soma().leak.gbar = nstomho(2.0, self.somaarea)

        else:
            raise ValueError('Species %s or species-modelType %s is not recognized for SGC cells' % (species, modelType))

        self.status['species'] = species
        self.status['modelType'] = modelType
        self.cell_initialize(showinfo=False)
        if not silent:
            print 'set cell as: ', species
            print ' with Vm rest = %f' % self.vm0

    def adjust_na_chans(self, soma, gbar = 1000., debug=False):
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
        if 'ihsgcApical' in self.mechanisms:
            self.ix['ihsgcApical'] = self.soma().ihsgcApical.gh*(V - self.soma().ihsgcApical.eh)
        if 'ihsgcBasalMiddle' in self.mechanisms:
            self.ix['ihsgcBasalMiddle'] = self.soma().ihsgcBasalMiddle.gh*(V - self.soma().ihsgcBasalMiddle.eh)
        if 'leak' in self.mechanisms:
            self.ix['leak'] = self.soma().leak.gbar*(V - self.soma().leak.erev)
#        print self.status['name'], self.status['type'], V, self.ix
        return np.sum([self.ix[i] for i in self.ix])

