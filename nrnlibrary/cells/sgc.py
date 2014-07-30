from neuron import h
import neuron as nrn
from ..util import nstomho
import numpy as np
from .cell import Cell
from .. import synapses

__all__ = ['SGC']


class SGC(Cell):
    """
    Spiral ganglion cell model
    """
    def __init__(self, nach='jsrna', ttx=False, debug=False, species='guineapig', type='bm'):
        super(SGC, self).__init__()

        if type == None:
            type = 'bm'  # types are: a (apical), bm (basal middle)
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'type': type, 'ttx': ttx, 'name': 'SGC'}

        self.i_test_range=(-0.3, 0.3, 0.02)

        soma = h.Section(name="SGC_Soma_%x" % id(self)) # one compartment of about 29000 um2

        soma.nseg = 1

        self.mechanisms = [nach, 'klt', 'kht', 'leak']
        if type == 'a':
            self.mechanisms.append('ihsgcApical')
        elif type == 'bm':
            self.mechanisms.append('ihsgcBasalMiddle')
        else:
            raise ValueError ('Type %s not know for SGC model' % type)
        for mech in self.mechanisms:
            soma.insert(mech)
        soma.ek = self.e_k
        soma().leak.erev = self.e_leak

        self.add_section(soma, 'soma')
        self.species_scaling(silent=True, species=species, type=type)  # set the default type I-c  cell parameters
        self.get_mechs(soma)
#        self.cell_initialize()
        if debug:
            print "<< SGC: Spiral Ganglion Cell created >>"

    def species_scaling(self, silent=True, species='guineapig', type='a'):
        soma = self.soma
        if species == 'mouse':
            self.set_soma_size_from_Cm(12.0)
            self.adjust_na_chans(soma, gbar=350.)
            soma().kht.gbar = nstomho(58.0, self.somaarea)
            soma().klt.gbar = nstomho(80.0, self.somaarea)
                # nstomho(200.0, somaarea) * scalefactor
            if type == 'a':
                soma().ihsgcApical.gbar = nstomho(3.0, self.somaarea)
                soma().ihsgcApical.eh = -41
            elif type == 'bm':
                soma().ihsgcBasalMiddle.gbar = nstomho(3.0, self.somaarea)
                soma().ihsgcBasalMiddle.eh = -41
            else:
                raise ValueError('Ihsgc type %s not recognized for species %s' % (species, type))
            soma().leak.gbar = nstomho(2.0, self.somaarea)

        elif species == 'guineapig':
            # guinea pig data from Rothman and Manis, 2003, type II
            self.set_soma_size_from_Cm(12.0)
            self.adjust_na_chans(soma, gbar=1000.)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().klt.gbar = nstomho(200.0, self.somaarea)
                # nstomho(200.0, somaarea) * scalefactor
            if type == 'a':
                soma().ihsgcApical.gbar = nstomho(3.0, self.somaarea)
                soma().ihsgcApical.eh = -41
            elif type == 'bm':
                soma().ihsgcBasalMiddle.gbar = nstomho(3.0, self.somaarea)
                soma().ihsgcBasalMiddle.eh = -41
            else:
                raise ValueError('Ihsgc type %s not recognized for species %s' % (species, type))
            soma().leak.gbar = nstomho(2.0, self.somaarea)

        else:
            raise ValueError('Species %s or species-type %s is not recognized for SGC cells' % (species, type))

        self.status['species'] = species
        self.status['type'] = type
        self.cell_initialize(showinfo=True)
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

    def make_terminal(self, pre_sec, post_sec, **kwds):
        from .. import cells
        post_cell = cells.cell_from_section(post_sec)
        #
        # set parameters according to the target cell type
        #
        
        if isinstance(post_cell, cells.Bushy):
            nzones, delay = 100, 0
        elif isinstance(post_cell, cells.TStellate):
            nzones, delay = 1, 0
        elif isinstance(post_cell, cells.DStellate):
            nzones, delay = 1, 0
        else:
            raise NotImplementedError("Cannot connect SGC to cell type %s" % 
                                      type(post_cell))
        
        return synapses.StochasticTerminal(pre_sec, post_sec, nzones=nzones, delay=delay)
        
