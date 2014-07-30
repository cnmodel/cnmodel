from neuron import h
import numpy as np
import neuron as nrn

from .cell import Cell
from .. import synapses
from ..util import nstomho

__all__ = ['TStellate', 'TStellateNav11', 'TStellateFast'] 


class TStellate(Cell):
    """
    VCN T-stellate base model.
    Rothman and Manis, 2003abc (Type I-c, Type I-t)
    """
    def __init__(self, nach='na', ttx=False, debug=False, species='guineapig', type=None):
        """
        initialize a planar stellate (T-stellate) cell, using the default parameters for guinea pig from
        R&M2003, as a type I cell.
        Modifications to the cell can be made by calling methods below. These include:
            Converting to a type IA model (add transient K current) (species: guineapig-TypeIA).
            Changing "species" to mouse or cat (scales conductances)
        """
        super(TStellate, self).__init__()
        if type == None:
            type = 'I-c'
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'type': type, 'ttx': ttx, 'name': 'TStellate'}

        self.i_test_range=(-0.15, 0.15, 0.01)

        soma = h.Section(name="TStellate_Soma_%x" % id(self)) # one compartment of about 29000 um2

        soma.nseg = 1

        if nach in ['nacn', 'na']:
            soma.insert('na')
        elif nach == 'nav11':
            soma.insert('nav11')
        elif nach == 'jsrna':
            soma.insert('jsrna')
        else:
            raise ValueError('Sodium channel %s in type 1 cell not known' % nach)
        self.mechanisms = ['kht', 'ka', 'ihvcn', 'leak', nach]
        for mech in self.mechanisms:
            soma.insert(mech)
        # soma.insert("kht")
        # soma.insert('ka')
        # soma.insert('ihvcn')
        # soma.insert('leak')
        soma.ek = self.e_k
        soma.ena = self.e_na
        soma().ihvcn.eh = self.e_h
        soma().leak.erev = self.e_leak
        self.add_section(soma, 'soma')
        self.species_scaling(silent=True, species=species, type=type)  # set the default type I-c  cell parameters
        self.get_mechs(soma)
        self.cell_initialize()
        if debug:
                print "<< T-stellate: JSR Stellate Type 1 cell model created >>"

    def species_scaling(self, species='guineapig', type='I-c', silent=True):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        """
        soma = self.soma
        if species == 'mouse' and type == 'I-c':
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
        elif species == 'guineapig' and type == 'I-c':  # values from R&M 2003, Type I
            self.set_soma_size_from_Cm(12.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().ka.gbar = nstomho(0.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(0.5, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            soma().leak.erev = -65.0
            self.axonsf = 0.5
        elif species == 'guineapig' and type =='I-t':
            # guinea pig data from Rothman and Manis, 2003, type It
            self.set_soma_size_from_Cm(12.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(80.0, self.somaarea)
            soma().ka.gbar = nstomho(65.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(0.5, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            soma().leak.erev = -65.0
            self.axonsf = 0.5
        elif species == 'cat' and type == 'I-c':  # a cat is a big guinea pig Type I
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
        self.status['type'] = type
        self.cell_initialize(showinfo=True)
        if not silent:
            print 'set cell as: ', species
            print ' with Vm rest = %f' % self.vm0

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
            raise ValueError("Dstellate setting Na channels: channel %s not known" % nach)


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

    def make_psd(self, pre_sec, post_sec, terminal, **kwds):
        from .. import cells
        pre_cell = cells.cell_from_section(pre_sec)
        if isinstance(pre_cell, cells.SGC):
            return synapses.GluPSD(pre_sec, post_sec, terminal, 
                                   ampa_gmax=4600.,
                                   nmda_ampa_ratio = 1.28,
                                   )
        else:
            raise TypeError("Cannot make PSD for %s => %s" % 
                            (pre_cell.__class__.__name__, 
                             self.__class__.__name__))



class TStellateNav11(Cell):
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
        print ("T-STELLATE ROTHMAN",
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

class TStellateFast(Cell):
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

