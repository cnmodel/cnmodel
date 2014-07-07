from neuron import h
import neuron as nrn
from ..pynrnutilities import nstomho

from .cell import Cell

__all__ = ['TStellate', 'TStellateNav11', 'TStellateFast'] 

class TStellate(Cell):
    """
    VCN T-stellate model.
    
    """
    def __init__(self, debug=False, ttx=False, message=None,
        species='guinea pig', nav11=False):
        soma = h.Section() # one compartment of about 29000 um2
        v_potassium = -80       # potassium reversal potential
        v_sodium = 50           # sodium reversal potential

        cm = 1.0
        scalefactor = 1.0 # This determines the relative size of the cell
        rinsf = 1.0           # input resistance adjustment (also current...)
        if species == 'guinea pig':
            totcap = scalefactor * 12.0 # cap in pF for cell
        elif species == 'mouse':
            totcap = scalefactor * 25.0 # adjusted...
        else:
            raise ValueError('Unrecognized species in call to tstellate_rothman')

        effcap = totcap # sometimes we change capacitance - that's effcap
        somaarea = totcap * 1E-6 / cm # pf -> uF, cm = 1uf/cm^2 nominal
        lstd = 1E4 * ((somaarea / 3.14159) ** 0.5) # convert from cm to um

        soma.nseg = 1
        soma.diam = lstd
        soma.L = lstd

        seg = soma
        seg.insert("kht")
        if nav11 is False:
            seg.insert('na')
            seg.ena = 50
        else:
            seg.insert('nav11')
            seg.ena = 50
        seg.insert('ka')
        seg.insert('ihvcn')
        seg.insert('leak')
        seg.ek = -84
        s = soma()
        if species == 'guinea pig':
            gnamax = 1000.0
            s.kht.gbar = nstomho(150.0, somaarea) * scalefactor
            s.ka.gbar = nstomho(0.0, somaarea) * scalefactor
            s.ihvcn.gbar = nstomho(0.5, somaarea) * scalefactor
            s.ihvcn.eh = -43 # Rodrigues and Oertel, 2006
            #print 'ih vcn vh: %f ' % (s.ihvcn.vh)
            s.leak.gbar = nstomho(2.0, somaarea) * scalefactor
            vm0 = -63.9
        elif species == 'mouse':
            gnamax = 800.0
            s.kht.gbar = nstomho(250.0, somaarea) * scalefactor
            s.ka.gbar = nstomho(0.0, somaarea) * scalefactor
            s.ihvcn.gbar = nstomho(18.0, somaarea) * scalefactor
            s.ihvcn.eh = -43 # Rodrigues and Oertel, 2006
            s.leak.gbar = nstomho(8.0, somaarea) * scalefactor
            # yields input resistance of 74.2 Mohm measured from -60 to -70 mV
            vm0 = -60.0
        else:
            raise ValueError('species not recognized in tstellate_rothman')

        if not ttx:
            if nav11 is False:
                s.na.gna = nstomho(gnamax, somaarea) * scalefactor
            else:
                s.nav11.gbar = nstomho(1800.0, somaarea) * scalefactor
                s.nav11.vsna = 4.3 # was 8
        else:
            if nav11 is False:
                s.na.gna = 0.0
            else:
                s.nav11.gbar = 0.0

        if debug:
            if message is None:
                print "<< T-stellate: JSR Stellate Type 1 cell model created >>"
            else:
                print message
        
        self.soma = soma

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
        self.soma = soma
        self.dendrites = dendrites

class TStellateFast(Cell):
    """ 
    VCN t-stellate model based on Rothman and Manis 2003, but with fast sodium
    channel 
    """
    def __init__(self, debug=False, ttx=False, message=None, dend=False):
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
        
        self.soma = soma
