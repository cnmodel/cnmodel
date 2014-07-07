from neuron import h
import neuron as nrn
from ..pynrnutilities import nstomho

from .cell import Cell

__all__ = ['Bushy', 'BushyWithAxon'] 

class Bushy(Cell):
    def __init__(self, debug=False, ttx=False, message=None, nach='jsrna',
                 species='mouse', axon=False, dendrite=False,
                 newModFiles=False, pump=False):
        super(Bushy, self).__init__()
        
        v_potassium = -70       # potassium reversal potential
    #    v_sodium = 55           # sodium reversal potential
    #    v_chloride = -20
        c_m = 1.0
        R_a = 150
        scalefactor = 1.0 # This determines the relative size of the cell
    #    rinsf = 1.0           # input resistance adjustment (also current...)
        if species == 'guineapig-bushy-II' or species == 'guineapig-bushy-II-I':
            totcap = scalefactor * 12.0 # cap in pF for cell
    #        axonsf = 0.57
    #        isdist = 0
        if species == 'mouse':
            totcap = scalefactor * 26.0 # cap in pF for cell from Cao et al., 2007
    #        axonsf = 0.57
    #        isdist = 0
        if species == 'cat':
            totcap = scalefactor * 35.0
    #        axonsf = 1.0
    #        isdist = 0

    #    effcap = totcap # sometimes we change capacitance - that's effcap
        somaarea = totcap * 1E-6 / c_m # pf -> uF, cm = 1uf/cm^2 nominal
        lstd = 1E4 * ((somaarea / 3.14159) ** 0.5) # convert from cm to um

        soma = h.Section() # one compartment of about 29000 um2
        soma.nseg = 1
        soma.diam = lstd
        soma.L = lstd

        if newModFiles:
            soma.insert('klt2')
            soma.insert('kht2')

        if pump:
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
        else:
            soma.insert('klt')
            soma.insert('kht')

        if nach == 'jsrna':
            soma.insert('jsrna')
        elif nach == 'nav11':
            soma.insert('nav11')
        else:
            soma.insert('na')
        #soma.insert('ka')
        soma.insert('ihvcn')
        #soma.insert('hcno')
        soma.insert('leak')
        soma.ek = v_potassium
        gnabar = nstomho(1000.0, somaarea) * scalefactor
        if ttx:
            gnabar = 0.0
        if nach == 'jsrna':
            soma().jsrna.gna = gnabar
            soma.ena = 50
        elif nach == 'nav11':
            soma().nav11.gnatbar = gnabar * 0.5
            soma.ena = 50
            soma().nav11.vsna = 4.3
            print "bushy using inva11"
        else:
            soma().na.gnabar = gnabar
            soma.ena = 50

        vm0 = self.bushy_species_scaling(soma, species, newModFiles, somaarea, scalefactor, debug)

        maindend = None
        secdend = None
        if dendrite:
            print 'Adding dendrite to Bushy model'
            section = h.Section
            maindend = section(cell=soma)
            maindend.connect(soma)
            maindend.nseg = 10
            maindend.L = 100.0
            maindend.diam = 2.5
            maindend.insert('klt')
            maindend.insert('ihvcn')
            if newModFiles:
                maindend().klt2.gbar = soma().klt2.gbar / 2.0
            else:
                maindend().klt.gbar = soma().klt.gbar / 2.0

            maindend().ihvcn.gbar = soma().ihvcn.gbar / 2.0

            maindend.cm = c_m
            maindend.Ra = R_a
            nsecd = range(0, 5)
            secdend = []
            for ibd in nsecd:
                secdend.append(section(cell=soma))
            for ibd in nsecd:
                secdend[ibd].connect(maindend)
                secdend[ibd].diam = 1.0
                secdend[ibd].L = 15.0
                secdend[ibd].cm = c_m
                secdend[ibd].Ra = R_a
            #h.topology()

        if debug:
            if message is None:
                print "<< bushy: JSR bushy cell model created >>"
            else:
                print message
                
        self.soma = soma
        self.maindend = maindend
        self.secdend = secdend

    @staticmethod
    def bushy_species_scaling(soma, species, newModFiles, somaarea, scalefactor, debug):
        if species == 'mouse':
            # use conductance levels from Cao et al.,  J. Neurophys., 2007.
            #if debug:
            print 'Mouse bushy cell'
            #s.nav11.gbar = gnabar*0.3
            if newModFiles:
                soma().kht2.gbar = nstomho(58.0, somaarea) * scalefactor
                soma().klt2.gbar = nstomho(80.0, somaarea) * scalefactor
            else:
                soma().kht.gbar = nstomho(58.0, somaarea) * scalefactor
                soma().klt.gbar = nstomho(80.0, somaarea) * scalefactor
            #s.ka.gkabar = nstomho(0.0, somaarea) * scalefactor
            soma().ihvcn.gbar = nstomho(30.0, somaarea) * scalefactor
            #s.hcno.gbar = 0.0
            soma().leak.gbar = nstomho(2.0, somaarea) * scalefactor
            vm0 = -63.6
        if species == 'guineapig-bushy-II':
            # guinea pig data from Rothman and Manis, 2003, type II
            if debug:
                print 'Guinea pig bushy cell type II'
            if newModFiles:
                soma().kht2.gbar = nstomho(150.0, somaarea) * scalefactor
                soma().klt2.gbar = nstomho(200.0, somaarea) * scalefactor
            else:
                soma().kht.gbar = nstomho(150.0, somaarea) * scalefactor
                soma().klt.gbar = nstomho(200.0, somaarea) * scalefactor
            #s.ka.gkabar = nstomho(0.0, somaarea) * scalefactor
            soma().ihvcn.gbar = nstomho(20.0, somaarea) * scalefactor
            #s.hcno.gbar = 0.0
            soma().leak.gbar = nstomho(2.0, somaarea) * scalefactor
            vm0 = -63.6
        if species == 'guineapig-bushy-II-I':
            # guinea pig data from Rothman and Manis, 2003, type II=I
            if debug:
                print 'Guinea pig bushy cell type II-I'
            if newModFiles:
                soma().kht2.gbar = nstomho(150.0, somaarea) * scalefactor
                soma().klt2.gbar = nstomho(35.0, somaarea) * scalefactor
            else:
                soma().kht.gbar = nstomho(150.0, somaarea) * scalefactor
                soma().klt.gbar = nstomho(35.0, somaarea) * scalefactor
            #s.ka.gkabar = nstomho(0.0, somaarea) * scalefactor
            soma().ihvcn.gbar = nstomho(3.5, somaarea) * scalefactor
            #s.hcno.gbar = 0.0
            soma().leak.gbar = nstomho(2.0, somaarea) * scalefactor
            vm0 = -63.8
        return vm0


# def bushy_addAxon(debug=False, ttx=False, message=None, nach=None):
# 
#     nnodes = range(5)
# #    print nnodes
#     axnode = []
#     internode = []
#     section = h.Section
#     initsegment = section(cell=soma)
#     initsegment.connect(soma)
#     #for i in range(nnodes):
#     #    axnode.append(Section(cell=soma))
#     #    internode.append(Section(cell=soma))
#     for i in nnodes:
#         axnode.append(Section(cell=soma))
#         internode.append(Section(cell=soma))
# 
#     axnode[0].connect(initsegment)
#     for i in nnodes:
#         internode[i].connect(axnode[i])
#         if i < nnodes[-1]:
#             axnode[i + 1].connect(internode[i])
# # create an initial segment
#     ninitseg = 21
#     initsegment.nseg = ninitseg
#     initsegment.diam = 4.0 * axonsf
#     initsegment.L = 36.0 * axonsf
#     initsegment.cm = cm
#     initsegment.Ra = Ra
#     initsegment.insert('na')
#     initsegment.insert('kht')
#     initsegment.insert('klt')
#     initsegment.insert('ihvcn')
#     initsegment.insert('leak')
#     gnamax = nstomho(6000.0, somaarea) * scalefactor
#     gnamin = 0.0 * gnamax
# 
#     gnastep = (gnamax - gnamin) / ninitseg
#     ip = 0
#     for inseg in initsegment:
#         ina = gnamin + ip * gnastep
#         print 'seg %d ina = %9.6f' % (ip, ina)
#         inseg.na.gbar = ina
#         inseg.klt.gbar = 0.2 * nstomho(200.0, somaarea) * scalefactor
#         inseg.kht.gbar = nstomho(150.0, somaarea) * scalefactor
#         inseg.ihvcn.gbar = 0.0 * nstomho(20.0, somaarea) * scalefactor
#         inseg.leak.gbar = nstomho(2.0, somaarea) * scalefactor
#         inseg.ena = v_sodium
#         inseg.ek = v_potassium
#         ip = ip + 1
# 
#     for i in nnodes:
#         axnode[i] = loadaxnodes(axnode[i], somaarea, scalefactor)
#         internode[i] = loadinternodes(internode[i], somaarea, scalefactor)
# 
#    # h.topology()
#     vm0 = -63.6
#     if not runQuiet:
#         if message is None:
#             print "<< bushy: JSR bushy cell model with axon created >>"
#         else:
#             print message
#     return(soma, [initsegment, axnode, internode])

class BushyWithAxon(Cell):
    def __init__(self, debug=False, ttx=False, message=None, nach=None):
        v_potassium = -70       # potassium reversal potential
        v_sodium = 55           # sodium reversal potential
        Ra = 150
        cm = 1.0

        species = 'mouse'
        if species is 'cat':
            axonsf = 1.0
            isdist = 0
        if species is 'mouse':
            axonsf = 0.57
            isdist = 0
        scalefactor = axonsf # This determines the relative size of the cell
    #    rinsf = 1.0           # input resistance adjustment (also current...)
        totcap = scalefactor * 12.0 # cap in pF for cell
        effcap = totcap # sometimes we change capacitance - that's effcap
        somaarea = totcap * 1E-6 / cm # pf -> uF, cm = 1uf/cm^2 nominal
        lstd = 1E4 * sqrt(somaarea / 3.14159) # convert from cm to um

        nnodes = range(5)
    #    print nnodes
        axnode = []
        internode = []
        Section = h.Section
        soma = Section() # one compartment of about 29000 um2
        initsegment = Section(cell=soma)
        initsegment.connect(soma)
        #for i in range(nnodes):
        #    axnode.append(Section(cell=soma))
        #    internode.append(Section(cell=soma))
        for i in nnodes:
            axnode.append(Section(cell=soma))
            internode.append(Section(cell=soma))

        axnode[0].connect(initsegment)
        for i in nnodes:
            internode[i].connect(axnode[i])
            if i < nnodes[-1]:
                axnode[i + 1].connect(internode[i])

        soma.nseg = 1
        soma.diam = lstd
        soma.L = lstd
        soma.Ra = Ra
        soma.cm = cm
        soma.insert('klt')
        soma.insert('kht')
        soma.insert('na')
        soma.insert('ihvcn')
        soma.insert('leak')
        soma.ena = v_sodium
        soma.ek = v_potassium
    #    print dir(soma)
        for s in soma:
            if ttx is False:
                s.na.gbar = 0. * nstomho(1000.0, somaarea) * scalefactor
            else:
                s.na.gbar = 0.0
            s.kht.gbar = nstomho(150.0, somaarea) * scalefactor
            s.klt.gbar = nstomho(200.0, somaarea) * scalefactor
            s.ihvcn.gbar = nstomho(20.0, somaarea) * scalefactor
            s.leak.gbar = nstomho(2.0, somaarea) * scalefactor

    # create an initial segment
        ninitseg = 21
        initsegment.nseg = ninitseg
        initsegment.diam = 4.0 * axonsf
        initsegment.L = 36.0 * axonsf
        initsegment.cm = cm
        initsegment.Ra = Ra
        initsegment.insert('na')
        initsegment.insert('kht')
        initsegment.insert('klt')
        initsegment.insert('ihvcn')
        initsegment.insert('leak')
        gnamax = nstomho(6000.0, somaarea) * scalefactor
        gnamin = 0.0 * gnamax

        gnastep = (gnamax - gnamin) / ninitseg
        ip = 0
        for inseg in initsegment:
            ina = gnamin + ip * gnastep
            print 'seg %d ina = %9.6f' % (ip, ina)
            inseg.na.gbar = ina
            inseg.klt.gbar = 0.2 * nstomho(200.0, somaarea) * scalefactor
            inseg.kht.gbar = nstomho(150.0, somaarea) * scalefactor
            inseg.ihvcn.gbar = 0.0 * nstomho(20.0, somaarea) * scalefactor
            inseg.leak.gbar = nstomho(2.0, somaarea) * scalefactor
            inseg.ena = v_sodium
            inseg.ek = v_potassium
            ip = ip + 1

        for i in nnodes:
            axnode[i] = loadaxnodes(axnode[i], somaarea, scalefactor)
            internode[i] = loadinternodes(internode[i], somaarea, scalefactor)

    # h.topology()
        vm0 = -63.6
        if not runQuiet:
            if message is None:
                print "<< bushy: JSR bushy cell model with axon created >>"
            else:
                print message
                
        self.soma = soma
        self.initsegment = initsegment
        self.axnode = axnode
        self.internode = internode


    def loadaxnodes(axnode, somaarea, scalefactor, nodeLength=2.5,
        nodeDiameter=2.0):
        v_potassium = -80       # potassium reversal potential
        v_sodium = 50           # sodium reversal potential
        Ra = 150
        cm = 1.0
        axnode.nseg = 1
        axnode.L = nodeLength * scalefactor
        axnode.diam = nodeDiameter * scalefactor
        axnode.Ra = Ra
        axnode.cm = cm
        axnode.insert('na')
        axnode.insert('kht')
        axnode.insert('klt')
        axnode.insert('leak')
        axnode.insert('ihvcn')
        for ax in axnode:
        #    print "   segment: ", ax
            ax.na.gbar = nstomho(1000.0, somaarea) * scalefactor
            ax.kht.gbar = nstomho(150.0, somaarea) * scalefactor
            ax.klt.gbar = nstomho(200.0, somaarea) * scalefactor
            ax.ihvcn.gbar = 0
            ax.leak.gbar = nstomho(2.0, somaarea) * scalefactor
            ax.ena = v_sodium
            ax.ek = v_potassium
        return(axnode)


    def loadinternodes(internode, somaarea, scalefactor, internodeLength=1000,
        internodeDiameter=10):
        v_potassium = -80       # potassium reversal potential
        v_sodium = 50           # sodium reversal potential
        Ra = 150
        cm = 1.0

        internode.nseg = 20
        internode.L = internodeLength * scalefactor
        internode.diam = internodeDiameter * scalefactor
        internode.Ra = Ra
        internode.cm = 0.002
        internode.insert('na')
        internode.insert('kht')
        internode.insert('leak')
        for inno in internode:
            inno.leak.gbar = nstomho(0.002, somaarea) * scalefactor
            inno.na.gbar = 0 * nstomho(500.0, somaarea) * scalefactor
            inno.kht.gbar = 0 * nstomho(150.0, somaarea) * scalefactor
            inno.ek = v_potassium
            inno.ena = v_sodium
            inno.leak.e = -80
        return(internode)


