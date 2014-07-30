import numpy as np
from neuron import h


class PSD(object):
    """
    Base class for postsynaptic density mechanisms, possibly including cleft. 
    May accept either NetCon or pointer inputs from a Terminal, and directly
    modifies the membrane potential and/or ion concentrations of the 
    postsynaptic cell.    
    """
    def __init__(self, pre_sec, post_sec, terminal):
        pass
    
class GluPSD(PSD):
    """
    Glutamatergic PSD with ionotropic AMPA / NMDA receptors
    """
    def __init__(self, pre_sec, post_sec, terminal,
                 ampa_gmax,
                 nmda_ampa_ratio,
                 message=None, debug=False,
                 gvar=0, eRev=0,
                 nmda_ratio=1.0, identifier=0):
        """ This routine generates the synaptic connections from one presynaptic
            input onto a postsynaptic cell.
            Each connection is a stochastic presynaptic synapse ("presynaptic") with
            depression and facilitation.
            Each synapse can have  multiple releasesites "releasesites" on the
            target cell, as set by "NRZones".
            Each release site releases transmitter using "cleftXmtr"
            Each release site's transmitter is in turn attached to a PSD at each ending ("psd")
            Each psd can have a different conductance centered about the mean of
            gmax, according to a gaussian distribution set by gvar.
            
        Notes:
        
        *ampa_gmax* should be provided as the maximum *measured* AMPA conductance;
        this will be automatically corrected for the maximum open probability of
        the AMPA mechanism.
        
        *nmda_ampa_ratio* should be the ratio nmda/ampa Po measured at +40 mV.
        """
        from .. import cells
        self.AN_Po_Ratio = 23.2917 # ratio of open probabilities for AMPA and NMDAR's at peak currents
        self.AMPA_Max_Po = 0.44727
        self.NMDARatio = 0.0
        
        self.pre_cell = cells.cell_from_section(pre_sec)
        self.post_cell = cells.cell_from_section(post_sec)

        # get AMPA gmax corrected for max open probability
        gmax = ampa_gmax / self.AMPA_Max_Po
        
        relzone = terminal.relsite
        n_rzones = terminal.n_rzones
        
        #
        # Create cleft mechanisms
        # 
        clefts = []
        for k in range(0, n_rzones):
            cl = h.cleftXmtr(0.5, sec=post_sec)
            clefts.append(cl) # cleft
        
        # and then make a set of postsynaptic receptor mechanisms
        #        print 'PSDTYPE: ', psdtype
        (psd, psdn, par, parn) = template_iGluR_PSD(post_sec, nReceptors=n_rzones,
                                                    nmda_ratio=nmda_ratio)
        
        # Connect terminal to psd (or cleft)
        for k in range(0, n_rzones):
            # Note: cleft kinetics is implemented in the AMPA mechanism
            relzone.setpointer(relzone._ref_XMTR[k], 'XMTR', psd[k])
            
            # Note: NMDA has no cleft mechanism, but it has a slow response that
            # would not be strongly affected by the relatively fast cleft kinetics.
            relzone.setpointer(relzone._ref_XMTR[k], 'XMTR', psdn[k]) # include NMDAR's as well at same release site
            
            v = 1.0 + gvar * np.random.standard_normal()
            psd[k].gmax = gmax * v # add a little variability - gvar is CV of amplitudes
            psd[k].Erev = eRev # set the reversal potential
            
            # also adjust the nmda receptors at the same synapse
            psdn[k].gmax = gmax * v
            psdn[k].Erev = eRev
        
        par = list(par)
        psd.extend(psdn)
        par.extend(parn)
        if message is not None:
            print message
                
        self.psd = psd
        self.clefts = clefts
        self.par = par

        # adjust NMDA receptors to match postsynaptic cell
        k = 0
        kNMDA = -1
        kAMPA = -1
        for p in self.psd:
            if p.hname().find('NMDA', 0, 6) >= 0:
                p.gNAR = nmda_ampa_ratio * self.AN_Po_Ratio * self.NMDARatio
                p.vshift = 0
                if kNMDA == -1:
                    kNMDA = k # save the first instance where we have an NMDA receptor
            else:
                if kAMPA == -1: # not NMDA, so get AMPA 
                    kAMPA = k
            k = k + 1

        self.kNMDA = kNMDA
        self.kAMPA = kAMPA



class GlyPSD(PSD):
    """
    Glycinergic PSD
    """
    def __init__(self, pre_sec, post_sec, terminal,
                 ampa_gmax,
                 nmda_ampa_ratio,
                 message=None, debug=False,
                 gvar=0, eRev=0,
                 nmda_ratio=1.0, identifier=0):
        """ This routine generates the synaptic connections from one presynaptic
            input onto a postsynaptic cell.
            Each connection is a stochastic presynaptic synapse ("presynaptic") with
            depression and facilitation.
            Each synapse can have  multiple releasesites "releasesites" on the
            target cell, as set by "NRZones".
            Each release site releases transmitter using "cleftXmtr"
            Each release site's transmitter is in turn attached to a PSD at each ending ("psd")
            Each psd can have a different conductance centered about the mean of
            gmax, according to a gaussian distribution set by gvar.
            
            
        Notes:
        
        *ampa_gmax* should be provided as the maximum *measured* AMPA conductance;
        this will be automatically corrected for the maximum open probability of
        the AMPA mechanism.
        
        *nmda_ampa_ratio* should be the ratio nmda/ampa contuctance measured at +40 mV.
        """
        from .. import cells
        #self.AN_Po_Ratio = 23.2917 # ratio of open probabilities for AMPA and NMDAR's at peak currents
        #self.AMPA_Max_Po = 0.44727
        #self.NMDARatio = 0.0
        
        self.pre_cell = cells.cell_from_section(pre_sec)
        self.post_cell = cells.cell_from_section(post_sec)

        # get AMPA gmax corrected for max open probability
        #gmax = ampa_gmax / self.AMPA_Max_Po

        self.psdType = None
        
        #if cellname not in ['bushy', 'MNTB', 'stellate']:
            #raise TypeError
            #exit()
        #if debug:
            #print "\nTarget cell  = %s, psdtype = %s" % (cellname, psdtype)
        #printParams(stochastic_pars)
        glyslowPoMax = 0.162297  # thse were measured from the kinetic models in Synapses.py, as peak open P for the glycine receptors
        glyfastPoMax = 0.038475  # also later verified, same numbers...
        if self.psdType == 'glyfast':
            gmax /= glyfastPoMax  # normalized to maximum open probability for this receptor
        if self.psdType == 'glyslow':
            gmax /= glyslowPoMax  # normalized to max open prob for the slow receptor.
        
        #        print "Stochastic syn: j = %d of n_fibers = %d n_rzones = %d\n" % (j, n_fibers, n_rzones)
        relzone = terminal.relsite
        n_rzones = terminal.n_rzones
        
        #
        # Create cleft mechanisms
        # 
        clefts = []
        for k in range(0, n_rzones):
            cl = h.cleftXmtr(0.5, sec=post_sec)
            clefts.append(cl) # cleft
        
        # and then make a set of postsynaptic receptor mechanisms
        #        print 'PSDTYPE: ', psdtype
        if self.psdType == 'glyslow':
            (psd, par) = template_Gly_PSD_State_Gly6S(post_sec, nReceptors=n_rzones,
                                                        psdtype=self.psdType)
        elif self.psdType == 'glyfast':
            (psd, par) = template_Gly_PSD_State_PL(post_sec, nReceptors=n_rzones,
                                                    psdtype=self.psdType)
        elif self.psdType == 'glyGC':
            (psd, par) = template_Gly_PSD_State_GC(post_sec, nReceptors=n_rzones,
                                                    psdtype=self.psdType)
        elif self.psdType == 'glya5':
            (psd, par) = template_Gly_PSD_State_Glya5(post_sec, nReceptors=n_rzones,
                                                        psdtype=self.psdType)
        elif self.psdType == 'glyexp':
            (psd, par) = template_Gly_PSD_exp(post_sec, nReceptors=n_rzones,
                                                psdtype=self.psdType)
        else:
            print "**PSDTYPE IS NOT RECOGNIZED: [%s]\n" % (psdtype)
            exit()
        if debug:
            print 'pre_sec: ', pre_sec
        
        
        # Connect terminal to psd (or cleft)
        for k in range(0, n_rzones):
            if self.psdType in ['glyslow', 'glyfast', 'glyexp', 'glya5', 'glyGC']:
                relzone.setpointer(relzone._ref_XMTR[k], 'pre',
                                    clefts[k]) # connect the cleft to the release of transmitter
                clefts[k].preThresh = 0.1
                if isinstance(self.post_cell, cells.TStellate):
                    clefts[k].KV = 531.0 # set cleft transmitter kinetic parameters
                    clefts[k].KU = 4.17
                    clefts[k].XMax = 0.731
                elif isinstance(self.post_cell, cells.Bushy):
                    clefts[k].KV = 1e9 # really fast for bushy cells.
                    clefts[k].KU = 4.46
                    clefts[k].XMax = 0.733
                else:
                    print('Error in cell name, dying:   '), cellname
                    exit()
                clefts[k].setpointer(clefts[k]._ref_CXmtr, 'XMTR', psd[k]) #connect transmitter release to the PSD
            else:
                print "PSDTYPE IS NOT RECOGNIZED: [%s]\n" % (self.psdType)
                exit()
            v = 1.0 + gvar * np.random.standard_normal()
            psd[k].gmax = gmax * v # add a little variability - gvar is CV of amplitudes
            psd[k].Erev = eRev # set the reversal potential
        
        par = list(par)
        #if self.psdType == 'ampa': # include nmda receptors in psd
            #psd.extend(psdn)
            #par.extend(parn)
        if message is not None:
            print message
                
        self.psd = psd
        self.clefts = clefts
        self.par = par

        ## adjust NMDA receptors to match postsynaptic cell
        #k = 0
        #kNMDA = -1
        #kAMPA = -1
        #for p in self.psd:
            #if p.hname().find('NMDA', 0, 6) >= 0:
                #p.gNAR = nmda_ampa_ratio * self.AN_Po_Ratio * self.NMDARatio
                #p.vshift = 0
                #if kNMDA == -1:
                    #kNMDA = k # save the first instance where we have an NMDA receptor
            #else:
                #if kAMPA == -1: # not NMDA, so get AMPA 
                    #kAMPA = k
            #k = k + 1

        #self.kNMDA = kNMDA
        #self.kAMPA = kAMPA



# define the "Calyx" template.
# Each axon ends in (and thus includes) a calyx; therefore we make a list of axons
# based on the template.
# incoming argument: which axon (index) to address SGs
# Technically, this isn't the axon, but the terminal....

#def template_Calyx_Billup(debug=False, nzones=1, message=None):
    #calyx = cells.HH(message='  >> creating calyx_Billup')
    #coh = []
    #calyx.push()
    #for k in range(0, nzones):
        #coh.append(h.COH2(0.5, calyx))
    #h.pop_section()
    #return (calyx, coh)


def template_iGluR_PSD(sec, nReceptors=1, debug=False, cellname=None, message=None, nmda_ratio=1):
    """
    Create an ionotropic Glutamate receptor "PSD"
    Each PSD has receptors for each active zone, which must be matched (connected) to presynaptic
    terminals. Each PSD recetpor consists of an AMPATRUSSELL and an NMDA_KAMPA receptor
    Inputs:
        sec: The template requires a segment to insert the receptors into
        nReceptors: The number of receptor sites to insert
        debug: flag for debugging (prints extra information)
        cellname: Bushy/MNTB/stellate: determines ampa receptor kinetics
        message: Not used.
        nmda_ratio: The relative conductance of the open NMDA receptors to the open AMPA receptors.
    Outputs:
        (psd, psdn, par, parn)
        psd is the list of PSDs that were created (AMPA)
        psdn is the list of NMDA PSDs (same number as psd, just the NMDARs)
        par: dictionary of AMPAR kinetics as inserted
        parn: NMDA ratio
    Side Effecdts: None
    """
    psd = []
    psdn = []
    sec.push()
    AN_Po_Ratio = 23.2917 # ratio of open probabilities for AMPA and NMDAR's at peak currents
    for k in range(0, nReceptors):
        psd.append(h.AMPATRUSSELL(0.5, sec)) # raman/trussell AMPA with rectification
        psdn.append(h.NMDA_Kampa(0.5, sec)) # Kampa state model NMDA receptors

        if cellname in ['bushy', 'MNTB']:
            psd[-1].Ro1 = 107.85
            psd[-1].Ro2 = 0.6193
            psd[-1].Rc1 = 3.678
            psd[-1].Rc2 = 0.3212
            psdn[-1].gNAR = 0.036 * AN_Po_Ratio * nmda_ratio # 0.36*AN_Po_Ratio*nmda_ratio
            #if k == 0:
            #    print "Bushy NMDAR set to %8.2f" % psdn[-1].gNAR
        if cellname == 'stellate':
            psd[-1].Ro1 = 39.25
            psd[-1].Ro2 = 4.40
            psd[-1].Rc1 = 0.667
            psd[-1].Rc2 = 0.237
            psd[-1].PA = 0.1
            psdn[-1].gNAR = 1 * AN_Po_Ratio * nmda_ratio

    h.pop_section()
    par = {'Ro1': ('r', psd[0].Ro1),
           'Ro2': ('r', psd[0].Ro2),
           'Rc1': ('r', psd[0].Rc1),
           'Rc2': ('r', psd[0].Rc2), }
    parn = {'gNAR': ('4', psdn[0].gNAR), }
    return (psd, psdn, par, parn)


# the following templates are a bit more complicated.
# The parameter names as defined in the model are returned
# but also those that are involved in the forward binding reactions are
# listed separately - this is to allow us to run curve fits that adjust
# only a subset of the parameters a  time - e.g., the rising phase, then
# the falling phase with the rising phase fixed.
# the dictionary selection is made by selectpars in glycine_fit.py.
#

def template_Gly_PSD_exp(sec, debug=False, nReceptors=2, cellname=None, message=None):
    psd = []
    sec.push()
    for k in range(0, nReceptors):
        psd.append(h.GLY2(0.5, sec))
    h.pop_section()
    par = ['alpha', 'beta']
    p = []
    for n in par:
        p.append(eval('psd[0].' + n)) # get default values from the mod file
    return (psd, par, p)


def template_Gly_PSD_State_Glya5(sec, debug=False, nReceptors=2, psdtype=None, message=None):
    psd = []
    sec.push()
    for k in range(0, nReceptors):
        psd.append(h.GLYa5(0.5, sec))
    h.pop_section()
    par = {'kf1': ('r', psd[0].kf1), # retreive values in the MOD file
           'kf2': ('r', psd[0].kf2),
           'kb1': ('r', psd[0].kb1),
           'kb2': ('r', psd[0].kb2),
           'a1': ('f', psd[0].a1),
           'b1': ('f', psd[0].b1),
           'a2': ('f', psd[0].a2),
           'b2': ('f', psd[0].b2)}
    return (psd, par)


def template_Gly_PSD_State_Gly6S(sec, debug=False, nReceptors=2, psdtype=None, message=None):
    psd = []
    sec.push()
    for k in range(0, nReceptors):
        psd.append(h.Gly6S(0.5, sec)) # simple using Trussell model 6 states with desens
        if not runQuiet:
            print 'Gly6S psdtype: ', psdtype
        if psdtype == 'glyslow': # fit on 8 March 2010, error =  0.164, max open: 0.155
            psd[-1].Rd = 1.177999
            psd[-1].Rr = 0.000005
            psd[-1].Rb = 0.009403
            psd[-1].Ru2 = 0.000086
            psd[-1].Ro1 = 0.187858
            psd[-1].Ro2 = 1.064426
            psd[-1].Ru1 = 0.028696
            psd[-1].Rc1 = 0.103625
            psd[-1].Rc2 = 1.730578
    h.pop_section()
    par = {'Rb': ('n', psd[0].Rb), # retrive values in the MOD file
           'Ru1': ('r', psd[0].Ru1),
           'Ru2': ('r', psd[0].Ru2),
           'Rd': ('f', psd[0].Rd),
           'Rr': ('f', psd[0].Rr),
           'Ro1': ('f', psd[0].Ro1),
           'Ro2': ('r', psd[0].Ro2),
           'Rc1': ('f', psd[0].Rc1),
           'Rc2': ('r', psd[0].Rc2)}
    return (psd, par)


def template_Gly_PSD_State_PL(sec, debug=False, nReceptors=2, cellname=None,
                              psdtype=None, message=None):
    psd = []
    sec.push()
    for k in range(0, nReceptors):
        psd.append(h.GLYaPL(0.5, sec)) # simple dextesche glycine receptors
        if not runQuiet:
            print 'PL psdtype: ', psdtype
        if psdtype == 'glyslow':
            psd[-1].a1 = 0.000451
            psd[-1].a2 = 0.220
            psd[-1].b1 = 13.27
            psd[-1].b2 = 6.845
            psd[-1].kon = 0.00555
            psd[-1].koff = 2.256
            psd[-1].r = 1.060
            psd[-1].d = 55.03
        if psdtype == 'glyfast': # fit from 3/5/2010. error =  0.174 maxopen = 0.0385
            psd[-1].a1 = 1.000476
            psd[-1].a2 = 0.137903
            psd[-1].b1 = 1.700306
            psd[-1].koff = 13.143132
            psd[-1].kon = 0.038634
            psd[-1].r = 0.842504
            psd[-1].b2 = 8.051435
            psd[-1].d = 12.821820
    h.pop_section()
    par = {'kon': ('r', psd[0].kon), # retrive values in the MOD file
           'koff': ('r', psd[0].koff),
           'a1': ('r', psd[0].a1),
           'b1': ('r', psd[0].b1),
           'a2': ('f', psd[0].a2),
           'b2': ('f', psd[0].b2),
           'r': ('f', psd[0].r),
           'd': ('f', psd[0].d)}
    return (psd, par)


def template_Gly_PSD_State_GC(sec, debug=False, nReceptors=2, psdtype=None, message=None):
    psd = []
    sec.push()
    for k in range(0, nReceptors):
        psd.append(h.GLYaGC(0.5, sec)) # simple dextesche glycine receptors
        if psdtype == 'glyslow':
            psd[-1].k1 = 12.81    # (/uM /ms)   : binding
            psd[-1].km1 = 0.0087#   (/ms)   : unbinding
            psd[-1].a1 = 0.0195#    (/ms)   : opening
            psd[-1].b1 = 1.138    #(/ms)    : closing
            psd[-1].r1 = 6.13    #(/ms) : desense 1
            psd[-1].d1 = 0.000462 # (/ms)   : return from d1
            psd[-1].r2 = 0.731# (/ms)     : return from deep state
            psd[-1].d2 = 1.65 #(/ms)     : going to deep state
            psd[-1].r3 = 3.83 #(/ms)     : return from deep state
            psd[-1].d3 = 1.806 #(/ms)     : going to deep state
            psd[-1].rd = 1.04 #(/ms)
            psd[-1].dd = 1.004 # (/ms)
    h.pop_section()
    par = {'k1': ('r', psd[0].k1), # retrive values in the MOD file
           'km1': ('r', psd[0].km1),
           'a1': ('r', psd[0].a1),
           'b1': ('r', psd[0].b1),
           'r1': ('f', psd[0].r1),
           'd1': ('f', psd[0].d1),
           'r2': ('f', psd[0].r2),
           'd2': ('f', psd[0].d2),
           'r3': ('f', psd[0].r3),
           'd3': ('f', psd[0].d3),
           'rd': ('f', psd[0].rd),
           'dd': ('f', psd[0].dd)}
    return (psd, par)




