import numpy as np
from neuron import h

from .psd import PSD


class GlyPSD(PSD):
    """Glycinergic PSD
    
    Parameters
    ----------
    section : Section
        The postsynaptic section in which to insert the receptor mechanism.
    terminal : Terminal
        The presynaptic Terminal instance
    params : dict
        Dictionary of kinetic parameters {'KV', 'KU', 'XMax'}
    psdType : str
        glyfast, glyslow, glyGC, glya5, or glyexp.
    """
    def __init__(self, section, terminal, params=None,
                 gmax=1000., psdType='glyfast',
                 message=None, debug=False,
                 gvar=0, eRev=-70):

        PSD.__init__(self, section, terminal)
        pre_sec = terminal.section
        post_sec = section
        
        from .. import cells
        
        params = {} if params is None else params
        
        self.pre_cell = cells.cell_from_section(pre_sec)
        self.post_cell = cells.cell_from_section(post_sec)

        self.psdType = psdType
        
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
            clefts.append(cl)
        
        # and then make a set of postsynaptic receptor mechanisms
        if self.psdType == 'glyslow':
            (psd, par) = self.template_Gly_PSD_State_Gly6S(nReceptors=n_rzones,
                                                        psdtype=self.psdType)
        elif self.psdType == 'glyfast':
            (psd, par) = self.template_Gly_PSD_State_PL(nReceptors=n_rzones,
                                                    psdtype=self.psdType)
        elif self.psdType == 'glyGC':
            (psd, par) = self.template_Gly_PSD_State_GC(nReceptors=n_rzones,
                                                    psdtype=self.psdType)
        elif self.psdType == 'glya5':
            (psd, par) = self.template_Gly_PSD_State_Glya5(nReceptors=n_rzones,
                                                        psdtype=self.psdType)
        elif self.psdType == 'glyexp':
            (psd, par) = self.template_Gly_PSD_exp(nReceptors=n_rzones,
                                                psdtype=self.psdType)
        else:
            print "**PSDTYPE IS NOT RECOGNIZED: [%s]\n" % (self.psdType)
            exit()
        if debug:
            print 'pre_sec: ', pre_sec
        
        
        # Connect terminal to psd (or cleft)
        self._cleft_netcons = []
        for k in range(0, n_rzones):
            pre_sec.push()
            netcon = h.NetCon(relzone._ref_XMTR[k], clefts[k], 0.1, 0.0, 1.0)
            self._cleft_netcons.append(netcon)
            h.pop_section()
            
            # set cleft transmitter kinetic parameters
            for pname, pval in params.items():
                setattr(clefts[k], pname, pval)

            clefts[k].setpointer(clefts[k]._ref_CXmtr, 'XMTR', psd[k]) #connect transmitter release to the PSD
            
            v = 1.0 + gvar * np.random.standard_normal()
            psd[k].gmax = gmax * v # add a little variability - gvar is CV of amplitudes
            #print 'GLY psd %s %d gmax=%f' % (self.psdType, k, gmax)
            psd[k].Erev = eRev # set the reversal potential
        
        par = list(par)
        #if self.psdType == 'ampa': # include nmda receptors in psd
            #psd.extend(psdn)
            #par.extend(parn)
        if message is not None:
            print message
                
        self.all_psd = psd
        self.clefts = clefts
        self.par = par


    # the following templates are a bit more complicated.
    # The parameter names as defined in the model are returned
    # but also those that are involved in the forward binding reactions are
    # listed separately - this is to allow us to run curve fits that adjust
    # only a subset of the parameters a  time - e.g., the rising phase, then
    # the falling phase with the rising phase fixed.
    # the dictionary selection is made by selectpars in glycine_fit.py.
    #

    def template_Gly_PSD_exp(self, debug=False, nReceptors=2, cellname=None, message=None):
        sec = self.section
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


    def template_Gly_PSD_State_Glya5(self, debug=False, nReceptors=2, psdtype=None, message=None):
        sec = self.section
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


    def template_Gly_PSD_State_Gly6S(self, debug=False, nReceptors=2, psdtype=None, message=None):
        sec = self.section
        psd = []
        sec.push()
        for k in range(0, nReceptors):
            psd.append(h.Gly6S(0.5, sec)) # simple using Trussell model 6 states with desens
            if debug:
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


    def template_Gly_PSD_State_PL(self, debug=False, nReceptors=2, cellname=None,
                                psdtype=None, message=None):
        sec = self.section
        psd = []
        sec.push()
        for k in range(0, nReceptors):
            psd.append(h.GLYaPL(0.5, sec)) # simple dextesche glycine receptors
            if debug:
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


    def template_Gly_PSD_State_GC(self, debug=False, nReceptors=2, psdtype=None, message=None):
        sec = self.section
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
