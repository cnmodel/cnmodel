__author__ = 'pbmanis'

"""
Decorator:
A class to insert biophysical mechanisms into a model.
This function attempts to automatically decorate a hoc-imported model set of sections
with appropriate conductances.

The class takes as input the object hf, which is an instance of morphology 
It also takes the cellType, a string that directs how conductances should be inserted.

"""

import string
import numpy as np
import cnmodel.util as nu
from cnmodel.util import Params

class Decorator():
    def __init__(self, cell, parMap=None, verify=False):


        print 'cell type: ', cell.type
        cellType = cell.type.lower().capitalize()
        self.channelInfo = Params(newCm=1.0,
                              newRa=100.0,  # changed 10/20/2007 to center in range
                              newg_leak=0.000004935,
                              eK_def=-85, eNa_def=50,
                              ca_init=70e-6,  # free calcium in molar
                              v_init=-80,  # mV
                              pharmManip={'TTX': False, 'ZD': False, 'Cd': False, 'DTX': False, 'TEA': False,
                                          'XE': False},
                              cellType=cellType,
                              modelType=cell.status['modelType'],
                              distanceMap=cell.hr.distanceMap,
                              parMap=parMap,
        )
        self.excludeMechs = [] # ['ihvcn', 'kht', 'klt', 'nav11']

        cell.channel_manager(modelType=cell.status['modelType'])
#        print 'Cell: \n', dir(cell)
#        print 'mechanisms: ', cell.hr.mechanisms
        # gmapper allows us tor remap the names of mechanisms and their conductance names, which may
        # vary in the mod files.
        # The versions in the mechanisms directory here have been systematized, but this
        # dictionary may help when adding other conductances.

        self.gmapper = {'nacn': 'gbar', 'kht': 'gbar', 'klt': 'gbar', 'leak': 'gbar',
                        'ihvcn': 'gbar', 'jsrna': 'gbar', 'nav11': 'gbar', 'nacncoop': 'gbar',
                        'hcnobo': 'gbar'}
        self.erev_mapper = {'nacn': 'ena', 'kht': 'ek', 'klt': 'ek', 'leak': 'erev',
                        'ihvcn': 'eh', 'jsrna': 'ena', 'nav11': 'ena', 'nacncoop': 'ena',
                        'hcnobo': 'eh'}
        self._biophys(cell, verify=verify)
        print 'Decorator: Model Decorated with channels (if this appears more than once per cell, there is a problem)'


    def _biophys(self, cell, verify=False):
        """
        Inputs: run parameter structure, model parameter structure
        verify = True to run through the inserted mechanisms and see that they are really there.
        Outputs: None
        Action: Channel insertion into model
        Side Effects:
            Sets conductances in every different kind of section
            Does not update any class variables (via self).

        original hoc code: Paul B. Manis, Ph.D.
        25 Sept. 2007
        Modified to use gca for HH formulation of calcium current
        14 Oct 2007
        converted for Python, 17 Oct 2012 (PB Manis)
        modified to use new hf hoc_reader class to access section types and mechanisms 10-14 Feb 2014 pbmanis
        """
       # createFlag = False
        cellType = self.channelInfo.cellType
        parMap = self.channelInfo.parMap
        dmap = self.channelInfo.distanceMap
        if self.channelInfo is None:
            raise Exception('biophys - no parameters or info passed!')
        if verify:
            print('Biophys: Inserting channels as if cell type is {:s} with modelType {:s}'
                 .format(cellType, self.channelInfo.modelType))
        # print dir(hf)
        # print 'hf section groupkeys: ', hf.sec_groups.keys()
        # print 'hf section groups: ', hf.sec_groups
        
        cell.hr.mechanisms = []
        for s in cell.hr.sec_groups.keys():
            sectype = self.remapSectionType(string.rsplit(s, '[')[0])
            if sectype not in cell.channelMap.keys():
                print 'encountered unknown section group type: %s  Not decorating' % sectype
                continue
            # print 'Biophys: Section type: ', sectype, 'from: ', s
            # print sectype
            # print 'channel mapping keys: ', self.cMan.channelMap.keys()
            for mech in cell.channelMap[sectype].keys():
                if mech not in self.gmapper.keys():
                    print 'Mechanism %s not found? ' % mech
                    continue
                if mech in self.excludeMechs:
                    continue
                # if mech in ['ihvcn']:
                #     print 'mechanism %s excluded' % mech
                #     continue
                if verify:
                    print('Biophys: section group: {:s}  insert mechanism: {:s} at {:.8f}'
                        .format(s, mech, cell.channelMap[sectype][mech]))
                if mech not in cell.hr.mechanisms:
                    cell.hr.mechanisms.append(mech)
                x = nu.Mechanism(mech)
                for sec in cell.hr.sec_groups[s]:
                    x.insert_into(cell.hr.get_section(sec))
                    if verify:
                        print '   inserting into section', sec
                    gbar = self.gbarAdjust(cell, sectype, mech, sec)  # map density by location/distance
                #print 'parmap: ', parMap, mech
                setup = ('%s_%s' % (self.gmapper[mech], mech))
                if parMap is not None and mech in parMap.keys():  # note, this allows parmap to have elements BESIDES mechanisms
                    if verify:
                        print 'parMap[mech]', mech, parMap[mech], gbar,
                    gbar = gbar * parMap[mech]
                    if verify:
                        print '  new gbar: ', gbar
                for sec in cell.hr.sec_groups[s]:  # set cpmdictamces///
                    setattr(cell.hr.get_section(sec), setup, gbar)  # set conductance magnitude
                    if hasattr(cell, 'channelErevMap'):  # may not always have this mapping
                        secobj = cell.hr.get_section(sec)  # get the NEURON section object
                        mechsinsec = cell.get_mechs(secobj)  # get list of mechanisms in this section
                        if mech in mechsinsec:  # confirm that the mechanism is really there
                            setrev = False  # We try two ways for different mechanisms - just flag it
                            try:
                                setattr(secobj, self.erev_mapper[mech], cell.channelErevMap[sectype][mech])
                                setrev = True
                                continue  # don't bother with second approach
                            except:
                                pass  # no error
                            try:
                                setattr(secobj(), self.erev_mapper[mech] + '_' + mech, cell.channelErevMap[sectype][mech])
                                setrev = True
                            except:
                                pass  # no error report
                            if not setrev:  # here is our error report - soft, not crash.
                                print ("Failed to set reversal potential in section %s for mechanism %s" % (sec, mech))
        if verify:
            self.channelValidate(cell)
        return cell


    def gbarAdjust(self, cell, sectype, mech, sec):
        gbar = cell.channelMap[sectype][mech]
        gbar_orig = gbar
        if sectype not in cell.distMap.keys():  # no map for this section type
            return gbar
        elif mech not in cell.distMap[sectype].keys():
            return gbar
        # mecanism exists in the distMap, so we will map gbar to distance from soma
        method = cell.distMap[sectype][mech]['gradient'] # grab the type
        gminf = cell.distMap[sectype][mech]['gminf']
        rate = cell.distMap[sectype][mech]['lambda']
        if method == 'flat':
            return gbar
        if sec in self.channelInfo.distanceMap.keys():
            dist = self.channelInfo.distanceMap[sec]
        else:  # the sec should be in the map, but there could be a coding error that would break that relationship
            raise NameError('gbarAdjust in channel_decorate.py: section %s not in distance map' % sec)
        if method == 'linear':  # rate is "half" point drop
            gbar = gbar - dist*(gbar-gminf)/(2*rate)
            if gbar < 0.:
                gbar = 0. # clip
        elif method in ['exp', 'expdown']:
            gbar = (gbar - gminf) * np.exp(-dist/rate) + gminf
        if gbar < 0.:
            gbar = 0.
        #print 'gbaradjust: orig/adj: ', gbar_orig, gbar, method, dist, sectype
        return gbar

    def channelValidate(self, cell, verify=False):
        """
         verify mechanisms insertions -
         go through all the groups, and find inserted conductances and their values
         print the results to the terminal
        """
        print '\nChannel Validation'
        print 'looking for sec_groups: ', cell.hr.sec_groups.keys()
        print 'channelmaps: ', cell.channelMap.keys()
        secstuff = {}
        for s in cell.hr.sec_groups.keys():
            sectype = self.remapSectionType(string.rsplit(s, '[')[0])
            if sectype not in cell.channelMap.keys():
                if sectype in ['undefined']:  # skip undefined sections
                    continue
                print 'Validation: encountered unknown section group type: %s  Cannot Validate' % sectype
                continue
#            print 'Validating Section: %s' % s
            for mech in cell.channelMap[sectype].keys():
                if mech not in self.gmapper.keys():
                    continue
                if mech in self.excludeMechs:
                    continue
                if verify:
                    print '\tSection: %-15ss  found mechanism: %-8ss at %.5f' % (s, mech, cell.channelMap[sectype][mech])
                x = nu.Mechanism(mech) # , {gmapper[mech]: self.channelMap[cellType][sectype][mech]})
                setup = ('%s_%s' % (self.gmapper[mech], mech))
                for sec in cell.hr.sec_groups[s]:
                    bar = getattr(cell.hr.get_section(sec), setup)
                    # print 'mech', mech
                    # print 'bar: ', bar
                    try:
                        Erev = getattr(cell.hr.get_section(sec), self.erev_mapper[mech])
                    except:
                        Erev=-999.
                    # print 'erev: ', Erev
                    if sec in secstuff.keys():
                        secstuff[sec] += ', g_%s = %g [%.1f]' % (mech, bar, Erev)
                    else:
                        secstuff[sec] = '(%10s) g_%-6s = %g [%.1f] ' % (sectype, mech, bar, Erev)
        if verify:
            for i, k in enumerate(secstuff.keys()):
                print '**%-20s ' % k, secstuff[k]
                

    def remapSectionType(self, sectype):
        if sectype in ['AXON_0']:
            sectype = 'axon'
        if sectype in ['initseg', 'initialsegment']:
            sectype = 'initseg'
        if sectype in ['dendscaled_0', 'dendscaled_1', 'dendscaled_2', 'dendrite']:
            sectype = 'dend'
        if sectype in ['apical_dendrite']:
            sectype = 'apic'
        return sectype
