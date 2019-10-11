from __future__ import print_function
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

        cellType = cell.celltype.lower()
        self.channelInfo = Params(newCm=1.0,
                              newRa=150.0,  # standard value
                              newg_leak=0.000004935,
                              eK_def=-85, eNa_def=50,
                              ca_init=70e-6,  # free calcium in molar
                              v_init=-80,  # mV
                              pharmManip={'TTX': False, 'ZD': False, 'Cd': False, 'DTX': False, 'TEA': False,
                                          'XE': False},
                              cellType=cell.status['cellClass'],
                              modelType=cell.status['modelType'],
                              modelName=cell.status['modelName'],
                              distanceMap=cell.hr.distanceMap,
                              parMap=parMap,
        )
        self.excludeMechs = [] # ['ihvcn', 'kht', 'klt', 'nav11']
        # print('modelType in dec: ', cell.status['modelType'])
        # print('modelName in dec is ', cell.status['modelName'])
        # print('cell type in dec: ', cellType)
        cell.channel_manager(modelName=cell.status['modelName'], modelType=cell.status['modelType'])
#        print 'Cell: \n', dir(cell)
#        print 'mechanisms: ', cell.hr.mechanisms
        # gmapper allows us tor remap the names of mechanisms and their conductance names, which may
        # vary in the mod files.
        # The versions in the mechanisms directory here have been systematized, but this
        # dictionary may help when adding other conductances.

        self.gbar_mapper = {'nacn': 'gbar', 'kht': 'gbar', 'klt': 'gbar', 'leak': 'gbar',
                        'ihvcn': 'gbar', 'jsrna': 'gbar', 'nav11': 'gbar', 'nacncoop': 'gbar',
                        'nabu': 'gbar',
                        'hcnobo': 'gbar'}
        self.erev_mapper = {'nacn': 'ena', 'kht': 'ek', 'klt': 'ek', 'leak': 'erev', 'nabu': 'ena',
                        'ihvcn': 'eh', 'jsrna': 'ena', 'nav11': 'ena', 'nacncoop': 'ena',
                        'hcnobo': 'eh'}
        self.vshift_mapper = {'nacn': None, 'kht': None, 'klt': None, 'leak': None,
                        'ihvcn': None, 'jsrna': None, 'nav11': 'vsna', 'nacncoop': 'vsna', 'nabu': 'vshift',
                        'hcnobo': None}
        self._biophys(cell, verify=verify)
        print('\033[1;31;40m Decorator: Model Decorated with channels (if this appears more than once per cell, there is a problem)\033[0m')


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
        # check to see if we already did this
        # createFlag = False
        cellType = self.channelInfo.cellType
        parMap = self.channelInfo.parMap
        dmap = self.channelInfo.distanceMap
        if self.channelInfo is None:
            raise Exception('biophys - no parameters or info passed!')
        if verify:
            print(f'Biophys: Inserting channels as if cell type is {cellType:s} with modelType {self.channelInfo.modelType:s}')
        
        cell.hr.mechanisms = []
        for s in list(cell.hr.sec_groups.keys()):
            sectype = self.remapSectionType(s.rsplit('[')[0])
            if sectype not in cell.channelMap.keys():
                raise ValueError(f'Encountered unknown section group type: {sectype:s}. Cannot complete decoration')
            
            # here we go through all themechanisms in the ionchannels table for this cell and compartment type
            # note that a mechanism may have multiple parameters in the table (gbar, vshft), so we:
            #   a. only insert the mechanism once
            #   b. only adjust the relevant parameter
            # if sectype == 'soma':
            #     print('soma challenmap: ', cell.channelMap[sectype])
            #     exit()
            # print('sectype: ', sectype)
            for mechname in list(cell.channelMap[sectype].keys()):
                mech = mechname.split('_')[0] # get the part before the _
                parameter = mechname.split('_')[1]  # and the part after
                if mech in self.excludeMechs:
                    continue
                # print('    *** ', mech, parameter)

                if mech not in self.gbar_mapper.keys():
                    raise ValueError(f'Mechanism {mech:s} was requested in decorator but is not found?')
                if verify:
                    print(f'Biophys: section group: {s:s}  insert mechanism: {mech:s} at {cell.channelMap[sectype][mech]:.8f}')
                if mech not in cell.hr.mechanisms:  
                    cell.hr.mechanisms.append(mech)  # just add the mechanism to our list
                x = nu.Mechanism(mech)
                if cell.hr.sec_groups[s] == set():
                    continue  # no sections of this type
                for sec in cell.hr.sec_groups[s]:  # insert into all the sections of this type (group)
                    try:
                        x.insert_into(cell.hr.get_section(sec))
                    except:
                        raise ValueError(f'Failed to insert mechanism: {mech:s}')  # fail if cannot insert.
                    if verify:
                        print('   Successfully inserted mechanism {mech:s} into section: {str(sec}:s)} ')
                
                gbar_setup = None
                gbar = 0.
                if parameter == 'gbar':
                    gbar = self.gbarAdjust(cell, sectype, mechname, sec)  # map density by location/distance
                    gbar_setup = ('%s_%s' % (self.gbar_mapper[mech], mech))  # map name into .mod file name
                    # if parMap is not None and mech in parMap.keys():  # note, this allows parmap to have elements BESIDES mechanisms
                    #     if verify:
                    #         print 'parMap[mech]', mech, parMap[mech], gbar,
                    #     gbar = gbar * parMap[mech]  # change gbar here...
                    #     if verify:
                    # print(f'####### new gbar: ', gbar)
                
                vshift_setup = None
                vshift = 0.
                # print('mapper: ', self.vshift_mapper[mech])
                if self.vshift_mapper[mech] is not None:
                    vshift_setup = ('%s_%s' % (self.vshift_mapper[mech], mech))  # map voltage shift
                    vshift = cell.channelMap[sectype]['%s_%s' % (mech, self.vshift_mapper[mech])]
                    # print("Channel Map: \n   ", cell.channelMap)
                    # print(f'*********   Shift add to mechanism {mech:s}: gbar={gbar:e}, vshift: {vshift:.6f}')
                
                cell.hr.h.Ra = self.channelInfo.newRa
                for sec in cell.hr.sec_groups[s]:  # now set conductances and other parameters as requested
                    cell.hr.get_section(sec).Ra = self.channelInfo.newRa  # set Ra here
                    if gbar_setup is not None:
                        setattr(cell.hr.get_section(sec), gbar_setup, gbar)  # set conductance magnitude
                        # print('gbar_setup: %s %s' % (sectype, gbar_setup), gbar)
                    # if m is not 'None':
                    #     print('param, vshift_setup: ', parameter, vshift_setup)
                    #     print('mapper, mech, vshift, sectype: ', self.vshift_mapper[mech], mech, vshift, sectype)
                    # exit()
                    if vshift_setup is not None:
                        # print(cell.hr.get_section(sec))
                        try:
                            setattr(cell.hr.get_section(sec), vshift_setup, vshift)  # set shift magnitude
                        except:
                            print(dir(cell.hr.get_section(sec)))
                            raise ValueError (f'cannot set mechanism attribute %s  ... %s ' % (vshift_setup, vshift))
                        # print('\033[1;31;40m Vshift set to : \033[0m', vshift, vshift_setup)
                    
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
                                raise ValueError('erev set failed')
                                pass  # no error
                            try:
                                setattr(secobj(), self.erev_mapper[mech] + '_' + mech, cell.channelErevMap[sectype][mech])
                                setrev = True
                            except:
                                raise ValueError('Erev2 set failed')
                                pass  # no error report
                            if not setrev:  # here is our error report - soft, not crash.
                                print ("Failed to set reversal potential in section %s for mechanism %s" % (sec, mech))
                        # if mech in mechinsec and mech in self.vshift_mapper.keys():
                        #     try:
                        #         setattr(secobj, self.vshift_mapper[mech], cell.channelVshiftMap[sectype][mech])
                        #     except:
                        #         raise ValueError('Failed to set vshift for mech: %s  sectpe: %s' % (mech, sectype[mech]))
                    
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
      #  print('sectype: %s  mech: %s method: %s  rate: %s' % (sectype, mech, method, rate))
        if method == 'flat':
            return gbar
        if sec in self.channelInfo.distanceMap.keys():
            dist = self.channelInfo.distanceMap[sec]
        else:  # the sec should be in the map, but there could be a coding error that would break that relationship
            raise NameError('gbarAdjust in channel_decorate.py: section %s not in distance map' % sec)
        if method == 'linear':  # rate is "half" point drop
           # print('doing linear, orig gbar: ', gbar)
            gbar = gbar - dist*(gbar-gminf)/rate
            if gbar < 0.:
                gbar = 0. # clip
        #    print('sec dist: %f  final gbar: %f' % (dist, gbar))
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
        print('\nChannel Validation')
        print('  Looking for sec_groups: ', sorted(cell.hr.sec_groups.keys()))
        print('  Available Channel Maps: ', sorted(cell.channelMap.keys()))
        secstuff = {}
        for s in list(cell.hr.sec_groups.keys()):
            sectype = self.remapSectionType(s.rsplit('[')[0])
            if sectype not in cell.channelMap.keys():
                if sectype in ['undefined']:  # skip undefined sections
                    continue
                print('\033[1;31;40m Validation: encountered unknown section group type: %s  Cannot Validate' % sectype)
                print('Cell morphology file: %s \033[0m' % cell.morphology_file)
                continue
#            print 'Validating Section: %s' % s
            for mech in list(cell.channelMap[sectype].keys()):
                if mech not in self.gbar_mapper.keys():
                    continue
                if mech in self.excludeMechs:
                    continue
                if verify:
                    print('\tSection: %-15ss  found mechanism: %-8ss at %.5f' % (s, mech, cell.channelMap[sectype][mech]))
                x = nu.Mechanism(mech) # , {gmapper[mech]: self.channelMap[cellType][sectype][mech]})
                setup = ('%s_%s' % (self.gbar_mapper[mech], mech))
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
                print('**%-20s ' % k, secstuff[k])
                

    def remapSectionType(self, sectype):
        if sectype in ['AXON_0']:
            sectype = 'axon'

        if sectype in ['initseg', 'initialsegment']:
            sectype = 'initialsegment'
        if sectype in ['dendscaled_0', 'dendscaled_1', 'dendscaled_2', 'dendrite']:
            sectype = 'dendrite'
        if sectype in ['apical_dendrite']:
            sectype = 'secondarydendrite'
        return sectype
