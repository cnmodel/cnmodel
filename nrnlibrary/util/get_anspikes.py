__author__ = 'pbmanis'
"""
ManageANSpikes is a class to read the output of the Zilany et al. 2009 AN model into
python, and provides services to access that data.

Basic usage is to create an instance of the class, and specify the data directory
if necessary.

You may then get the data in the format of a list using one of the "get" routines.
The data is pulled from the nearest CF or the specified CF.

"""
import os
import re
import numpy as np
from collections import OrderedDict
import scipy.io
import matplotlib.pyplot as MP


class ManageANSpikes():
    def __init__(self):
        # self.datadir = environment['HOME'] + '/Desktop/Matlab/ZilanyCarney-JASAcode-2009/'
        self.data_dir = os.environ['HOME'] + '/Desktop/Matlab/ANData/'
        self.data_read_flag = False
        self.dataType = 'RI'
        self.set_CF_map(4000, 38000, 25) # the default list
        self.all_AN = None

    def get_data_dir(self):
        return self.data_dir

    def set_data_dir(self, directory):
        if os.path.isdir(directory):
            self.data_dir = directory
        else:
            raise ValueError('ManageANSpikes.set_data_dir: Path %d is not a directory')

    def get_data_read(self):
        return self.data_read_flag

    def get_CF_map(self):
        return self.CF_map

    def set_CF_map(self, low, high, nfreq):
        self.CF_map = np.round(np.logspace(np.log10(low), np.log10(high), nfreq))

    def get_dataType(self):
        return self.dataType

    def set_dataType(self, datatype):
        if datatype in ['RI', 'PL']:  # only for recognized data types
            self.dataType = datatype
        else:
            raise ValueError ('get_anspikes.set_dataType: unrecognized type %s ' % datatype)

    def plot_RI_vs_F(self, freq=10000., spontclass='HS', display=False):
        cfd = {}
        if display:
            MP.figure(10)
        for i, fr in enumerate(self.CF_map):
            self.read_AN_data(freq=freq, CF=fr, spontclass=spontclass, ignoreflag=True)
            nsp = np.zeros(len(self.SPLs))  # this is the same regardless
            for i, db in enumerate(self.SPLs):
                spkl = self.combine_reps(self.spikelist[i])
                nsp[i] = len(spkl)#/self.SPLs[i]
            if display:
                MP.plot(self.SPLs, nsp)

    def plot_Rate_vs_F(self, freq=10000., spontclass='HS', display=False, SPL=0.):
        """
        assumes we have done a "read all"

        :param freq:
        :param spontclass:
        :param display:
        :param SPL:
        :return:
        """
        cfd = {}
        if display:
            MP.figure(10)
        for i, db in enumerate(self.SPLs):  # plot will have lines at each SPL
           # retrieve_from_all_AN(self, cf, SPL)   # self.read_AN_data(freq=freq, CF=fr, spontclass=spontclass)
            nsp = np.zeros(len(self.CF_map))  # this is the same regardless
            for j, fr in enumerate(self.CF_map):
                spikelist = self.retrieve_from_all_AN(fr, db)
                spkl = self.combine_reps(spikelist)
                nsp[j] = len(spkl)/len(spikelist)
            if display:
                MP.plot(self.CF_map, nsp)

    def getANatFandSPL(self, spontclass='MS', freq=10000., CF=None, SPL=0):
        """
        getANatFandSPL:
        ====================
        Get the AN data at a particular frequency and SPL. Note that the tone freq is specified,
        but the CF of the fiber might be different. If CF is None, we try to get the closest
        fiber data to the tone Freq.
        Other arguments are the spontaneous rate and the SPL level.
        The data must exist, or we epically fail.

        :param spontclass: 'HS', 'MS', or 'LS' for high, middle or low spont groups (1,2,3 in Zilany et al model)
        :param freq: The stimulus frequency, in Hz
        :param CF:  The 'characteristic frequency' of the desired AN fiber, in Hz
        :param SPL: The sound pressure level, in dB SPL for the stimulus
        :return: an array of nReps of spike times.
        """
        if CF is None:
            closest = np.argmin(np.abs(self.CF_map-freq))  # find closest to the stim freq
        else:
            closest = np.argmin(np.abs(self.CF_map-CF))  # find closest to the stim freq
        CF = self.CF_map[closest]
        print 'closest f: ', CF
        self.read_AN_data(spontclass=spontclass, freq=10000., CF=CF)
        return(self.get_AN_at_SPL(SPL))

    def get_AN_at_SPL(self, spl=None):
        """
        get_AN_at_SPL:
        ====================
        grabs the AN data for the requested SPL for the data set currently loaded
        The data must exist, or we epically fail.

        :param spl: sound pressure level, in dB SPL
        :return: spike trains (times) as a list of nReps numpy arrays.
        """
        if not self.data_read_flag:
            print 'getANatSPL: No data read yet'
            return None
        try:
            k = int(np.where(self.SPLs == spl)[0])
        except (ValueError, TypeError) as e:
            print 'get_anspikes::getANatSPL: spl=%6.1f not in list of spls:' % (spl), self.SPLs
            exit()  # no match
        # now clean up array to keep it iterable upon return
        spkl = [[]]*len(self.spikelist[k])
        for i in xrange(len(self.spikelist[k])):
            try:
                len(self.spikelist[k][i])
                spkl[i] = self.spikelist[k][i]
            except:
                spkl[i] = [np.array(self.spikelist[k][i])]
        return spkl  # returns all the trials in the list

    def read_all_ANdata(self, freq=10000., CFList=None, spontclass='HS', stim='BFTone'):
        """
        read_all_ANdata:
        ================
        Reads a bank of AN data, across frequency and intensity, for a given stimulus frequency
        Assumptions: the bank of data is consistent in terms of nReps and SPLs

        :param freq: Tone stimulus frequency (if stim is Tone or BFTone)
        :param CFList: A list of the CFs to read
        :param spontclass: 'HS', 'MS', or 'LS', corresponding to spont rate groups
        :param stim: Stimulus type - Tone/BFTone, or Noise
        :return: Nothing. The data are stored in an array accessible directly or through a selector function
        """
        self.all_AN = OrderedDict([(f, None) for f in CFList])  # access through dictionary with keys as freq
        for cf in CFList:
            self.read_AN_data(freq=freq, CF=cf, spontclass=spontclass, stim=stim, setflag=False)
            self.all_AN[cf] = self.spikelist
        self.data_read_flag = True  # inform and block

    def retrieve_from_all_AN(self, cf, SPL):
        """
        retrieve_from_all_AN:
        ====================+

        :param cf:
        :param SPL:
        :return: spike list (all nreps) at the cf and SPL requested, for the loaded stimulus set

        """
        if not self.data_read_flag or self.all_AN is None:
            print 'get_anspikes::retrieve_from_all_AN: No data read yet: '
            exit()

        ispl = int(np.where(self.SPLs == float(SPL))[0])
        icf = self.all_AN.keys().index(cf)
        spikelist = self.all_AN[cf][ispl]
        spkl = [[]]*len(spikelist)
        #print len(spikelist)
        for i in xrange(len(spikelist)): # across trials
            try:
                len(spikelist[i])
                spkl[i] = spikelist[i]
            except:
                spkl[i] = [np.array(spikelist[i])]
        #print 'cf: %f  spl: %d  nspk: %f' % (cf, SPL, len(spkl[i]))
        return(spikelist)

    def read_AN_data(self, freq=10000, CF=5300, spontclass = 'HS', display=False,
                     stim='BFTone', setflag=True, ignoreflag = True):
        """
        read responses of auditory nerve model of Zilany et al. (2009).
        display = True plots the psth's for all ANF channels in the dataset
        Request response to stimulus at Freq, for fiber with CF
        and specified spont rate
        This version is for rate-intensity runs March 2014.
        Returns:

        """
        if not ignoreflag:
            assert self.data_read_flag == False
            # print 'Each instance of ManageANSPikes is allowed ONE data set to manage'
            # print 'This is a design decision to avoid confusion about which data is in the instance'
            # exit()
        if stim in ['BFTone', 'Tone']:
            fname = '%s_F%06.3f_CF%06.3f_%2s.mat' % (self.dataType, freq/1000.0, CF/1000., spontclass)
        elif stim=='Noise':
            fname = '%s_Noise_CF%06.3f_%2s.mat' % (self.dataType, CF/1000., spontclass)

        #print 'Reading: %s' % (fname)
        try:
            mfile = scipy.io.loadmat(os.path.join(self.data_dir, fname), squeeze_me = True)
        except IOError:
            print 'get_anspikes::read_AN_data: Failed to find data file %s' % (fname)
            print 'Corresponding to Freq: %f  CF: %f  spontaneous rate class: %s' % (freq, CF, spontclass)
            exit()

        n_spl = len(mfile[self.dataType])
        if display:
            n = int(np.sqrt(n_spl))+1
            mg = np.meshgrid(range(n), range(n))
            mg = [mg[0].flatten(), mg[1].flatten()]
        spkl = [[]]*n_spl
        spl = np.zeros(n_spl)
        for k in range(n_spl):
            spkl[k] = mfile[self.dataType]['data'][k] # get data for one SPL
            spl[k] = mfile[self.dataType]['SPL'][k]  # get SPL for this set of runs
            #print spkl[k].shape
            if display:
                self.display(spkl[k], k, n,  mg)

        if display:
            MP.show()
        self.spikelist = spkl  # save these so we can have a single point to parse them
        self.SPLs = spl
        self.n_reps = mfile[self.dataType]['nrep']
        if setflag:
            self.data_read_flag = True
        #return spkl, spl, mfile['RI']['nrep']

    def getANatSPL(self, spl=None):
        if not self.dataRead:
            print 'getANatSPL: No data read yet'
            return None
        try:
            k = int(np.where(self.SPLs == spl)[0])
        except (ValueError, TypeError) as e:
            print 'get_anspikes::getANatSPL: spl=%6.1f not in list of spls:' % (spl), self.SPLs
            exit()  # no match
        # now clean up array to keep it iterable upon return
        spkl = [[]]*len(self.spikelist[k])
        for i in xrange(len(self.spikelist[k])):
            try:
                len(self.spikelist[k][i])
                spkl[i] = self.spikelist[k][i]
            except:
                spkl[i] = [np.array(self.spikelist[k][i])]
        return spkl  # returns all the trials in the list

    def get_AN_info(self):
        """
        get_AN_info:
        ============
        :return: Dictionary of nReps and SPLs that are in the current data set (instance).
        """
        if not self.data_read_flag:
            print 'getANatSPL: No data read yet'
            return None
        return {'nReps': self.n_reps, 'SPLs': self.SPLs}

    def read_cmrr_data(self, P, signal = 'S0', display=False):
        """
        read responses of auditory nerve model of Zilany et al. (2009).
        The required parameters are passed via the class P.
        This includes:
            the SN (s2m) is the signal-to-masker ratio to select
            the mode is 'CMR', 'CMD' (deviant) or 'REF' (reference, no flanking bands),
            'Tone', or 'Noise'
            the Spont rate group
        signal is either 'S0' (signal present) or 'NS' (no signal)
        display = True plots the psth's for all ANF channels in the dataset
        Returns:
        tuple of (spikelist and frequency list for modulated tones)
        """
        if P.modF == 10:
            datablock = '%s_F4000.0' % (P.mode)
        else:
            datablock = '%s_F4000.0_M%06.1f' % (P.mode, P.modF) # 100 Hz modulation directory
        #print P.mode
        fs=os.listdir(self.datadir+datablock)
        s_sn = "SN%03d" % (P.s2m)
        fntemplate = '(\S+)_%s_%s_%s_%s.mat' % (s_sn, P.mode, signal, P.SR)
        p = re.compile(fntemplate)
        fl = [re.match(p, file) for file in fs]
        fl = [f.group(0) for f in fl if f != None] # returns list of files matching...
        fr = [float(f[1:7]) for f in fl] # frequency list
        if display:
            self.makeFig()
        i = 0
        j = 0
        spkl = [[]]*len(fl)
        n = int(np.sqrt(len(fl)))
        mg = np.meshgrid(range(n), range(n))
        mg = [mg[0].flatten(), mg[1].flatten()]
        for k, fi in enumerate(fl):
            mfile = scipy.io.loadmat(self.datadir + datablock + '/' + fi, squeeze_me = True)
            spkl[k] = mfile['Cpsth']
            if display:
                self.display(spkl[k], k, n,  mg)

        if display:
            MP.show()
        return(spkl, fr)

    def make_fig(self):
        self.fig = MP.figure(1)

    def combine_reps(self, spkl):
        """
        combine_reps:
        =============
        Just turns the spike list into one big linear sequence.
        :param spkl: a spike train (nreps of numpy arrays)
        :return: all the spikes in one long array.
        """
         #print spkl
        allsp = np.array(self.flatten(spkl))
        #print allsp.shape
        if allsp.shape == (): # nothing to show
            return
        #print allsp
        spks = []
        for i in range(len(allsp)):
            if allsp[i] == ():
               continue
            if isinstance(allsp[i], float):
                spks.append(allsp[i])
            else:
                spks.extend(np.array(allsp[i]))
        return spks

    def display(self, spkl, k, n, mg):
        #print spkl
        spks = self.combine_reps(spkl)
        if spks != []:
            spks = np.sort(spks)
            MP.hist(spks, 100)
        return

    def flatten(self, x):
        """flatten(sequence) -> list

        Returns a single, flat list which contains all elements retrieved
        from the sequence and all recursively contained sub-sequences
        (iterables).
        """
        result = []
        if isinstance(x, float) or len([x]) <= 1:
            return x
        for el in x:
            if hasattr(el, "__iter__") and not isinstance(el, basestring):
                result.extend(self.flatten(el))
            else:
                result.append(el)
        return result

if __name__ == '__main__':
    import pylibrary.Params as Params

    modes = ['CMR', 'CMD', 'REF']
    sr_types = ['H', 'M', 'L']
    sr_type = sr_types[0]
    sr_names = {'H': 'HS', 'M': 'MS', 'L': 'LS'} # spontaneous rate groups (AN fiber input selection)
    # define params for reading/testing. This is a subset of params...
    P = Params.Params(mode = modes[0], s2m = 0, n_rep=0, start_rep=0,
                      sr = sr_names['L'],
                      modulation_freq = 10, DS_phase = 0, dataset='test',
                      sr_types=sr_types, sr_names=sr_names, fileVersion=1.5,
                     )

    manager = ManageANSpikes()  # create instance of the manager
    test='RI'

    if test == 'all':
        manager.read_all_ANdata(freq=10000., CFList=manager.CF_map, spontclass='HS', stim='BFTone')
        #spkl = manager.retrieve_from_all_AN(manager.CF_map[0], 40.)
        #spks = manager.combine_reps(spkl)
        manager.plot_Rate_vs_F(freq=manager.CF_map[0], display=True, spontclass='HS')
        MP.show()

    if test == manager.dataType:
        spikes = manager.get_AN_at_F_and_SPL(spontclass = 'MS', freq=10000., CF=None, SPL=50)
        spks = manager.combine_reps(spikes)
        if spks != []:
            print 'spikes.'
            manager.plot_RI_vs_F(freq=10000., display=True, spontclass='MS')
            MP.figure(11)
            spks = np.sort(spks)
            MP.hist(spks, 100)
            MP.show()
