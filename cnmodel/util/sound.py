"""
Tools for generating auditory stimuli. 
"""
from __future__ import division
import numpy as np
import scipy
import scipy.io.wavfile
import resampy

def create(type, **kwds):
    """ Create a Sound instance using a key returned by Sound.key().
    """
    cls = globals()[type]
    return cls(**kwds)


class Sound(object):
    """
    Base class for all sound stimulus generators.
    """
    def __init__(self, duration, rate=100e3, **kwds):
        """
        Parameters
        ----------
        duration: float (no default):
            duration of the stimulus, in seconds
        
        rate : float (default: 100000.)
            sample rate for sound generation
        
        """
        self.opts = {'rate': rate, 'duration': duration}
        self.opts.update(kwds)
        self._time = None
        self._sound = None

    @property
    def sound(self):
        """ 
        :obj:`array`: The generated sound array, expressed in Pascals.
        """
        if self._sound is None:
            self._sound = self.generate()
        return self._sound
        
    @property
    def time(self):
        """
        :obj:`array`: The array of time values, expressed in seconds.
        """
        if self._time is None:
            self._time = np.linspace(0, self.opts['duration'], self.num_samples)
        return self._time

    @property 
    def num_samples(self):
        """ 
        int: The number of samples in the sound array. 
        """
        return 1 + int(self.opts['duration'] * self.opts['rate'])
    
    @property
    def dt(self):
        """
        float: the sample period (time step between samples).
        """
        return 1.0 / self.opts['rate']
    
    @property
    def duration(self):
        """
        float: The duration of the sound
        """
        return self.opts['duration']

    def key(self):
        """ 
        The sound can be recreated using ``create(**key)``.
        :obj:`dict`: Return dict of parameters needed to completely describe this sound.
        """
        k = self.opts.copy()
        k['type'] = self.__class__.__name__
        return k

    def measure_dbspl(self, tstart, tend):
        """ 
        Measure the sound pressure for the waveform in a window of time
        
        Parameters
        ----------
        tstart :
            time to start spl measurement (seconds).
        
        tend :
            ending time for spl measurement (seconds).
        
        Returns
        -------
        float : The measured amplitude (dBSPL) of the sound from tstart to tend

        """
        istart = int(tstart * self.opts['rate'])
        iend = int(tend * self.opts['rate'])
        return pa_to_dbspl(self.sound[istart:iend].std())

    def generate(self):
        """
        Generate and return the sound output. This method is defined by subclasses.
        """
        raise NotImplementedError()

    def __getattr__(self, name):
        if 'opts' not in self.__dict__:
            raise AttributeError(name)
        if name in self.opts:
            return self.opts[name]
        else:
            return object.__getattr__(self, name)


class TonePip(Sound):
    """ Create one or more tone pips with cosine-ramped edges.
    
    Parameters
    ----------
    rate : float
        Sample rate in Hz
    duration : float
        Total duration of the sound
    f0 : float or array-like
        Tone frequency in Hz. Must be less than half of the sample rate.
    dbspl : float
        Maximum amplitude of tone in dB SPL. 
    pip_duration : float
        Duration of each pip including ramp time. Must be at least 
        2 * ramp_duration.
    pip_start : array-like
        Start times of each pip
    ramp_duration : float
        Duration of a single ramp period (from minimum to maximum). 
        This may not be more than half of pip_duration.
    
    """
    def __init__(self, **kwds):
        reqdWords = ['rate', 'duration', 'f0', 'dbspl', 'pip_duration', 'pip_start', 'ramp_duration']
        for k in reqdWords:
            if k not in kwds.keys():
                raise TypeError("Missing required argument '%s'" % k)
        if kwds['pip_duration'] < kwds['ramp_duration'] * 2:
            raise ValueError("pip_duration must be greater than (2 * ramp_duration).")
        if kwds['f0'] > kwds['rate'] * 0.5:
            raise ValueError("f0 must be less than (0.5 * rate).")
        Sound.__init__(self, **kwds)
        
    def generate(self):
        """
        Call to compute the tone pips
        
        Returns
        -------
        array : 
            generated waveform
        
        """
        o = self.opts
        return piptone(self.time, o['ramp_duration'], o['rate'], o['f0'], 
                       o['dbspl'], o['pip_duration'], o['pip_start'])


class FMSweep(Sound):
    """ Create an FM sweep with either linear or logarithmic rates, 
    of a specified duration between two frequencies.
    
    Parameters
    ----------
    rate : float
        Sample rate in Hz
    duration : float
        Total duration of the sweep
    start : float
        t times of each pip
    freqs : list 
        [f0, f1]: the start and stop times for the sweep
    ramp : string
        valid input for type of sweep (linear, logarithmic, etc)
    dbspl : float
        Maximum amplitude of pip in dB SPL.
    """
    def __init__(self, **kwds):
        for k in ['rate', 'duration', 'start', 'freqs', 'ramp', 'dbspl']:
            if k not in kwds:
                raise TypeError("Missing required argument '%s'" % k)
        Sound.__init__(self, **kwds)
        
    def generate(self):
        """
        Call to actually compute the the FM sweep
                
        Returns
        -------
        array : 
            generated waveform
        
        """
        o = self.opts
        return fmsweep(self.time, o['start'], o['duration'],
                         o['freqs'], o['ramp'], o['dbspl'])

class NoisePip(Sound):
    """ One or more noise pips with cosine-ramped edges.
    
    Parameters
    ----------
    rate : float
        Sample rate in Hz
    duration : float
        Total duration of the sound
    seed : int >= 0
        Random seed
    dbspl : float
        Maximum amplitude of tone in dB SPL. 
    pip_duration : float
        Duration of each pip including ramp time. Must be at least 
        2 * ramp_duration.
    pip_start : array-like
        Start times of each pip
    ramp_duration : float
        Duration of a single ramp period (from minimum to maximum). 
        This may not be more than half of pip_duration.
        
    """
    def __init__(self, **kwds):
        for k in ['rate', 'duration', 'dbspl', 'pip_duration', 'pip_start', 'ramp_duration', 'seed']:
            if k not in kwds:
                raise TypeError("Missing required argument '%s'" % k)
        if kwds['pip_duration'] < kwds['ramp_duration'] * 2:
            raise ValueError("pip_duration must be greater than (2 * ramp_duration).")
        if kwds['seed'] < 0:
            raise ValueError("Random seed must be integer > 0")
        
        Sound.__init__(self, **kwds)
        
    def generate(self):
        """
        Call to compute the noise pips
        
        Returns
        -------
        array : 
            generated waveform
        
        """
        o = self.opts
        return pipnoise(self.time, o['ramp_duration'], o['rate'],
                        o['dbspl'], o['pip_duration'], o['pip_start'], o['seed'])

class ClickTrain(Sound):
    """ One or more clicks (rectangular pulses).
    
    Parameters
    ----------
    rate : float
        Sample rate in Hz
    dbspl : float
        Maximum amplitude of click in dB SPL. 
    click_duration : float
        Duration of each click including ramp time. Must be at least 
        1/rate.
    click_starts : array-like
        Start times of each click
    """
    def __init__(self, **kwds):
        for k in ['rate',  'duration', 'dbspl', 'click_duration', 'click_starts']:
            if k not in kwds:
                raise TypeError("Missing required argument '%s'" % k)
        if kwds['click_duration'] < 1./kwds['rate']:
            raise ValueError("click_duration must be greater than sample rate.")
        if len(kwds['click_starts']) < 1:
            raise ValueError("Click Train needs at least one click.")
            
        Sound.__init__(self, **kwds)
        
    def generate(self):
        o = self.opts
        out = []
        print('generate click')
        print(o['click_starts'])
        for i, start in enumerate(o['click_starts']):
            print(i, start)
            if i == 0:
                out = click(self.time, o['rate'], 
                        o['dbspl'], o['click_duration'], start)
            else:
                out += click(self.time, o['rate'], 
                        o['dbspl'], o['click_duration'], start)
        print('click out: ', out)
        return out
    

class SAMNoise(Sound):
    """ One or more gaussian noise pips with cosine-ramped edges.
    
    Parameters
    ----------
    rate : float
        Sample rate in Hz
    duration : float
        Total duration of the sound
    seed : int >= 0
        Random seed
    dbspl : float
        Maximum amplitude of pip in dB SPL. 
    pip_duration : float
        Duration of each pip including ramp time. Must be at least 
        2 * ramp_duration.
    pip_start : array-like
        Start times of each pip
    ramp_duration : float
        Duration of a single ramp period (from minimum to maximum). 
        This may not be more than half of pip_duration.
    fmod : float
        SAM modulation frequency
    dmod : float
        Modulation depth
    """
    def __init__(self, **kwds):
        parms = ['rate', 'duration', 'seed', 'pip_duration', 
                 'pip_start', 'ramp_duration', 'fmod', 'dmod', 'seed']
        for k in parms:
            if k not in kwds:
                raise TypeError("Missing required argument '%s'" % k)
        if kwds['pip_duration'] < kwds['ramp_duration'] * 2:
            raise ValueError("pip_duration must be greater than (2 * ramp_duration).")
        if kwds['seed'] < 0:
            raise ValueError("Random seed must be integer > 0")
        
        Sound.__init__(self, **kwds)
        
    def generate(self):
        """
        Call to compute the SAM noise
        
        Returns
        -------
        array : 
            generated waveform
        
        """
        o = self.opts
        o['phaseshift'] = 0.
        return modnoise(self.time, o['ramp_duration'], o['rate'], o['f0'], 
                       o['pip_duration'], o['pip_start'], o['dbspl'],
                       o['fmod'], o['dmod'], 0., o['seed'])


class SAMTone(Sound):
    """ SAM tones with cosine-ramped edges.
    
    Parameters
    ----------
    rate : float
        Sample rate in Hz
    duration : float
        Total duration of the sound
    f0 : float or array-like
        Tone frequency in Hz. Must be less than half of the sample rate.
    dbspl : float
        Maximum amplitude of tone in dB SPL. 
    pip_duration : float
        Duration of each pip including ramp time. Must be at least 
        2 * ramp_duration.
    pip_start : array-like
        Start times of each pip
    ramp_duration : float
        Duration of a single ramp period (from minimum to maximum). 
        This may not be more than half of pip_duration.
    fmod : float
        SAM modulation frequency, Hz
    dmod : float
        Modulation depth, %
        
    """
    def __init__(self, **kwds):

        for k in ['rate', 'duration', 'f0', 'dbspl', 'pip_duration', 'pip_start',
                  'ramp_duration', 'fmod', 'dmod']:
            if k not in kwds:
                raise TypeError("Missing required argument '%s'" % k)
        if kwds['pip_duration'] < kwds['ramp_duration'] * 2:
            raise ValueError("pip_duration must be greater than (2 * ramp_duration).")
        if kwds['f0'] > kwds['rate'] * 0.5:
            raise ValueError("f0 must be less than (0.5 * rate).")
        
        Sound.__init__(self, **kwds)
        
    def generate(self):
        """
        Call to compute a SAM tone
        
        Returns
        -------
        array : 
            generated waveform
        
        """
        o = self.opts
        basetone = piptone(self.time, o['ramp_duration'], o['rate'], o['f0'], 
                       o['dbspl'], o['pip_duration'], o['pip_start'])
        return sinusoidal_modulation(self.time, basetone, o['pip_start'], o['fmod'], o['dmod'], 0.)


def pa_to_dbspl(pa, ref=20e-6):
    """ Convert Pascals (rms) to dBSPL. By default, the reference pressure is
    20 uPa.
    """
    return 20 * np.log10(pa / ref)


def dbspl_to_pa(dbspl, ref=20e-6):
    """ Convert dBSPL to Pascals (rms). By default, the reference pressure is
    20 uPa.
    """
    return ref * 10**(dbspl / 20.0)


class SAMNoise(Sound):
    """ One or more gaussian noise pips with cosine-ramped edges, sinusoidally modulated.
    
    Parameters
    ----------
    rate : float
        Sample rate in Hz
    duration : float
        Total duration of the sound
    seed : int >= 0
        Random seed
    dbspl : float
        Maximum amplitude of pip in dB SPL. 
    pip_duration : float
        Duration of each pip including ramp time. Must be at least 
        2 * ramp_duration.
    pip_start : array-like
        Start times of each pip
    ramp_duration : float
        Duration of a single ramp period (from minimum to maximum). 
        This may not be more than half of pip_duration.
    fmod : float
        SAM modulation frequency
    dmod : float
        Modulation depth
    
    Returns
    -------
    array :
        waveform
    
    """
    def __init__(self, **kwds):
        for k in ['rate', 'duration', 'seed', 'pip_duration', 'pip_start', 'ramp_duration', 
                    'fmod', 'dmod']:
            if k not in kwds:
                raise TypeError("Missing required argument '%s'" % k)
        if kwds['pip_duration'] < kwds['ramp_duration'] * 2:
            raise ValueError("pip_duration must be greater than (2 * ramp_duration).")
        if kwds['seed'] < 0:
            raise ValueError("Random seed must be integer > 0")
        
        Sound.__init__(self, **kwds)
        
    def generate(self):
        o = self.opts
        basenoise = pipnoise(self.time, o['ramp_duration'], o['rate'],
                        o['dbspl'], o['pip_duration'], o['pip_start'], o['seed'])
        return sinusoidal_modulation(self.time, basenoise, o['pip_start'], o['fmod'], o['dmod'], 0.)


class ComodulationMasking(Sound):
    """
    Make a stimulus for comodulation masking release.
    Note the parameter names are shortened so that the SGC generated filename
    in cnmodel fits within system limits.
    
    Parameters
    ----------
    rate : float
        sample rate, in Hz
    
    dur : float
        entire waveform duration in seconds
    
    pipst : float
        time to start the test tone pips (seconds)
        array, such as [0.25, 0.35, 0.45]
    
    pipdu : float
        duration of the test (target) tone pips

    maskst : float
        time to start the masker tones pips (seconds)
        array, such as [0.1]
    
    maskdu : float
        duration of the masker tone pips

    rf : float
        rise/fall of the pips 
        
    f0 : float (kHz)
        Center frequency for the target tone, in kHz
    
    db : float
        on-target masker and flankinb band intensity
        In dB SPL (re 0.00002 dynes/cm2)

    s2n : float
        signal re masker, in dbspl

    fmod : float
        amplitude modulation frequency, in Hz
    
    dmod : float
        amplitude modulation depth, in %
 
    fltype : string
        Flanking type:
        One of:
            'MultiTone' : multiple tones, with phase, spacing and # of bands specified as below
            'NBNoise' : the flanking stimulus is made up of narrow band noises (not implemented)
            'None' : no flanking sounds (just on-target stimuli)
    
    flspc : float
        Spacing of flanking bands in octaves from the center frequency, f0
    
    flgap : int
        gap around the cf in flspc (e.g., 1 would skip the first adjacent band)
    
    flph : string
        One of:
            'Comodulated': all of the flanking tones are comodulated in phase 
                with the target at the same frequency and depth
            'Codgh' : 'Grose and Hall codeviant':
                The flanking bands have the same amplitude and frequency
                as the target, but the phase of each band is different. Phases are
                calculated so that all the bands wrap around 2*pi
            'Codvw' : 'Verhey and Winter codeviant':
                The flanking bands have the same amplitude and frequency
                as the target, but are out of phase with the on-frequency masker
            'Random': The phases are selected at random. This is probably best only
                used when there are a large number of flanking bands.
    
    flspl : float
        Flanking signal, in dbspl. 
    
    flN : int
        Number of flanking bands on either side of f0, spaced accordingly.
    
    Returns
    -------
    array :
        waveform
    
    """
    def __init__(self, **kwds):
        # print (kwds)
        for k in ['rate', 'duration', 'pipdu', 'pipst', 'rf', 'maskst', 'maskdu', 
                 # general:
                 'f0', 'dbspl', 's2n', 'fmod', 'dmod',
                 # flankers:
                 'flgap', 'fltype', 'flspc', 'flph', 'flN']:
            if k not in kwds:
                raise TypeError("Missing required argument '%s'" % k)
        if 'flspl' not in kwds:
            kwds['flspl'] = kwds['dbspl']
        # if 'mask_spl' not in kwds:
        #      kwds['mask_spl'] = kwds['dbspl']
        # if kwds['mask_spl'] is None:
        #     kwds['mask_spl'] = 0.
        if kwds['flspl'] is None:
            raise ValueError()

        Sound.__init__(self, **kwds)
    
    def generate(self):
        
        o = self.opts
        # start with center tone
        onfreqmasker = piptone(self.time, o['rf'], o['rate'], o['f0'],
                       o['flspl'], o['maskdu'], o['maskst'])
        tardelay = 0. # 1.5/o['fmod']  # delay by one and one half cycles (no target in first dip)
        target = piptone(self.time, o['rf'], o['rate'], o['f0'],
                       o['dbspl']+o['s2n'], o['pipdu']-tardelay, [p + tardelay for p in o['pipst']])
        # target = shape_signal(onfreqmasker, self.time, o['rf'], o['rate'], o['f0'],
        #                o['dbspl']+o['s2n'], o['pipdu']-tardelay, [p + tardelay for p in o['pipst']])
        if (o['dbspl']+o['s2n']) <= 0.:
            target = np.zeros_like(target)
        onfreqmasker = sinusoidal_modulation(self.time, onfreqmasker, o['maskst'],
            o['fmod'], o['dmod'], 0.)
        # print("o['dbspl']+o['s2n']: ", o['dbspl']+o['s2n'])
        # print('np.max(target): ', np.max(target))

        # target = sinusoidal_modulation(self.time, target, [p + tardelay for p in o['pip_start']],
        #                o['fmod'], o['dmod'], 0.)
        self.onmask = onfreqmasker
        self.target = target
        if o['fltype'] not in ['None', 'Ref', 'NBN', 'Tone']:
            raise ValueError('Unknown flanking_type: %s' % o['fltype'])
        if o['fltype'] in ['NBN']:
            raise ValueError('Flanking type "NBNoise" is not yet implemented')
        elif o['fltype'] in ['None']:
            return (onfreqmasker+target)/2.0  # scaling...
        elif o['fltype'] in ['Tone']:
            nband = int(o['flN'])
            gap = int(o['flgap'])
            octspace = o['flspc']
            f0 = o['f0']
            flankfs = [f0*(2**(octspace*(k+1+gap))) for k in range(nband)]
            flankfs.extend([f0/((2**(octspace*(k+1+gap)))) for k in range(nband)])
            flankfs = sorted(flankfs)
            flanktone = [[]]*len(flankfs)
            for i, fs in enumerate(flankfs):
                flanktone[i] = piptone(self.time, o['rf'], o['rate'], flankfs[i],
                               o['flspl'], o['maskdu'], o['maskst'])
        elif o['fltype'] in ['None', 'Ref']:
            return (onfreqmasker+target)/2.0  # scaling...
        if o['fltype'] == 'NBN':
            raise ValueError('Flanking type nbnoise not yet implemented')
        
        ph = 0.
        if o['flph'] == 'Comod':
                ph = np.zeros(len(flankfs))
        elif o['flph'] == 'Codvw':  # verhey and winter: just 180 out of phase
                ph = np.pi*np.ones(len(flankfs))
        elif o['flph'] in ['Codgh', 'Codev']:  # with precession: Grose and Hall 89
                ph = 2.0*np.pi*np.arange(-o['flN'], o['flN']+1, 1)/o['flN']
        elif o['flph'] == 'Random':
                ph = 2.0*np.pi*np.arange(-o['flN'], o['flN']+1, 1)/o['flN']
                raise ValueError('Random flanking phases not implemented')
        else:
            raise ValueError('Masker Phase pattern of type %s is not implemented' % o['flph'])
        #print(('flanking phases: ', ph))
        #print (len(flanktone))
        #print(('flanking freqs: ', flankfs))
        for i, fs in enumerate(flankfs):
            flanktone[i] = sinusoidal_modulation(self.time, flanktone[i],
                    o['maskst'], o['fmod'], o['dmod'], ph[i])
            if i == 0:
                maskers = flanktone[i]
            else:
                maskers = maskers + flanktone[i]
        signal = (onfreqmasker+maskers+target)/(o['flN']+2)
        return signal


class DynamicRipple(Sound):
    def __init__(self, **kwds):
        for k in ['rate', 'duration']:
            if k not in kwds:
                raise TypeError("Missing required argument '%s'" % k)
        # if kwds['pip_duration'] < kwds['ramp_duration'] * 2:
        #     raise ValueError("pip_duration must be greater than (2 * ramp_duration).")
        import DMR
        self.dmr = DMR.DMR()
        Sound.__init__(self, **kwds)
    
    def generate(self):
        """
        Call to compute a dynamic ripple stimulus
        
        Returns
        -------
        array :
           
           generated waveform
        """
        o = self.opts
        self.dmr.set_params(Fs=o['rate'], duration=o['duration']+1./o['rate'])
        self.dmr.make_waveform()
        self._time = self.dmr.vTime # get time from the generator, not linspace
        return(self.dmr.vStim)


class SpeechShapedNoise(Sound):
    """
    Adapted from http://www.srmathias.com/speech-shaped-noise/
    """
    def __init__(self, **kwds):
        for k in ['rate', 'duration', 'waveform', 'samplingrate']:
            if k not in kwds:
                raise TypeError("Missing required argument '%s'" % k)
        # if kwds['pip_duration'] < kwds['ramp_duration'] * 2:
        #     raise ValueError("pip_duration must be greater than (2 * ramp_duration).")
        Sound.__init__(self, **kwds)
    
    def generate(self):
        o = self.opts
        print('opts: ', o)
        ssn, t = make_ssn(o['rate'], o['duration'], o['waveform'].sound, o['samplingrate'])
        self._time = t  # override time array because we read a wave file
        # if self.opts['duration'] == 0:
        #     self.opts['duration'] = np.max(t) - 1./o['rate']
        return ssn


class RandomSpectrumShape(Sound):
    """
    Random Spectral Shape stimuli
    log-spaced tones
    Amplitudes adjusted in groups of 4 or 8 (amp_group_size)
    Amplitude SD (amp_sd)
    Frequency range (octaves above and below f0) (octaves)
    spacing (fraction of octave: e.g, 1/8 or 1/64 as 8 or 64) (spacing)
    
    Generates one sample
    
    Young and Calhoun, 2005
    Yu and Young, 2000
    """
    def __init__(self, **kwds):
        for k in ['rate', 'duration', 'f0', 'dbspl', 'pip_duration', 'pip_start',
                  'ramp_duration', 'amp_group_size', 'amp_sd', 'spacing', 'octaves']:
            if k not in kwds:
                raise TypeError("Missing required argument '%s'" % k)
        if kwds['pip_duration'] < kwds['ramp_duration'] * 2:
            raise ValueError("pip_duration must be greater than (2 * ramp_duration).")
        if kwds['f0'] > kwds['rate'] * 0.5:
            raise ValueError("f0 must be less than (0.5 * rate).")
        
        Sound.__init__(self, **kwds)
    
    def generate(self):
        o = self.opts
        octaves = o['octaves']
        lowf = o['f0']/octaves
        highf = o['f0']*octaves
        freqlist = np.logspace(np.log2(lowf), np.log2(highf), num=o['spacing']*octaves*2, endpoint=True, base=2)
        amplist = np.zeros_like(freqlist)
        db = o['dbspl']
        # assign amplitudes
        if db == None:
            db = 80.
        groupsize = o['amp_group_size']
        for i in range(0, len(freqlist), groupsize):
            if o['amp_sd'] > 0.:
                a = np.random.normal(scale=o['amp_sd'])
            else:
                a = 0.
            amplist[i:i+groupsize] = 20.0*np.log10(a + db)
        for i in range(len(freqlist)):
            wave = piptone(self.time, o['ramp_duration'], o['rate'], freqlist[i],
                    amplist[i], o['pip_duration'], o['pip_start'], pip_phase=np.pi*2.0*np.random.rand())
            if i == 0:
                result = wave
            else:
                result = result + wave
        return result/(np.sqrt(np.mean(result**2.0))) # scale by rms level
        

class ReadWavefile(Sound):
    """ Read a .wav file from disk, possibly converting the sample rate and the scale
    for use in driving the auditory nerve fiber model.
    
    Parameters
    ----------
    wavefile : str
        name of the .wav file to read
    rate : float
        Sample rate in Hz (waveform will be resampled to this rate)
    channel: int (default: 0)
        If wavefile has 2 channels, select 0 or 1 for the channel to read
    dbspl : float or None
        If specified, the wave file is scaled such that its overall dBSPL
        (measured from RMS of the entire waveform) is equal to this value.
        Either ``dbspl`` or ``scale`` must be specified.
    scale : float or None
        If specified, the wave data is multiplied by this value to yield values in dBSPL. 
        Either ``dbspl`` or ``scale`` must be specified.
    delay: float (default: 0.)
        Silent delay time to start sound, in s. Allows anmodel and cells to run to steady-state
    maxdur : float or None (default: None)
        If specified, maxdur defines the total duration of generated waveform to return (in seconds).
        If None, the generated waveform duration will be the sum of any delay value and
        the duration of the waveform from the wavefile.
    
    Returns
    -------
    array :
        waveform
    
    """
    def __init__(self, wavefile, rate, channel=0, dbspl=None, scale=None, delay=0., maxdur=None):
        if dbspl is not None and scale is not None:
            raise ValueError('Only one of "dbspl" or "scale" can be set')
        duration = 0.  # forced because of the way num_samples has to be calculated first
        if delay < 0.:
            raise ValueError('ReadWavefile: delay must be > 0., got: %f' % delay)
        if maxdur is not None and maxdur < 0:
            raise ValueError('ReadWavefile: maxdur must be None or > 0., got: %f' % maxdur)
        Sound.__init__(self, duration, rate, wavefile=wavefile, channel=channel,
                        dbspl=dbspl, scale=scale,
                        maxdur=maxdur, delay=delay)

    def generate(self):
        """
        Read the wave file from disk, clip duration, resample if necessary, and scale
        
        Returns
        -------
        array :   generated waveform
        """
        [fs_wav, stimulus] = scipy.io.wavfile.read(self.opts['wavefile']) # raw is a numpy array of integer, representing the samples
        if len(stimulus.shape) > 1 and stimulus.shape[1] > 0:
            stimulus = stimulus[:,self.opts['channel']]  # just use selected channel
        fs_wav = float(fs_wav)
        maxdur = self.opts['maxdur']
        delay = self.opts['delay']
        delay_array = np.zeros(int(delay*fs_wav))  # build delay array (may have 0 length)
        if maxdur is None:
            maxdur = delay + len(stimulus)/fs_wav  # true total length
        maxpts = int(maxdur * fs_wav)
        stimulus = np.hstack((delay_array, stimulus))[:maxpts]

        if self.opts['rate'] != fs_wav:
            stimulus = resampy.resample(stimulus, fs_wav, self.opts['rate'])
        self.opts['duration'] = (stimulus.shape[0]-1)/self.opts['rate'] # compute the duration, match for linspace calculation used in time.
        self._time = None
        self.time   # requesting time should cause recalulation of the time
        if self.opts['dbspl'] is not None:
            rms = np.sqrt(np.mean(stimulus**2.0))  # find rms of the waveform
            stimulus = dbspl_to_pa(self.opts['dbspl'] ) * stimulus / rms # scale into Pascals
        if self.opts['scale'] is not None:
            stimulus = stimulus * self.opts['scale']
        return stimulus


def sinusoidal_modulation(t, basestim, tstart, fmod, dmod, phaseshift):
    """
    Impose a sinusoidal amplitude-modulation on the input waveform.
    For dmod=100%, the envelope max is 2, the min is 0; for dmod = 0, the max and min are 1
    maintains equal energy for all modulation depths.
    Equation from Rhode and Greenberg, J. Neurophys, 1994 (adding missing parenthesis) and
    Sayles et al. J. Physiol. 2013
    The envelope can be phase shifted (useful for co-deviant stimuli).
    
    Parameters
    ----------
    t : array
        array of waveform time values (seconds)
    basestim : array
        array of waveform values that will be subject to sinulsoidal envelope modulation
    tstart : float
        time at which the base sound starts (modulation starts then, with 0 phase crossing)
        (seconds)
    fmod : float
        modulation frequency (Hz)
    dmod : float
        modulation depth (percent)
    phaseshift : float
        modulation phase shift (starting phase, radians)
    
    """

    env = (1.0 + (dmod/100.0) * np.sin((2.0*np.pi*fmod*(t-tstart)) + phaseshift - np.pi/2)) # envelope...
    return basestim*env


def modnoise(t, rt, Fs, F0, dur, start, dBSPL, FMod, DMod, phaseshift, seed):
    """
    Generate an amplitude-modulated noise with linear ramps.
    
    Parameters
    ----------
    t : array
        array of waveform time values
    rt : float
        ramp duration
    Fs : float
        sample rate
    F0 : float
        tone frequency
    dur : float
        duration of noise
    start : float
        start time for noise
    dBSPL : float
        sound pressure of stimulus
    FMod : float
        modulation frequency
    DMod : float
        modulation depth percent
    phaseshift : float
        modulation phase
    seed : int
        seed for random number generator
    
    Returns
    -------
    array :
        waveform
    
    """
    irpts = int(rt * Fs)
    mxpts = len(t)+1
    pin = pipnoise(t, rt, Fs, dBSPL, dur, start, seed)
    env = (1 + (DMod/100.0) * np.sin((2*np.pi*FMod*t) - np.pi/2 + phaseshift)) # envelope...

    pin = linearramp(pin, mxpts, irpts)
    env = linearramp(env, mxpts, irpts)
    return pin*env


def linearramp(pin, mxpts, irpts):
    """
    Apply linear ramps to *pin*.
    
    Parameters
    ----------
    pin : array
        input waveform to apply ramp to

    mxpts : int
        point in array to start ramp down
    
    irpts : int
        duration of the ramp

    Returns
    -------
    array :
        waveform
    
    
    Original (adapted from Manis; makeANF_CF_RI.m)::
    
        function [out] = ramp(pin, mxpts, irpts)
            out = pin;
            out(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts;
            out((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;
            return;
        end
    """
    out = pin.copy()
    r = np.linspace(0, 1, irpts)
    irpts = int(irpts)
    # print 'irpts: ', irpts
    # print len(out)
    out[:irpts] = out[:irpts]*r
    # print  out[mxpts-irpts:mxpts].shape
    # print r[::-1].shape
    out[mxpts-irpts-1:mxpts] = out[mxpts-irpts-1:mxpts] * r[::-1]
    return out


def pipnoise(t, rt, Fs, dBSPL, pip_dur, pip_start, seed):
    """
    Create a waveform with multiple sine-ramped noise pips. Output is in 
    Pascals.
    
    Parameters
    ----------
    t : array
        array of time values
    rt : float
        ramp duration 
    Fs : float
        sample rate
    dBSPL : float
        maximum sound pressure level of pip
    pip_dur : float
        duration of pip including ramps
    pip_start : float
        list of starting times for multiple pips
    seed : int
        random seed

    Returns
    -------
    array :
        waveform

    """
    rng = np.random.RandomState(seed)
    pin = np.zeros(t.size)
    for start in pip_start:
        # make pip template
        pip_pts = int(pip_dur * Fs) + 1
        pip = dbspl_to_pa(dBSPL) * rng.randn(pip_pts)  # unramped stimulus

        # add ramp
        ramp_pts = int(rt * Fs) + 1
        ramp = np.sin(np.linspace(0, np.pi/2., ramp_pts))**2
        pip[:ramp_pts] *= ramp
        pip[-ramp_pts:] *= ramp[::-1]
        
        ts = int(np.floor(start * Fs))
        pin[ts:ts+pip.size] += pip

    return pin
        
   
def piptone(t, rt, Fs, F0, dBSPL, pip_dur, pip_start):
    """
    Create a waveform with multiple sine-ramped tone pips. Output is in 
    Pascals.
    
    Parameters
    ----------
    t : array
        array of time values
    rt : float
        ramp duration 
    Fs : float
        sample rate
    F0 : float
        pip frequency
    dBSPL : float
        maximum sound pressure level of pip
    pip_dur : float
        duration of pip including ramps
    pip_start : float
        list of starting times for multiple pips

    Returns
    -------
    array :
        waveform

    """
    # make pip template
    pip_pts = int(pip_dur * Fs) + 1
    pip_t = np.linspace(0, pip_dur, pip_pts)
    pip = np.sqrt(2) * dbspl_to_pa(dBSPL) * np.sin(2*np.pi*F0*pip_t)  # unramped stimulus

    # add ramp
    ramp_pts = int(rt * Fs) + 1
    ramp = np.sin(np.linspace(0, np.pi/2., ramp_pts))**2
    pip[:ramp_pts] *= ramp
    pip[-ramp_pts:] *= ramp[::-1]
    
    # apply template to waveform
    pin = np.zeros(t.size)
    ps = pip_start
    if ~isinstance(ps, list):
        ps = [ps]
    for start in pip_start:
        ts = int(np.floor(start * Fs))
        pin[ts:ts+pip.size] += pip

    return pin

def shape_signal(signal, t, rt, Fs, F0, dBSPL, pip_dur, pip_start):
    """
    Create a waveform with multiple sine-ramped tone pips, based on the 
    signal (so output is *always* in phase with the reference signal waveform)
    Output is in Pascals.
    
    Parameters
    ----------
    t : array
        array of time values
    rt : float
        ramp duration (risetime)
    Fs : float
        sample rate
    F0 : float
        pip frequency
    dBSPL : float
        maximum sound pressure level of pip
    pip_dur : float
        duration of pip including ramps
    pip_start : float
        list of starting times for multiple pips

    Returns
    -------
    array :
        waveform

    """
    pip_t = t
    pip_pts = pip_t.shape[0]
    pip = np.sqrt(2) * dbspl_to_pa(dBSPL) * signal # referencestimulus

    # make envelope with cos2 rise-fall
    ramp_pts = int(rt * Fs) + 1
    env_pts = int(pip_dur * Fs)
    envelope = np.ones(env_pts)
    ramp = np.sin(np.linspace(0, np.pi/2., ramp_pts))**2
    envelope[:ramp_pts] *= ramp
    envelope[-ramp_pts:] *= ramp[::-1]
    
    # apply envelope template to waveform
    pin = np.zeros(t.size)
    ps = pip_start
    if ~isinstance(ps, list):
        ps = [ps]
    for start in pip_start:
        ts = int(np.floor(start * Fs))
        pin[ts:ts+envelope.size] += pip[ts:ts+envelope.size]*envelope

    return pin

def click(t, Fs, dBSPL, click_duration, click_start):
    """
    Create a waveform with one retangular click. Output is in 
    Pascals.
    
    Parameters
    ----------
    t : array
        array of time values
    Fs : float
        sample frequency (Hz)
    dspl : float
        maximum sound pressure level of c
    click_start : float (seconds)
        delay to the click train 
    click_duration : float (seconds)
        duration of the click
 
 
    Returns
    -------
    array :
        waveform

    """
    swave = np.zeros(t.size)
    amp = dbspl_to_pa(dBSPL)
    t0 = int(np.floor(click_start * Fs))
    td = int(np.floor(click_duration * Fs))
    print('to, td: ', t0, td)
    swave[t0:t0+td] = amp
    return swave
    # nclicks = len(click_starts)
    # for n in range(nclicks):
    #     t0s = click_starts[n]  # time for nth click
    #     t0 = int(np.floor(t0s * Fs))  # index
    #     if t0+td > t.size:
    #         raise ValueError('Clicks: train duration exceeds waveform duration')
    #     swave[t0:t0+td] = amp
    # return swave


def fmsweep(t, start, duration, freqs, ramp, dBSPL):
    """
    Create a waveform for an FM sweep over time. Output is in 
    Pascals.

    Parameters
    ----------
    t : array
        time array for waveform 
    start : float (seconds)
        start time for sweep
    duration : float (seconds)
        duration of sweep
    freqs : array (Hz)
        Two-element array specifying the start and end frequency of the sweep
    ramp : str
        The shape of time course of the sweep (linear, logarithmic)
    dBSPL : float
        maximum sound pressure level of sweep
    
    Returns
    -------
    array :
        waveform
    
    
    """
    
    sw = scipy.signal.chirp(t, freqs[0], duration, freqs[1],
        method=ramp, phi=0, vertex_zero=True)
    sw = np.sqrt(2) * dbspl_to_pa(dBSPL) * sw
    return sw


def make_ssn(rate, duration, sig, samplingrate):
        """
        Speech-shaped noise
        Adapted from http://www.srmathias.com/speech-shaped-noise/
        Created on Thu Jun 26 12:42:08 2014
        @author: smathias
        """
        # note rate is currently ignored...
        sig = np.array(sig).astype('float64')
        if rate != samplingrate:  # interpolate to the current system sampling rate from the original rate
            sig = np.interp(np.arange(0, len(sig)/rate, 1./rate),
                np.arange(0, len(sig)/samplingrate), 1./samplingrate)
        sig = 2*sig/np.max(sig)
        z, t = noise_from_signal(sig, rate, keep_env=True)
        return z, t


def noise_from_signal(x, fs=40000, keep_env=True):
    """
    Create a noise with same spectrum as the input signal.
    
    Parameters
    ----------
    x : array_like
        Input signal.
    fs : int
         Sampling frequency of the input signal. (Default value = 40000)
    keep_env : bool
         Apply the envelope of the original signal to the noise. (Default
         value = False)
    
    Returns
    -------
    ndarray
        Noise signal.
    """
    x = np.asarray(x)
    n_x = x.shape[-1]
    n_fft = next_pow_2(n_x)
    X = np.fft.rfft(x, next_pow_2(n_fft))
    # Randomize phase.
    noise_mag = np.abs(X) * np.exp(
        2. * np.pi * 1j * np.random.random(X.shape[-1]))
    noise = np.real(np.fft.irfft(noise_mag, n_fft))
    out = noise[:n_x]
    if keep_env:
        env = np.abs(scipy.signal.hilbert(x))
        [bb, aa] = scipy.signal.butter(6., 50. / (fs / 2.))  # 50 Hz LP filter
        env = scipy.signal.filtfilt(bb, aa, env)
        out *= env
    t = np.arange(0, (len(out))/fs, 1./fs)
    return out, t

