"""
Tools for generating auditory stimuli. 
"""
from __future__ import division
import numpy as np


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
        self.opts = {'rate': rate, 'duration': duration}
        self.opts.update(kwds)
        self._time = None
        self._sound = None

    @property
    def sound(self):
        """ The generated sound array expressed in Pascals.
        """
        if self._sound is None:
            self._sound = self.generate()
        return self._sound
        
    @property
    def time(self):
        """ The array of time values expressed in seconds.
        """
        if self._time is None:
            self._time = np.linspace(0, self.opts['duration'], self.num_samples)
        return self._time

    @property 
    def num_samples(self):
        """ The number of samples in the sound array.
        """
        return 1 + int(self.opts['duration'] * self.opts['rate'])
    
    @property
    def dt(self):
        """ The sample period (time step between samples).
        """
        return 1.0 / self.opts['rate']

    def key(self):
        """ Return dict of parameters needed to completely describe this sound.
        The sound can be recreated using ``Sound.create(**key)``.
        """
        k = self.opts.copy()
        k['type'] = self.__class__.__name__
        return k

    def measure_dbspl(self, tstart, tend):
        """ Return the measured amplitude (dBSPL) of the sound from tstart to tend
        (both specified in seconds). 
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
        if name in self.opts:
            return self.opts[name]
        else:
            return object.__getattr__(self, name)


class TonePip(Sound):
    """ One or more tone pips with cosine-ramped edges.
    
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
        for k in ['rate', 'duration', 'f0', 'dbspl', 'pip_duration', 'pip_start', 'ramp_duration']:
            if k not in kwds:
                raise TypeError("Missing required argument '%s'" % k)
        if kwds['pip_duration'] < kwds['ramp_duration'] * 2:
            raise ValueError("pip_duration must be greater than (2 * ramp_duration).")
        if kwds['f0'] > kwds['rate'] * 0.5:
            raise ValueError("f0 must be less than (0.5 * rate).")
        
        Sound.__init__(self, **kwds)
        
    def generate(self):
        o = self.opts
        return piptone(self.time, o['ramp_duration'], o['rate'], o['f0'], 
                       o['dbspl'], o['pip_duration'], o['pip_start'])
    

class NoisePip(Sound):
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
    """
    def __init__(self, **kwds):
        for k in ['rate', 'duration', 'seed', 'pip_duration', 'pip_start', 'ramp_duration']:
            if k not in kwds:
                raise TypeError("Missing required argument '%s'" % k)
        if kwds['pip_duration'] < kwds['ramp_duration'] * 2:
            raise ValueError("pip_duration must be greater than (2 * ramp_duration).")
        if kwds['seed'] < 0:
            raise ValueError("Random seed must be integer > 0")
        
        Sound.__init__(self, **kwds)
        
    def generate(self):
        o = self.opts
        return pipnoise(self.time, o['ramp_duration'], o['rate'],
                        o['dbspl'], o['pip_duration'], o['pip_start'], o['seed'])
    

def pa_to_dbspl(pa, ref=20e-6):
    """ Convert Pascals (rms) to dBSPL. By default, the reference pressure is
    20 uPa.
    """
    return 20 * np.log10(pa / ref)


def dbspl_to_pa(dbspl, ref=20e-6):
    """ Convert dBSPL to Pascals (rms). By default, the reference pressure is
    20 uPa.
    """
    return ref * 10**(dbspl / 20)


def modtone(t, rt, Fs, F0, dBSPL, FMod, DMod, phaseshift):
    """
    Generate an amplitude-modulated tone with linear ramps.
    
    t: waveform time values
    rt: ramp duration
    Fs: sample rate
    F0: tone frequency
    FMod : modulation frequency
    DMod : modulation depth percent
    phaseshift : modulation phase
    
    Original (adapted from Manis; makeANF_CF_RI.m): 
    
    function [pin, env] = modtone(t, rt, Fs, F0, dBSPL, FMod, DMod, phaseshift)
        % fprintf(1, 'Phase: %f\n', phaseshift)
        irpts = rt*Fs;
        mxpts = length(t);
        env = (1 + (DMod/100.0)*sin((2*pi*FMod*t)-pi/2+phaseshift)); % envelope...
        pin = sqrt(2)*20e-6*10^(dBSPL/20)*(sin((2*pi*F0*t)-pi/2).*env); % unramped stimulus

        pin = ramp(pin, mxpts, irpts);
        env = ramp(env, mxpts, irpts);
        %pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts;
        %pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;
        return
    end
    """
    irpts = rt * Fs
    mxpts = len(t)
    
    # TODO: is this envelope correct? For dmod=100, the envelope max is 2.
    # I would have expected something like  (dmod/100) * 0.5 * (sin + 1)
    env = (1 + (DMod/100.0) * np.sin((2*pi*FMod*t) - np.pi/2 + phaseshift)) # envelope...
    
    pin = (np.sqrt(2) * dbspl_to_pa(dBSPL)) * np.sin((2*pi*F0*t) - np.pi/2) * env # unramped stimulus
    pin = ramp(pin, mxpts, irpts)
    env = ramp(env, mxpts, irpts)
    return pin, env


def ramp(pin, mxpts, irpts):
    """
    Apply linear ramps to *pin*.
    
    Original (adapted from Manis; makeANF_CF_RI.m): 
    
    function [out] = ramp(pin, mxpts, irpts)
        out = pin;
        out(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts;
        out((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;
        return;
    end
    """
    out = pin.copy()
    r = np.linspace(0, 1, irpts)
    out[:irpts] *= r
    out[mxpts-irpts:mxpts] *= r[::-1]
    return out

def pipnoise(t, rt, Fs, dBSPL, pip_dur, pip_start, seed):
    """
    Create a waveform with multiple sine-ramped noise pips. Output is in 
    Pascals.
    
    t: array of time values
    rt: ramp duration 
    Fs: sample rate
    dBSPL: maximum sound pressure level of pip
    pip_dur: duration of pip including ramps
    pip_start: list of starting times for multiple pips
    seed: random seed
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
        
        ts = np.floor(start * Fs)
        pin[ts:ts+pip.size] += pip

    return pin
        
   
def piptone(t, rt, Fs, F0, dBSPL, pip_dur, pip_start):
    """
    Create a waveform with multiple sine-ramped tone pips. Output is in 
    Pascals.
    
    t: array of time values
    rt: ramp duration 
    Fs: sample rate
    F0: pip frequency
    dBSPL: maximum sound pressure level of pip
    pip_dur: duration of pip including ramps
    pip_start: list of starting times for multiple pips
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
    for start in pip_start:
        ts = np.floor(start * Fs)
        pin[ts:ts+pip.size] += pip

    return pin
    