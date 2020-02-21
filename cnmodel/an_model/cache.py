"""
Utilities for generating and caching spike trains from AN model.
"""

import logging
import os, sys, pickle
import numpy as np
from .wrapper import get_matlab, model_ihc, model_synapse, seed_rng
from ..util.filelock import FileLock

import cochlea
try:
    import cochlea
    HAVE_COCHLEA = True
except ImportError:
    HAVE_COCHLEA = False

_cache_version = 2
_cache_path = os.path.join(os.path.dirname(__file__), 'cache')
_index_file = os.path.join(_cache_path, 'index.pk')
_index = None


def get_spiketrain(cf, sr, stim, seed, verbose=False, **kwds):
    """ Return an array of spike times in response to the given stimulus.
    
    Arrays are automatically cached and may be returned from disk if 
    available. See generate_spiketrain() for a description of arguments.
    
    If the flag --ignore-an-cache was given on the command line, then spike 
    times will be regenerated and cached, regardless of the current cache 
    state.
    
    If the flag --no-an-cache was given on the command line, then the cache
    will not be read or written. This can improve overall performance if there
    is little chance the cache would be re-used.
    
    """
    filename = get_cache_filename(cf=cf, sr=sr, seed=seed, stim=stim, **kwds)
    subdir = os.path.dirname(filename)
    
    if not os.path.exists(subdir):
        try:
            os.mkdir(subdir)
        except OSError as err:
            # probably another process already created this directory
            # since we last checked
            pass
    
    with FileLock(filename):
        
        if '--ignore-an-cache' in sys.argv or '--no-an-cache' in sys.argv or not os.path.exists(filename):
            create = True
        else:
            create = False
            # try loading cached data
            try:
                data = np.load(open(filename, 'rb'))['data']
                if verbose:
                    logging.info("Loaded AN spike train from cache: %s", filename)
            except Exception:
                create = True
                sys.excepthook(*sys.exc_info())
                logging.error("Error reading AN spike train cache file; will "
                    "re-generate. File: %s", filename)

        if create:
            if verbose:
                logging.info("Generate new AN spike train: %s", filename)
            data = generate_spiketrain(cf, sr, stim, seed, **kwds)
            if '--no-an-cache' not in sys.argv:
                np.savez_compressed(filename, data=data)
            
    return data


def make_key(**kwds):
    """ Make a unique key used for caching spike time arrays.
    """
    # flatten any nested dicts
    for key in list(kwds.keys()):
        if isinstance(kwds[key], dict):
            val = kwds.pop(key)
            for k,v in val.items():
                kwds[key + '.' + k] = v
    # sort and convert to string
    for k in kwds:
        if k == 'click_starts':
            x = kwds[k]
            xd = np.mean(np.diff(x))
            kwds[k] = f"{x[0]:0.2f}-{x[-1]:0.2f}by{xd:0.2f}"
    kwds = list(kwds.items())
    kwds.sort()
    return '_'.join(['%s=%s' % kv for kv in kwds])


def get_cache_filename(cf, sr, seed, stim, **kwds):
    global _cache_path
    # print('**stim.key(): ', make_key(**stim.key()))
    subdir = os.path.join(_cache_path, make_key(**stim.key()))
    filename = make_key(cf=cf, sr=sr, seed=seed, **kwds)
    filename = os.path.join(subdir, filename) + '.npz'
    return filename


def generate_spiketrain(cf, sr, stim, seed, simulator=None, **kwds):
    """ Generate a new spike train from the auditory nerve model. Returns an 
    array of spike times in seconds.
    
    Parameters
    ----------
    cf : float
        Center frequency of the fiber to simulate
    sr : int
        Spontaneous rate group of the fiber: 
        0=low, 1=mid, 2=high.
    stim : Sound instance
        Stimulus sound to be presented on each repetition
    seed : int >= 0
        Random seed
    simulator : 'cochlea' | 'matlab' | None
        Specifies the auditory periphery simulator to use. If None, then a
        simulator will be automatically chosen based on availability.
    **kwds : 
        All other keyword arguments are given to model_ihc() and model_synapse()
        based on their names. These include 'species', 'nrep', 'reptime', 'cohc', 
        'cihc', and 'implnt'. 
        'simulator' is used to set the simulator ('matlab' or 'cochlea')
    """
    
    for k in ['pin', 'CF', 'fiberType', 'noiseType']:
        if k in kwds:
            raise TypeError("Argument '%s' is not allowed here." % k)
    
    ihc_kwds = dict(pin=stim.sound, CF=cf, nrep=1, tdres=stim.dt, 
                    reptime=stim.duration*2, cohc=1, cihc=1, species=1)
    syn_kwds = dict(CF=cf, nrep=1, tdres=stim.dt, fiberType=sr, noiseType=1, implnt=0)
    # copy any given keyword args to the correct model function
    for kwd in kwds:
        if kwd in ihc_kwds:
            ihc_kwds[kwd] = kwds.pop(kwd)
        if kwd in syn_kwds:
            syn_kwds[kwd] = kwds.pop(kwd)

    if simulator is None:
        simulator = detect_simulator()

    if len(kwds) > 0:
        raise TypeError("Invalid keyword arguments: %s" % list(kwds.keys()))
    
    if simulator in ['MATLAB', 'matlab']:
        seed_rng(seed)
        vihc = model_ihc(_transfer=False, **ihc_kwds) 
        m, v, psth = model_synapse(vihc, _transfer=False, **syn_kwds)
        psth = psth.get().ravel()
        times = np.argwhere(psth).ravel()
        return times * stim.dt
    elif (simulator in ['cochlea']) and HAVE_COCHLEA:
        fs = int(0.5+1./stim.dt)  # need to avoid roundoff error
        srgrp = [0,0,0] # H, M, L (but input is 1=L, 2=M, H = 3)
        srgrp[2-sr] = 1
        sp = cochlea.run_zilany2014(
                stim.sound,
                fs=fs,
                anf_num=srgrp,
                cf=cf,
                seed=seed,
                species='cat')
        return np.array(sp.spikes.values[0])
    else:  # it remains possible to have a typo.... 
        raise ValueError("anmodel/cache.py: Simulator must be specified as either MATLAB or cochlea; found <%s> of type %s (cochlea? %r)"
            % (simulator, type(simulator), HAVE_COCHLEA))


def detect_simulator():
    """Return the name of any available auditory periphery model.
    
    Return 'cochlea' if the Rudnicki cochlea model can be imported. 
    If not, return 'matlab' if the Zilany model can be accessed via MATLAB.
    If not, raise an exception.
    """
    try:
        import cochlea
        simulator = 'cochlea'
    except ImportError:
        get_matlab()
        simulator = 'matlab'
    return simulator
