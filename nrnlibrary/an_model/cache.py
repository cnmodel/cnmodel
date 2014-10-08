"""
Utilities for generating and caching spike trains from AN model.
"""

import os, pickle
import numpy as np
from .wrapper import get_matlab, model_ihc, model_synapse, seed_rng

_cache_version = 1
_cache_path = os.path.join(os.path.dirname(__file__), 'cache')
_index_file = os.path.join(_cache_path, 'index.pk')
_index = None


def get_spiketrain(cf, sr, stim, seed, **kwds):
    """ Return an array of spike times in response to the given stimulus.
    
    Arrays are automatically cached and may be returned from disk if 
    available. See generate_spiketrain() for a description of arguments.
    """
    index = cache_index()
    keydata = dict(cf=cf, sr=sr, stim=stim.key(), seed=seed, **kwds)
    key = make_key(**keydata)
    data = None
    
    # Load data from cache if possible
    if key in index:
        data_file = index[key]['file']
        if os.path.exists(data_file):
            data = np.load(open(data_file, 'rb'))['data']
   
    # Generate new data if needed
    if data is None:
        data = generate_spiketrain(cf, sr, stim, seed, **kwds)
        
        # store new data to cache
        subdir = os.path.join(_cache_path, make_key(**stim.key()))
        filename = make_key(cf=cf, sr=sr, seed=seed, **kwds)
        filename = os.path.join(subdir, filename) + '.npz'
        
        if not os.path.exists(subdir):
            os.mkdir(subdir)
        np.savez_compressed(filename, data=data)
        
        # update the index
        keydata['file'] = filename
        index[key] = keydata
        save_index()
    
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
    kwds = list(kwds.items())
    kwds.sort()
    return '_'.join(['%s=%s' % kv for kv in kwds])


def cache_index():
    """ Return an index that gives the file name for each stored cache entry.
    """
    global _index, _cache_path, _cache_version, _index_file
    if _index is None:
        if not os.path.isdir(_cache_path):
            os.mkdir(_cache_path)
            
        if os.path.isfile(_index_file):
            _index = pickle.load(open(_index_file, 'rb'))
            if _index['_cache_version'] != _cache_version:
                i = 0
                while True:
                    old_cache = _cache_path + '.old-%d' % i
                    if not os.path.exists(old_cache):
                        break
                os.rename(_cache_path, old_cache)
                print ("Cache version is too old; starting new cache. "
                       "(old cache is stored at %s)" % old_cache)
                return cache_index()
        else:
            _index = {'_cache_version': _cache_version}
    return _index


def save_index():
    """ Write the index to file
    """
    global _index_file, _index
    pickle.dump(_index, open(_index_file, 'wb'))


def generate_spiketrain(cf, sr, stim, seed, **kwds):
    """ Generate a new spike train from the auditory nerve model.
    
    Parameters
    ----------
    cf : float
        Center frequency of the fiber to simulate
    sr : int
        Spontaneous rate group of the fiber: 
        1=low, 2=mid, 3=high.
    stim : Sound instance
        Stimulus sound to be presented on each repetition
    seed : int >= 0
        Random seed
        
    All other keyword arguments are given to model_ihc() and model_synapse()
    based on their names. These include 'species', 'nrep', 'reptime', 'cohc', 
    'cihc', and 'implnt'. 
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
    
    if len(kwds) > 0:
        raise TypeError("Invalid keyword arguments: %s" % list(kwds.keys()))
    
    seed_rng(seed)
    vihc = model_ihc(_transfer=False, **ihc_kwds) 
    m, v, psth = model_synapse(vihc, _transfer=False, **syn_kwds)
    psth = psth.get()
    
    return np.argwhere(psth).ravel()
