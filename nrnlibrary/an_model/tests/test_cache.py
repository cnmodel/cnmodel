import os, tempfile
from multiprocessing import Pool
import nrnlibrary.an_model.cache as cache
import nrnlibrary.util.sound as sound
from nrnlibrary import an_model

def new_cache():
    # Do all cache testing on a temporary cache directory
    cdir = tempfile.mkdtemp()
    cache._index_file = os.path.join(cdir, 'index.pk')
    cache._cache_path = os.path.join(cdir, 'cache')


def test_cache():
    new_cache()
    stim = sound.TonePip(rate=100e3, duration=0.01, f0=4000, dbspl=80,
                         ramp_duration=0.002, pip_duration=0.004, 
                         pip_start=[0.001])
    spikes1 = an_model.get_spiketrain(cf=1000, sr=2, seed=12345678, stim=stim)
    spikes2 = an_model.get_spiketrain(cf=1000, sr=2, seed=12345678, stim=stim)
    assert all(spikes1 == spikes2)


def modify_cache(x):
    index = cache.cache_index()
    index[x] = x
    cache.save_index()
    index = cache.cache_index(reload=True)


def test_locking():
    # Make sure file locking works correctly.
    new_cache()  # note that subprocesses will all inherit this new cache 
    p = Pool(10)
    p.map(modify_cache, range(10))
    index = cache.cache_index()
    for i in range(10):
        assert index[i] == i


if __name__ == '__main__':
    #test_cache()
    test_locking()
