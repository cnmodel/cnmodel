import os, tempfile
from multiprocessing import Pool
import cnmodel.an_model.cache as cache
import cnmodel.util.sound as sound
from cnmodel import an_model

def new_cache():
    # Do all cache testing on a temporary cache directory
    cdir = tempfile.mkdtemp()
    cache._cache_path = cdir


def run_sim(*args):
    cf = 1000
    sr = 2
    seed = 12345678
    stim = sound.TonePip(rate=100e3, duration=0.01, f0=4000, dbspl=80,
                         ramp_duration=0.002, pip_duration=0.004, 
                         pip_start=[0.001])
    spikes = an_model.get_spiketrain(cf=cf, sr=sr, seed=seed, stim=stim)
    cf = cache.get_cache_filename(cf, sr, seed, stim)
    mtime = os.stat(cf).st_mtime
    return cf, mtime, spikes


def test_cache():
    new_cache()
    cfile1, mtime1, spikes1 = run_sim()
    cfile2, mtime2, spikes2 = run_sim()
    assert cfile1 == cfile2
    assert mtime1 == mtime2
    assert all(spikes1 == spikes2)


def test_parallel():
    # Make sure file locking works correctly.
    return
    new_cache()  # note that subprocesses will all inherit this new cache 
    p = Pool(10)
    results = p.map(run_sim, range(10))
    for i in range(1, 10):
        assert results[0][0] == results[i][0]
        assert results[0][1] == results[i][1]
        assert all(results[0][2] == results[i][2])


if __name__ == '__main__':
    #test_cache()
    test_locking()
