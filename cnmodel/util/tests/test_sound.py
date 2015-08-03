import numpy as np
from cnmodel.util import sound

def test_conversions():
    pa = np.array([3990.5, 20, 0.3639, 2e-5])
    db = np.array([ 166, 120,  85.2,   0])
        
    assert np.allclose(sound.pa_to_dbspl(pa), db, atol=0.1, rtol=0.002)
    assert np.allclose(sound.dbspl_to_pa(db), pa, atol=0.1, rtol=0.002)
    
    
def test_tonepip():
    rate = 100000
    dur = 0.1
    ps = 0.01
    rd = 0.02
    pd = 0.08
    db = 60
    s1 = sound.TonePip(rate=rate, duration=dur, f0=5321, dbspl=db, 
                       pip_duration=pd, pip_start=[ps], ramp_duration=rd)

    # test array sizes
    assert s1.sound.size == s1.time.size == int(dur * rate) + 1

    # test for consistency
    assert np.allclose([s1.sound.min(), s1.sound.mean(), s1.sound.max()], 
                       [-0.028284253158247834, -1.0954891976168953e-10, 0.028284270354167296])
    
    # test that we got the requested amplitude
    assert np.allclose(s1.measure_dbspl(ps+rd, ps+pd-rd), db, atol=0.1, rtol=0.01)
    
    # test for quiet before and after pip
    assert np.all(s1.sound[:(ps*rate)-1] == 0)
    assert np.all(s1.sound[((ps+pd)*rate)+1:] == 0)
    
    # test the sound can be recreated from its key
    key = s1.key()
    s2 = sound.create(**key)
    assert np.all(s1.time == s2.time)
    assert np.all(s1.sound == s2.sound)


def test_noisepip():
    rate = 100000
    dur = 0.1
    ps = 0.01
    rd = 0.02
    pd = 0.08
    db = 60
    s1 = sound.NoisePip(rate=rate, duration=dur, seed=184724, dbspl=db, 
                        pip_duration=pd, pip_start=[ps], ramp_duration=rd)

    # test array sizes
    assert s1.sound.size == s1.time.size == int(dur * rate) + 1

    # test for consistency
    assert np.allclose([s1.sound.min(), s1.sound.mean(), s1.sound.max()], 
                       [-0.082260796003197786, -0.00018484322982972046, 0.069160217220832404])
    
    # test that we got the requested amplitude
    assert np.allclose(s1.measure_dbspl(ps+rd, ps+pd-rd), db, atol=0.1, rtol=0.01)
    
    # test for quiet before and after pip
    assert np.all(s1.sound[:(ps*rate)-1] == 0)
    assert np.all(s1.sound[((ps+pd)*rate)+1:] == 0)
    
    # test the sound can be recreated from its key
    key = s1.key()
    s2 = sound.create(**key)
    # also test new seed works, and does not affect other sounds
    key['seed'] += 1
    s3 = sound.create(**key)
    s3.sound  # generate here to advance rng before generating s2
    
    assert np.all(s1.time == s2.time)
    assert np.all(s1.sound == s2.sound)
    start = (ps * rate) + 1
    end = ((ps+pd) * rate) - 1
    assert not np.any(s1.sound[start:end] == s3.sound[start:end])


    