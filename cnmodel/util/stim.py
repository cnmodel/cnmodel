import numpy as np


def make_pulse(stim):
    """
    Generate a pulse train for current / voltage command. Returns a tuple.
    
    
    Parameters
    ----------
    stim : dict
        Holds parameters that determine stimulus shape:
        
        * delay : time before first pulse
        * Sfreq : frequency of pulses
        * dur : duration of one pulse or main pulse
        * predur : duration of prepulse (default should be 0 for no prepulse)
        * amp : pulse amplitude
        * preamp : amplitude of prepulse
        * PT : delay between end of train and test pulse (0 for no test)
        * NP : number of pulses
        * hold : holding level (optional)
        * dt : timestep

    Returns
    -------
    w : stimulus waveform
    maxt : duration of waveform
    tstims : index of each pulse in the train
    """
    defaults = {
        'delay': 10,
        'Sfreq': 50,
        'dur': 100,
        'predur': 0.,
        'post': 50.,
        'amp': None,
        'preamp': 0.,
        'PT': 0,
        'NP': 1,
        'hold': 0.0,
        'dt': None,
        }
    for k in stim:
        if k not in defaults:
            raise Exception("Stim parameter '%s' not accepted." % k)
    defaults.update(stim)
    stim = defaults
    for k,v in stim.items():
        if v is None:
            raise Exception("Must specify stim parameter '%s'." % k)
    
    dt = stim['dt']
    delay = int(np.floor(stim['delay'] / dt))
    ipi = int(np.floor((1000.0 / stim['Sfreq']) / dt))
    pdur = int(np.floor(stim['dur'] / dt))
    posttest = int(np.floor(stim['PT'] / dt))
    ndur = 5
    if stim['predur'] > 0.:
        predur = int(np.floor(stim['predur'] / dt))
    else:
        predur = 0.
    if stim['PT'] == 0:
        ndur = 1

    maxt = dt * (delay + predur + (ipi * (stim['NP'] + 3)) +
        posttest + pdur * ndur)
    hold = stim.get('hold', None)
    
    w = np.zeros(int(np.floor(maxt / dt)))
    if hold is not None:
        w += hold
    
    #   make pulse
    tstims = [0] * int(stim['NP'])
    for j in range(0, int(stim['NP'])):
        prestart = delay
        start = int(prestart + predur + j*ipi)
        if predur > 0.:
            w[prestart:prestart+predur] = stim['preamp']
        w[start:start + pdur] = stim['amp']
        tstims[j] = start
    if stim['PT'] > 0.0:
        for i in range(start + posttest, start + posttest + pdur):
            w[i] = stim['amp']
    w = np.append(w, 0.)
    maxt = maxt + dt
    return(w, maxt, tstims)

