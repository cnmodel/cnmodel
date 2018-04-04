__author__ = 'pbmanis'
"""
ccstim
Generate current-clamp (or voltage-clamp) stimulus waveforms from a dictionary
used for vectory play modes in current clamp
(prior version was called 'makestim')

Can generate several types of pulses

"""

import numpy as np


def ccstim(stim, dt, pulsetype='square'):
    """
    Create stimulus pulse waveforms of different types.
    
    Parameters
    ----------
    stim : dict
        a dictionary with keys [required]
        delay (delay to start of pulse train, msec [all]
        duration: duration of pulses in train, msec [all]
        Sfreq: stimulus train frequency (Hz) [timedSpikes]
        PT: post-train test delay [all]
        NP: number of pulses in the train [timedSpikes, exp]
        amp: amplitude of the pulses in the train [all]
        hypamp: amplitude of prehyperpolarizing pulse [hyp]
        hypdur: duration of prehyperpolarizing pulse [hyp]
        spikeTimes" times of spikes [timedSpikes]

    dt : time (microseconds) [required]
        step time, in microseconds. Required parameter

    pulsetype : string (default: 'square')
        Type of pulse to generate: one of square, hyp, timedspikes or exp
        square produces a train of "square" (retangular) pulses levels 0 and ampitude
        hyp is like square, but precedes the pulse train with a single prepulse
        of hypamp and hypdur
        timedspikes is like square, excpet the pulses are generated at times specified
        in the spikeTimes key in the stim dictionary
        exp: pulses with an exponential decay.
    
    TO DO: 
        add pulsetypes, including sine wave, rectified sine wave, etc.

    Returns
    -------
    list containing [waveform (numpy array),
                    maxtime(float),
                    timebase (numpy array)]
    """
    
    assert dt is not None
    assert 'delay' in stim.keys()
    delay = int(np.floor(stim['delay']/dt))
    if pulsetype in ['square', 'hyp', 'exp']:
        ipi = int(np.floor((1000.0/stim['Sfreq'])/dt))
    pdur = int(np.floor(stim['duration']/dt))
    posttest = int(np.floor(stim['PT']/dt))
    if pulsetype not in ['timedSpikes']:
        NP = int(stim['NP'])
    else:
        NP = len(stim['spikeTimes'])
    tstims = [0]*NP
    if pulsetype == 'hyp':
        assert 'hypamp' in stim.keys()
        assert 'hypdur' in stim.keys()
        hypdur = int(np.floor(stim['hypdur']/dt))
        delay0 = delay # save original delay


    if pulsetype in ['square', 'hyp']:
        maxt = dt*(stim['delay'] + (ipi*(NP+2)) + posttest + pdur*2)
        if pulsetype == 'hyp':
            maxt = maxt + dt*stim['hypdur']
            delay = delay + hypdur
        w = np.zeros(int(np.floor(maxt/dt)))
        for j in range(0, NP):
            t = (delay + j *ipi)*dt
            w[delay+ipi*j:delay+(ipi*j)+pdur] = stim['amp']
            tstims[j] = delay+ipi*j
        if stim['PT'] > 0.0:
            send = delay+ipi*j
            for i in range(send+posttest, send+posttest+pdur):
                w[i] = stim['amp']
        if pulsetype == 'hyp': # fill in the prepulse now
            for i in range(delay0, delay0+hypdur):
                w[i] = stim['hypamp']

    if pulsetype == 'timedSpikes':
        maxt = np.max(stim['spikeTimes'])+stim['PT']+stim['duration']*2
        w = np.zeros(int(np.floor(maxt/dt)))
        for j in range(len(stim['spikeTimes'])):
            st = delay + int(np.floor(stim['spikeTimes'][j]/dt))
            t = st*dt
            w[st:st+pdur] = stim['amp']
            tstims[j] = st

        if stim['PT'] > 0.0:
            for i in range(st+posttest, st+posttest+pdur):
                w[i] = stim['amp']

    if pulsetype == 'exp':
        maxt = dt*(stim['delay'] + (ipi*(NP+2)) + posttest + pdur*2)
        w = np.zeros(int(np.floor(maxt/dt)))

        for j in range(0, NP):
            for i in range(0, len(w)):
                if delay+ipi*j+i < len(w):
                    w[delay+ipi*j+i] = w[delay+ipi*j+i] + stim['amp']*((1.0-np.exp(-i/(pdur/3.0)))*
                    np.exp(-(i-(pdur/3.0))/pdur))
            tstims[j] = delay+ipi*j
        if stim['PT'] > 0.0:
            send = delay+ipi*j
            for i in range(send+posttest, len(w)):
                w[i] += stim['amp']*(1.0-np.exp(-i/(pdur/3.0)))*np.exp(-(i-(pdur/3.0))/pdur)

    return(w, maxt, tstims)
