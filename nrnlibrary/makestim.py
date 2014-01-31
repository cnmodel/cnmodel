__author__ = 'pbmanis'

import numpy as np


def makestim(stim, pulsetype = 'square', dt=None):
    """
        Create the stimulus pulse waveform.
        Inputs:
            pulsetype: one of square, timedspikes or exp
                square produces a train of "square" (retangular) pulses levels 0 and ampitude
                timedspikes is like square, excpet the pulses are generated at times specified
                    in the spikeTimes key in the stim dictionary
                exp: pulses with an exponential decay.
            stim: a dictionary with keys:
                delay (delay to start of pulse train, msec
                duration: duration of pulses in train, msec
                SFreq: stimulus train frequency (Hz)
                PT: post-train test delay
                NP: number of pulses in the train
                amp: amplitude of the pulses in the train

        Outputs:
        list containing [waveform (numpy array),
                        maxtime(float),
                        timebase (numpy array)]
        Side Effects: None
    """
    assert dt is not None
    assert 'delay' in stim.keys()
    delay = int(np.floor(stim['delay']/dt))
    ipi = int(np.floor((1000.0/stim['Sfreq'])/dt))
    pdur = int(np.floor(stim['dur']/dt))
    posttest = int(np.floor(stim['PT']/dt))
    NP = int(stim['NP'])

    #   make pulse
    tstims = [0]*NP
    if pulsetype == 'square':
        maxt = dt*(stim['delay'] + (ipi*(NP+2)) + posttest + pdur*2)
        w = np.zeros(np.floor(maxt/dt))
        for j in range(0, NP):
            t = (delay + j *ipi)*dt
            w[delay+ipi*j:delay+(ipi*j)+pdur] = stim['amp']
            tstims[j] = delay+ipi*j
        if stim['PT'] > 0.0:
            send = delay+ipi*j
            for i in range(send+posttest, send+posttest+pdur):
                w[i] = stim['amp']

    if pulsetype == 'timedSpikes':
        maxt = np.max(stim['spikeTimes'])+stim['PT']+stim['dur']*2
        w = np.zeros(np.floor(maxt/dt))
        print 'makestim: spike times max: ', maxt
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

        for j in range(0, NP):
            for i in range(0, len(w)):
                if delay+ipi*j+i < len(w):
                    w[delay+ipi*j+i] += stim['amp']*(1.0-np.exp(-i/(pdur/3.0)))*np.exp(-(i-(pdur/3.0))/pdur)
            tstims[j] = delay+ipi*j
        if stim['PT'] > 0.0:
            send = delay+ipi*j
            for i in range(send+posttest, len(w)):
                w[i] += stim['amp']*(1.0-np.exp(-i/(pdur/3.0)))*np.exp(-(i-(pdur/3.0))/pdur)

    return(w, maxt, tstims)
