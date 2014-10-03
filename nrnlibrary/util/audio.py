"""
Tools for generating auditory stimuli. 

Adapted from Manis; makeANF_CF_RI.m

"""
from __future__ import division
import numpy as np

def modtone(t, rt, Fs, F0, dBSPL, FMod, DMod, phaseshift):
    """
    Apply modulation and ramps to a waveform.
    
    FMod : modulation frequency
    DMod : modulation depth percent
    phaseshift : modulation phase
    """
    """
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
    
    pin = (np.sqrt(2) * 20e-6 * 10**(dBSPL/20.)) * np.sin((2*pi*F0*t) - np.pi/2) * env # unramped stimulus
    pin = ramp(pin, mxpts, irpts)
    env = ramp(env, mxpts, irpts)
    return pin, env


def ramp(pin, mxpts, irpts):
    """
    Apply linear ramps to *pin*.
    """
    """
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

def pipnoise(t, rt, Fs, dBSPL, pip_dur, pip_start):
    """
    function [pin, mask] = pipnoise(t, rt, Fs, dBSPL, pip_dur, pip_start)
        % note factor of 2 in the scaling here: this is to make the signal level
        % match the "100%" modulation of the modtone above, with the same SPL.
        irpts = rt*Fs;
        mxpts = length(t);
        pin = sqrt(2)*20e-6*10^(dBSPL/20)*randn(size(t))*2; % unramped stimulus
        mask = zeros(size(t));
        for i = 1:length(pip_start)
            ts = floor(pip_start(i)*Fs);
            te = floor((pip_start(i)+pip_dur)*Fs);
            if te > mxpts
                te = mxpts;
            end;
            mask(ts:te) = 1;
            mask(ts:ts+irpts-1)=sin(0.5*pi*((0:(irpts-1))/irpts)).^2;
            mask((te-irpts):te)=sin(pi/2-0.5*pi*((0:(irpts))/irpts)).^2;

        end
        pin = pin .* mask;
        return;
    end
    """
    irpts = rt * Fs
    mxpts = len(t)
    pin = np.sqrt(2) * 20e-6 * 10^(dBSPL/20) * np.random.randn(len(t)) * 2  # unramped stimulus
    mask = np.zeros(len(t))
    for i in range(len(pip_start)):
        ts = np.floor(pip_start[i] * Fs)
        te = np.floor((pip_start[i] + pip_dur) * Fs)
        if te > mxpts:
            te = mxpts
        mask[ts:te] = 1
        ramp = np.linspace(0, 1, irpts)
        mask[ts:ts+irpts-1] = np.sin(0.5*np.pi * ramp)**2
        mask[te-irpts:te] = np.sin(pi/2 - 0.5*np.pi * ramp)**2
    pin *= mask
    return pin
        
   
def piptone(t, rt, Fs, F0, dBSPL, pip_dur, pip_start)
    """
    function [pin, mask] = piptone(t, rt, Fs, F0, dBSPL, pip_dur, pip_start)
        % note factor of 2 in the scaling here: this is to make the signal level
        % match the "100%" modulation of the modtone above, with the same SPL.
        irpts = rt*Fs;
        mxpts = length(t);
        pin = sqrt(2)*20e-6*10^(dBSPL/20)*(sin((2*pi*F0*t)-pi/2))*2; % unramped stimulus
        mask = zeros(size(t));
        for i = 1:length(pip_start)
            ts = floor(pip_start(i)*Fs);
            te = floor((pip_start(i)+pip_dur)*Fs);
            if te > mxpts
                te = mxpts;
            end;
            mask(ts:te) = 1;
            mask(ts:ts+irpts-1)=sin(0.5*pi*((0:(irpts-1))/irpts)).^2;
            mask((te-irpts):te)=sin(pi/2-0.5*pi*((0:(irpts))/irpts)).^2;

        end
        pin = pin .* mask;
        return;
    end
    """
    irpts = rt * Fs
    mxpts = len(t)
    pin = np.sqrt(2) * 20e-6 * 10^(dBSPL/20) * (np.sin((2*pi*F0*t) - np.pi/2))*2  # unramped stimulus
    mask = np.zeros(len(t))
    for i in range(len(pip_start)):
        ts = np.floor(pip_start[i] * Fs)
        te = np.floor((pip_start[i] + pip_dur) * Fs)
        if te > mxpts:
            te = mxpts
        mask[ts:te] = 1
        ramp = np.linspace(0, 1, irpts)
        mask[ts:ts+irpts-1] = np.sin(0.5*np.pi * ramp)**2
        mask[te-irpts:te]  = np.sin(pi/2 - 0.5*np.pi * ramp)**2

    pin *= mask
    return pin
    