from scipy import interpolate
import numpy as np


def find_point(x, y, peakindex, val, direction='left', limits=None):
    """
    Given a waveform defined by *x* and *y* arrays, return the first time
    at which the waveform crosses (y[peakindex] * val). The search begins at
    *peakindex* and proceeds in *direction*. 
    
    Optionally, *limits* may specify a smaller search region in the form of
    (t0, t1, dt).
    """
    #F = interpolate.UnivariateSpline(x, y, s=0) # declare function
    #To find x at y then do:
    istart = 0
    iend = len(y)
    if limits is not None:
        istart = int(limits[0] / limits[2])
        iend = int(limits[1] / limits[2])
    yToFind = y[peakindex] * val
    if direction == 'left':
        yreduced = np.array(y[istart:peakindex]) - yToFind
        try:
            Fr = interpolate.UnivariateSpline(x[istart:peakindex], yreduced, s=0)
        except:
            print 'find_point: insufficient time points for analysis'
            print 'arg lengths:', len(x[istart:peakindex]), len(yreduced)
            print 'istart, peakindex: ', istart, peakindex
            print 'ytofine: ', yToFind
            raise
            res = float('nan')
            return (res)
        res = Fr.roots()
        if len(res) > 1:
            res = res[-1]
    else:
        yreduced = np.array(y[peakindex:iend]) - yToFind
        try:
            Fr = interpolate.UnivariateSpline(x[peakindex:iend], yreduced, s=0)
        except:
            print 'find_point: insufficient time points for analysis?'
            print 'arg lengths:', len(x[peakindex:iend]), len(yreduced)
            raise
            res = float('nan')
            return (res)
        res = Fr.roots()
        if len(res) > 1:
            res = res[0]
            #pdb.set_trace()
    try:
        res.pop()
    except:
        pass
    if not res: # tricky - an empty list is False, but does not evaluate to False
        res = float('nan') # replace with a NaN
    else:
        res = float(res) # make sure is just a simple number (no arrays)
    return (res)
