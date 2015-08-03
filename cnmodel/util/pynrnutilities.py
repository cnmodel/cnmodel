#!/usr/bin/python
#
# utilities for NEURON, in Python
# Module neuron for cnmodel

import numpy as np
import numpy.ma as ma # masked array
import re, sys, gc, collections

import neuron


_mechtype_cache = None
def all_mechanism_types():
    """Return a dictionary of all available mechanism types.
    
    Each dictionary key is the name of a mechanism and each value is
    another dictionary containing information about the mechanism::
    
        mechanism_types = {
            'mech_name1': {
                'point_process': bool,
                'artificial_cell': bool,
                'netcon_target': bool,
                'has_netevent': bool,
                'internal_type': int,
                'globals': {name:size, ...},
                'parameters': {name:size, ...},
                'assigned': {name:size, ...},
                'state': {name:size, ...},
            },
            'mech_name2': {...},
            'mech_name3': {...},
            ...
        }
    
    * point_process: False for distributed mechanisms, True for point 
        processes and artificial cells.
    * artificial_cell: True for artificial cells, False otherwise
    * netcon_target: True if the mechanism can receive NetCon events
    * has_netevent: True if the mechanism can emit NetCon events
    * internal_type: Integer specifying the NEURON internal type index of 
        the mechanism
    * globals: dict of the name and vector size of the mechanism's global
        variables
    * parameters: dict of the name and vector size of the mechanism's 
        parameter variables
    * assigned: dict of the name and vector size of the mechanism's 
        assigned variables
    * state: dict of the name and vector size of the mechanism's state
        variables


    Note: The returned data structure is cached; do not modify it.
        
    For more information on global, parameter, assigned, and state 
    variables see:
    http://www.neuron.yale.edu/neuron/static/docs/help/neuron/nmodl/nmodl.html
    """
    global _mechtype_cache
    if _mechtype_cache is None:
        _mechtype_cache = collections.OrderedDict()
        mname = neuron.h.ref('')
        # Iterate over two mechanism types (distributed, point/artificial)
        for i in [0, 1]:
            mt = neuron.h.MechanismType(i)
            nmech = int(mt.count())
            # Iterate over all mechanisms of this type
            for j in range(nmech):
                mt.select(j)
                mt.selected(mname)
                
                # General mechanism properties
                name = mname[0]  # convert hoc string ptr to python str
                
                desc = {
                    'point_process': bool(i),
                    'netcon_target': bool(mt.is_netcon_target(j)),
                    'has_netevent': bool(mt.has_net_event(j)),
                    'artificial_cell': bool(mt.is_artificial(j)),
                    'internal_type': int(mt.internal_type()),
                }
                
                # Collect information about 4 different types of variables
                for k,ptype in [(-1, 'globals'), (1, 'parameters'), 
                                (2, 'assigned'), (3, 'state')]:
                    desc[ptype] = {} # collections.OrderedDict()
                    ms = neuron.h.MechanismStandard(name, k)
                    for l in range(int(ms.count())):
                        psize = ms.name(mname, l)
                        pname = mname[0]  # parameter name
                        desc[ptype][pname] = int(psize)
                
                # Assemble everything in one place
                _mechtype_cache[name] = desc
            
    return _mechtype_cache


def reset():
    """Introspect the NEURON kernel to verify that no objects are left over
    from previous simulation runs.
    """
    # Release objects held by an internal buffer
    # See https://www.neuron.yale.edu/phpBB/viewtopic.php?f=2&t=3221
    neuron.h.Vector().size()    
    
    # Make sure nothing is hanging around in an old exception or because of
    # reference cycles 
    sys.exc_clear()
    gc.collect()
    
    remaining = []
    
    # No sections left
    n = len(list(neuron.h.allsec()))
    if n > 0:
        remaining.append((n, 'Section'))
        
    n = len(neuron.h.List('NetCon'))
    if n > 0:
        remaining.append((n, 'NetCon'))
    
    # No point processes or artificial cells left
    for name, typ in all_mechanism_types().items():
        if typ['artificial_cell'] or typ['point_process']:
            n = len(neuron.h.List(name))
            if n > 0:
                remaining.append((n, name))
        
    if len(remaining) > 0:
        msg = ("Cannot reset--old objects have not been cleared: %s" % 
               ', '.join(['%d %s' % rem for rem in remaining]))
        raise RuntimeError(msg)



# routine to convert conductances from nS as given elsewhere
#   to mho/cm2 as required by NEURON 1/28/99 P. Manis
# units: nano siemens, soma area in um^2
#
def nstomho(ns, somaarea, refarea = None):
    if refarea == None:
        return 1E-9*float(ns)/float(somaarea)
    else:
        return 1e9*float(ns)/float(refarea)

def mho2ns(mho, somaarea):
    return float(mho)*somaarea/1E-9

def spherearea(dia):
    """
    given diameter in microns, return sphere area in cm2
    """
    r = dia*1e-4 # convert to cm
    return(4*np.pi*r**2)

def get_sections(h):
    """
    go through all the sections and find the names of the sections and all of their
    parts (ids). Returns a dict, of sec: [id0, id1...]

    """
    secnames = {}
    resec = re.compile('(\w+)\[(\d*)\]')
    for sec in h.allsec():
        g = resec.match(sec.name())
        if g.group(1) not in secnames.keys():
            secnames[g.group(1)] = [int(g.group(2))]
        else:
            secnames[g.group(1)].append(int(g.group(2)))
    return secnames

def all_objects():
    """ Return a dict of all objects known to NEURON. 
    
    Keys are 'Section', 'Segment', 'Mechanism', 'Vector', 'PointProcess', 
    'NetCon', ...
    """
    objs = {}
    objs['Section'] = list(h.all_sec())
    objs['Segment'] = []
    for sec in objs['Section']:
        objs['Segment'].extend(list(sec.allseg()))
    objs['PointProcess'] = []
    for seg in objs['Segment']:
        objs['PointProcess'].extend(list(seg.point_processes()))
        
    return objs
    


def alpha(alpha=0.1, delay=1, amp = 1.0, tdur = 50., dt=0.010):
    tvec = np.arange(0, tdur , dt)
    aw = np.zeros(tvec.shape)
    i = 0
    for t in tvec:
        if t > delay:
            aw[i] = amp * (t-delay)*(1./alpha) * np.exp(-(t-delay)/alpha) # alpha waveform time course
        else:
            aw[i] = 0.
        i += 1
    return(aw, tvec)

def syns(alpha=0.1, rate=10, delay=0, dur=50, amp=1.0, dt=0.020, N=1, mindur = 120, makewave=True):
    """	Calculate a poisson train of alpha waves
    with mean rate rate, with a delay and duration (in mseco) dt in msec.
    N specifies the number of such independent waveforms to sum """
    deadtime = 0.7
    if dur + delay < mindur:
        tvec = np.arange(0.0, mindur , dt)
    else:
        tvec = np.arange(0.0, dur+delay , dt)
    npts = len(tvec)
    ta = np.arange(0.0, 20.0, dt)
    aw = ta * alpha* np.exp(-ta/alpha)/alpha # alpha waveform time course
    spt = [[]]*N # list of spike times
    wave = np.array([]) # waveform
    sptime=[]
    for j in range(0,N):
        done = False
        t=0.0
        nsp = 0
        while not done:
            a = np.random.sample(1)
            if t < delay:
                t = delay
                continue
            if t >= delay and t <= (delay+dur):
                ti = -np.log(a)/(rate/1000.0) # convert to exponential distribution with rate
                if ti < deadtime:
                    continue
                t = t + ti # running time
            if t > delay+dur:
                done = True
                continue
            if nsp is 0:
                sptime = t
                nsp = nsp+1
            else:
                sptime = np.append(sptime, t)
                nsp = nsp+1
        if j is 0:
            wavej = np.zeros(len(tvec))
        for i in range(0,len(sptime)):
            st = int(sptime[i]/dt)
            wavej[st] = wavej[st] + 1
        spt[j] = sptime

        if makewave:
            w = np.convolve(wavej, aw/max(aw))*amp
            if len(w) < npts:
                w = np.append(w, np.zeros(npts-len(w)))
            if len(w) > npts:
                w = w[0:npts]
            if j is 0:
                wave = w
            else:
                wave = wave + w
    return (spt, wave, tvec, N)


def an_syn(alpha=0.1, spont=10, driven=100, delay=50, dur=100, post = 20,
           amp=0.1, dt=0.020, N=1, makewave=True):
    # constants for AN:
    deadtime = 0.7 # min time between spikes, msec
    trise = 0.2 # rise rate, ms
    tfall =  0.5 # fall rate, ms
    rss = driven/1000.0 # spikes/millisecond
    rr = 3 * rss # transient driven rate
    rst = rss
    taur = 3 # rapid decay, msec
    taust = 10
    ton = delay # msec
    stim_end = ton + dur
    trace_end = stim_end + post
    tvec = np.arange(0.0, trace_end , dt) # dt is in msec, so tvec is in milliseconds
    ta = np.arange(0.0, 20.0, dt)
    aw = ta * alpha* np.exp(-ta/alpha)/alpha # alpha waveform time course
    spt = [[]]*N # list of spike times
    wave = [[]]*N # waveform
    for j in range(0,N): # for each
        done = False
        sptime = []
        qe = 0
        nsp = 0
        i = int(0)
        if spont <= 0:
            q = 1e6
        else:
            q = 1000.0/spont # q is in msec (spont in spikes/second)
        t = 0.0
        while not done:
            a = np.random.sample(1)
            if t < ton:
                if spont <= 0:
                    t = ton
                    continue
                ti = - (np.log(a)/(spont/1000.0)) # convert to msec
                if ti < deadtime: # delete intervals less than deadtime
                    continue
                t = t + ti
                if t > ton: # if the interval would step us to the stimulus onset
                    t = ton # then set to to the stimulus onset
                    continue
            if t >= ton and t < stim_end:
                if t > ton:
                    rise = 1.0 - np.exp(-(t-ton)/trise)
                else:
                    rise = 1.0
                ra = rr*np.exp(-(t-ton)/taur)
                rs = rst*np.exp(-(t-ton)/taust)
                q = rise*(ra+rs+rss)
                ti = - np.log(a)/(q + spont/1000) # random.negexp(1000/q)
                if ti < deadtime:
                    continue
                t = t + ti
                if t > stim_end: # only include interval if it falls inside window
                    t = stim_end
                    continue
            if t >= stim_end and t <= trace_end:
                if spont <= 0.0:
                    t = trace_end
                    continue
                if qe is 0: # have not calculated the new qe at end of stimulus
                    rise = 1.0 - np.exp(-(stim_end-ton)/trise)
                    ra = rr*np.exp(-(stim_end-ton)/taur)
                    rs = rst*np.exp(-(stim_end-ton)/taust)
                    qe = rise*(ra+rs+rss) # calculate the rate at the end
                fall = np.exp(-(t-stim_end)/tfall)
                q = qe * fall
                ti = -np.log(a)/(q + spont/1000.) # keeps rate from falling below spont rate
                if ti < deadtime:
                    continue
                t = t + ti
            if t >= trace_end:
                done = True
                continue
            # now add the spike time to the list
            if nsp is 0:
                sptime = t
                nsp = nsp+1
            else:
                sptime = np.append(sptime, t)
                nsp = nsp+1
        # end of for loop on i
        if j is 0:
            wavej = np.zeros(len(tvec))
        for i in range(0,len(sptime)):
            st = int(sptime[i]/dt)
            wavej[st] = wavej[st] + 1
        spt[j] = sptime
        npts = len(tvec)
        if makewave:
            w = np.convolve(wavej, aw/max(aw))*amp
            wave[j] = w[0:npts]
    return (spt, wave, tvec, N)

def findspikes(t, v, thresh):
    """ findspikes identifies the times of action potential in the trace v, with the
    times in t. An action potential is simply timed at the first point that exceeds
    the threshold. """
    tm = np.array(t)
    s0 = np.array(v) > thresh # find points above threshold
    dsp = tm[s0]
    sd = np.append(True, np.diff(dsp) > 1.0) # find first points of spikes
    if len(dsp) > 0:
        sp = dsp[sd]
    else:
        sp = []
    return(sp) # list of spike times.

def measure(mode, x, y, x0, x1):
    """ return the mean and standard deviation of y in the window x0 to x1
    """
    xm = ma.masked_outside(x, x0, x1)
    ym = ma.array(y, mask = ma.getmask(xm))
    if mode == 'mean':
        r1 = ma.mean(ym)
        r2 = ma.std(ym)
    if mode == 'max':
        r1 = ma.max(ym)
        r2 = 0
    if mode == 'min':
        r1 = ma.min(ym)
        r2 = 0
    if mode == 'median':
        r1 = ma.median(ym)
        r2 = 0
    if mode == 'p2p': # peak to peak
        r1 = ma.ptp(ym)
        r2 = 0
    return(r1, r2)

def mask(x, xm, x0, x1):
    xmask = ma.masked_outside(xm, x0, x1)
    xnew = ma.array(x, mask=ma.getmask(xmask))
    return(xnew.compressed())

def vector_strength(spikes, freq):
    """
    Calculate vector strength and related parameters from a spike train, for the specified frequency
    :param spikes: Spike train, in msec.
    :param freq: Stimulus frequency in Hz
    :return: a dictionary containing:
        r: vector strength
        n: number of spikes
        R: Rayleigh coefficient
        p: p value (is distribution not flat?)
        ph: the circularized spike train over period of the stimulus freq, freq, in radians
        d: the "dispersion" computed according to Ashida et al., 2010, etc.
    """
    per = 1e3/freq # convert from Hz to period in msec
    ph = 2*np.pi*np.fmod(spikes, per)/(per) # convert to radians within a cycle
    c = np.sum(np.cos(ph))**2
    s = np.sum(np.sin(ph))**2
    vs = (1./len(ph))*np.sqrt(c+s)  # standard vector strength computation
    n = len(spikes)
    R = n*vs  # Raleigh coefficient
    Rp = np.exp(-n*vs*vs)  # p value for n > 50 (see Ashida et al. 2010).
    d = np.sqrt(2.*(1-vs))/(2*np.pi*freq)
    return{'r': vs, 'n': n, 'R': R, 'p': Rp, 'ph': ph, 'd': d}

def isi_cv2(splist, binwidth=1, t0=0, t1=300, tgrace=25):
    """ compute the cv and regularity according to Young et al., J. Neurophys, 60: 1, 1988.
        Analysis is limited to isi's starting at or after t0 but before t1, and ending completely
        before t1 + tgrace(to avoid end effects). t1 should correspond to the
        the end of the stimulus
        VERSION using dictionary for cvisi
    """
    cvisit = np.arange(0, t1, binwidth) # build time bins
    cvisi = {} # isi is dictionary, since each bin may have different length
    for i in range(0, len(splist)): # for all the traces
        isit = splist[i] # get the spike times for this trial [1:-1]
        if len(isit) <= 1: # need 2 spikes to get an interval
            continue
        isib = np.floor(isit[0:-2]/binwidth) # discreetize
        isii = np.diff(splist[i]) # isis.
        for j in range(0, len(isib)): # loop over possible start time bins
            if isit[j] < t0 or isit[j] > t1 or isit[j+1] > t1+tgrace: # start time and interval in the window
                continue
            if isib[j] in cvisi:
                print "spike in bin: %d" % (isib[j])
                cvisi[isib[j]] = np.append(cvisi[isib[j]], isii[j]) # and add the isi in that bin
            else:
                cvisi[isib[j]] = isii[j]	# create it
    cvm = np.array([]) # set up numpy arrays for mean, std and time for cv analysis
    cvs = np.array([])
    cvt = np.array([])
    for i in cvisi.keys(): # for each entry (possible bin)
        c = [cvisi[i]]
        s = c.shape
        print c
        if len(s) > 1 and s[1] >= 3: # require 3 spikes in a bin for statistics
            cvm = np.append(cvm, np.mean(c))
            cvs = np.append(cvs, np.std(c))
            cvt = np.append(cvt, i*binwidth)
    return(cvisit, cvisi, cvt, cvm, cvs)

def isi_cv(splist, binwidth=1, t0=0, t1=300, tgrace=25):
    """ compute the cv and regularity according to Young et al., J. Neurophys, 60: 1, 1988.
        Analysis is limited to isi's starting at or after t0 but before t1, and ending completely
        before t1 + tgrace(to avoid end effects). t1 should correspond to the
        the end of the stimulus
        Version using a list of numpy arrays for cvisi
    """
    cvisit = np.arange(0, t1, binwidth) # build time bins
    cvisi = [[]]*len(cvisit)
    for i in range(0, len(splist)): # for all the traces
        if len(splist[i]) < 2: # need at least 2 spikes
            continue
        isib = np.floor(splist[i][0:-2]/binwidth) # begining spike times for each interval
        isii = np.diff(splist[i]) # associated intervals
        for j in range(0, len(isib)): # loop over spikes
            if splist[i][j] < t0 or splist[i][j] > t1 or splist[i][j+1] > t1+tgrace: # start time and interval in the window
                continue
            cvisi[int(isib[j])] = np.append(cvisi[int(isib[j])], isii[j]) # and add the isi in that bin
    cvm = np.array([]) # set up numpy arrays for mean, std and time for cv analysis
    cvs = np.array([])
    cvt = np.array([])
    for i in range(0, len(cvisi)): # for each entry (possible bin)
        c = cvisi[i]
        if len(c) >= 3: # require 3 spikes in a bin for statistics
            cvm = np.append(cvm, np.mean(c))
            cvs = np.append(cvs, np.std(c))
            cvt = np.append(cvt, i*binwidth)
    return(cvisit, cvisi, cvt, cvm, cvs)

if __name__ == "__main__":

    test = 'isicv'

    if test == 'isicv':
        """ this test is not perfect. Given an ISI, we calculate spike times
        by drawing from a normal distribution whose standard deviation varies
        with time, from 0 (regular) to 1 (irregular). As coded, the standard
        deviation never reaches the target value because spikes fall before or
        after previous spikes (thus reducing the stdev). Nonetheless, this shows
        that the CV calculation works correctly. """
        nreps = 500
        # cv will be 0 for first 50 msec, 0.5 for next 50 msec, and 1 for next 50 msec
        d=[[]]*nreps
        isi = 5.0 # mean isi
# we create 100 msec of data where the CV goes from 0 to 1
        maxt = 100.0
        for i in range(nreps):
            for j in range(int(maxt/isi)+1):
                t = float(j) * isi
                sd = float(j)/isi
                if sd == 0.0:
                    d[i] = np.append(d[i], t)
                else:
                    d[i] = np.append(d[i], np.random.normal(t, sd, 1) )
            for j in range(1, 10): # add more intervals at the end
                te = t + float(j)*isi
                d[i] = np.append(d[i], np.random.normal(te, sd, 1))
            d[i] = np.sort(d[i])
#			print d[i]
#			print diff(d[i])
        sh = np.array([])
        for i in  range(len(d)):
            sh = np.append(sh, np.array(d[i]))
        (hist, bins) = np.histogram(sh, bins = 250, range=(0, 250), new=True)
        if len(bins) > len(hist):
            bins=bins[0:len(hist)]

        pl.figure(1)
        pl.subplot(2,2,1)
        pl.plot(bins, hist)

        (cvisit, cvisi, cvt, cvm, cvs) = isi_cv(d, binwidth = 0.5, t0=0, t1=100, tgrace = 25)
        order = np.argsort(cvt)
        cvt = cvt[order]
        cvs = cvs[order]
        cvm = cvm[order]
        pl.subplot(2,2,2)
        pl.plot(cvt, cvm)
        pl.hold(True)
        pl.plot(cvt, cvs)
        pl.subplot(2,2,4)
        pl.plot(cvt, cvs/cvm)
        pl.show()

    if test == 'measure':
        x = np.arange(0, 100, 0.1)
        s = np.shape(x)
        y = np.random.randn(s[0])
        for i in range(0,4):
            print "\ni is : %d" % (i)
            x0 = i*20
            x1 = x0 + 20
            (r0, r1) = measure('mean', x, y, x0, x1)
            print 'mean: %f   std: %f [0, 20]' % (r0, r1)
            (r0, r1) = measure('max', x, y, x0, x1)
            print 'max: %f   std: %f [0, 20]' % (r0, r1)
            (r0, r1) = measure('min', x, y, x0, x1)
            print 'min: %f   std: %f [0, 20]' % (r0, r1)
            (r0, r1) = measure('median', x, y, x0, x1)
            print 'median: %f   std: %f [0, 20]' % (r0, r1)
            (r0, r1) = measure('p2p', x, y, x0, x1)
            print 'peak to peak: %f   std: %f [0, 20]' % (r0, r1)

    if test == 'an_syn':
        (s,w,t,n) = an_syn(N=50, spont=50, driven=150, post=100, makewave=True)
        sh = np.array([])
        for i in  range(len(s)):
            sh = np.append(sh, np.array(s[i]))
        (hist, bins) = np.histogram(sh, bins = 250, range=(0, 250), new=True)
        if len(bins) > len(hist):
            bins=bins[0:len(hist)]

        import pylab as pl
        pl.figure(1)
        pl.subplot(2,2,1)
        pl.plot(bins, hist)

        pl.subplot(2,2,3)
        for i in range(len(w)):
            pl.plot(t, w[i])
            pl.hold = True
        (cvisit, cvisi, cvt, cvm, cvs) = isi_cv(s)
        order = np.argsort(cvt)
        cvt = cvt[order]
        cvs = cvs[order]
        cvm = cvm[order]
        pl.subplot(2,2,2)
        pl.plot(cvt, cvs/cvm)
        pl.show()

    if test == 'syns':
        (s,w,t,n) = syns(rate=20, delay=0, dur=100.0, N=5, makewave=True)
        sh = np.array([])
        for i in  range(len(s)):
            sh = np.append(sh, np.array(s[i]))
        (hist, bins) = np.histogram(sh, bins = 250, range=(0, 250), new=True)
        if len(bins) > len(hist):
            bins=bins[0:len(hist)]

        import pylab as pl
        pl.figure(1)
        pl.subplot(2,2,1)
        pl.plot(bins, hist)

        pl.subplot(2,2,3)
        pl.plot(t, w)
        pl.hold = True
        (cvisit, cvisi, cvt, cvm, cvs) = isi_cv(s)
        order = np.argsort(cvt)
        cvt = cvt[order]
        cvs = cvs[order]
        cvm = cvm[order]
        pl.subplot(2,2,2)
        pl.plot(cvt, cvs/cvm)
        pl.show()
