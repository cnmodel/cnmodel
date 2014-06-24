#!/usr/bin/env python

import numpy
import matplotlib.pylab as MP
from optparse import OptionParser


class Table():
	def __init__(self):
		self.F = 0.4 # release probability (constant; no facilitiation)
		self.k0 = 1/1.75 # /s, baseline recovery rate from depletion (slow rate)
		self.kmax = 1/0.025 # /s, maximal recovery rate from depletion (fast rate)

		self.td = 0.05 #  time constant for calcium-dependent recovery
		self.kd =  0.7 # affinity of fast recovery proces for calcium sensor

		self.ts = 0.015 # decay time constant of glutatme clearance
		self.ks = 1000 # affinity of receptor desensitization for glutatmate
		# The large value means no desense occurs (otherwise, ks should be about
		# 0.6)

		self.kf = 0.6 # affinity of facilitiation process
		self.tf = 0.01 # make facilitation VERY slow

		self.dD = 0.02 # sets Ca that drives recovery(Ca influx per AP)
		# 0.02 yields rate-dep recovery in 100-300 Hz
		self.dF = 0.02 # sets Ca that drives facilitation

		self.glu = 0.3
		return


	def bushy_epsc(self):
		""" data is average of 3 cells studied with recovery curves and individually fit """
		self.F = 0.29366
		self.k0 = 0.52313
		self.kmax = 19.33805
		self.td = 0.01516
		self.kd = 0.11283
		self.ts = 17.9122
		self.ks = 11.531
		self.kf = 17.78
		self.tf = 0.00975
		self.dD = 0.57771
		self.dF = 0.60364
		self.glu = 2.12827
		return(self)


	def stellate_epsc(self):
		""" data is average of 2 cells studied with recovery curves and individually fit """
		self.F = 0.30514
		self.k0 = 0.04045
		self.kmax = 38.4408
		self.td = 0.00549
		self.kd = 0.11911
		self.ts = 23.626
		self.ks = 21.06690
		self.kf = 12.24439
		self.tf = 0.01256
		self.dD = 2.07174
		self.dF = 2.08471
		self.glu = 8.298
		return(self)

	def bushy_ipsc(self):
		""" data is average of 4 cells studied with recovery curves and individually fit """
		self.F = 0.23382
		self.k0 = 0.67554
		self.kmax = 52.93832
		self.td = 0.08195
		self.kd = 0.28734
		self.ts = 0.17500
		self.ks = 4.57098
		self.kf = 16.21564
		self.tf = 0.12336
		self.dD = 2.21580
		self.dF = 1.17146
		self.glu = 1.90428
		return(self)

	def stellate_ipsc(self):
		""" data is average of 3 cells studied with recovery curves and individually fit """
		self.F = 0.23047
		self.k0 = 1.23636
		self.kmax = 45.34474
		self.td = 0.09809
		self.kd = 0.01183
		self.ts = 17.61450
		self.ks = 17.88618
		self.kf = 19.11424
		self.tf = 0.03228
		self.dD = 2.52072
		self.dF = 2.33317
		self.glu = 3.06948
		return(self)

def XuF(callmode = None, table=None, stim = None, recovery=None, plot=False, celltype = 'bushy_epsc', desense=False):
	"""
	# implementation eqs 1-4 of xu-f 2008///
	# or you could say, of Dittman et al, 2000
	# or Dittman and Regehr, 1998....
	#
	# call:
	# XuF no arguments just runs with some default parameters
	# XuF (mode, varargs):
	# if mode = 1, the table is passed in the second argument, and the plot is
	# generated.
	# if mode = 2, the table is in the second, the train in the third and the
	# recovery times in the 4th. For modeling/fitting
	# if mode = 3, the table is in the second, the train in the third and we
	# assume that there is no recovery data.
	# P. Manis 6/2008.
	# converted to Python 11/2009
	#
	"""

	colors = ['r','g','b', 'c', 'k']
# define constants: (Table 1, Yang and Xu-Friedman, 2008, modified)
	xout = []
	yout = []
	plotsvn = False
	if callmode == None:
		theTable = Table()
		if celltype == 'bushy_epsc':
			table = Table.bushy_epsc(theTable)
		if celltype == 'bushy_ipsc':
			table = Table.bushy_ipsc(theTable)
		if celltype == 'stellate_epsc':
			table = Table.stellate_epsc(theTable)
		if celltype == 'stellate_ipsc':
			table = Table.stellate_ipsc(theTable)
		freqlist = [50, 100,200,300] # in Hz
		modelist = 0 # if 0, use regular; otherwise use exponential.
		plotvsn = True # control plots
		traindur = 0.5 # sec
		recovery = numpy.array([10.0,20.,30.,40.,50.,100.,200.,300.,500.,1000.]) # recovery points
		recovery = recovery*0.001 # convert to seconds
	else:
		if callmode in (0,1):
			freqlist = [50,100,200,400] # in Hz
			if stim is not None:
				modelist = stim
			else:
				modelist = 0 # if 0, use regular; otherwise use exponential.
			if callmode is 1:
				plotvsn = True # control plots

			traindur = 0.5 # sec
			recovery = 0.001*numpy.array([10.0,20,30,40,50,100,200,300,500,1000]) # recovery points
		elif callmode == 2:# with recovery data
			xtb = stim # cell array of times for stimulation.
			xtr = recovery # and recoveries (matched)
			freqlist = numpy.zeros(len(xtb))
			for k in range(0, len(xtb)):
				freqlist[k] = 1.0/numpy.mean(numpy.diff(xtb[k]))
			modelist = 2
		elif callmode == 3: # no recovery
			xtb = stim # cell array of times for stimulation.
			xtr = [] # and recoveries (matched)
			freqlist = numpy.arange((0,len(xtb))*100)
			modelist = 2
		else:
			print ('XuF - mode not recognized\n')
			return(None)

	mode = modelist
	if plot is True:
		fig = MP.figure()
		p1 = fig.add_subplot(3,1,3)
		p2 = fig.add_subplot(6,1,1)
		p3 = fig.add_subplot(6,1,2)
		p4 = fig.add_subplot(6,1,3)
		p5 = fig.add_subplot(6,1,4)


	n = -1
	if mode == 1:
		ntrial = 10
	else:
		ntrial = 1
	for freq in freqlist:
		n = n + 1
		if mode == 0:
			mark_face = 'w' # regular are open symbols
		elif mode == 2:
			mark_face = 'k'
		else:
			mark_face = colors[n] # poisson are filled
		for nt in range(0,ntrial): # to allow averaging of poisson train trials
			# print "trial: %d" % nt
			if mode == 0:
				pt = (1.0/freq)*numpy.ones(freq*traindur) # pulse train, msec.
				pt = numpy.cumsum(pt)
				pt = 0.010+pt-pt[0] # make starts all the same.
			elif mode == 1:
				pt=[]
				pt[0] = numpy.randexpo(freq)
				k = 1
				while numpy.max(numpy.cumsum(pt)) < traindur:
					pt[k]=numpy.randexpo(freq) ##ok<AGROW>
					k = k + 1
				pt=numpy.cumsum(pt) # that's the stimulus train.
				pt = 0.010+pt-pt[0] # make starts all the same.
			elif mode == 2: # get from incoming data
				pt = xtb[n]
				if xtr is not []:
					rec = xtr[n]
				else:
					rec = 0

			else:
				pass

			ESum = numpy.array([])
			DSum = ESum
			SSum = ESum
			CaDiSum = ESum
			CaFiSum = ESum
			FnSum = ESum
			gluSum = ESum
			for j in range(0, len(recovery)):
				ptt = numpy.append(pt, numpy.amax(pt) + recovery[j]) # pick just one recovery pulse per train to prevent interactions
				D = numpy.ones(len(ptt)) # depletion (1-D)
				S = numpy.zeros(len(ptt)) # desens (1-S)
				E = numpy.zeros(len(ptt)) # EPSP
				glu = numpy.zeros(len(ptt))
				CaDi = numpy.zeros(len(ptt))
				CaFi = numpy.zeros(len(ptt))
				Fn = numpy.zeros(len(ptt))
				Fn[0] = table.F
				glu[0] = table.glu
				dt = 1 # just a very long time since the last pulse
				(D[0], CaDi[0], Fn[0], CaFi[0]) = dstep(dt, table, 1, 0.1,  table.F, 0.1)
				if desense:
					S[0] = table.ks/(table.ks + glu[0])
				else:
					S[0] = 1.0
				E[0] = Fn[0] * D[0]*S[0]
				for i in range(1, len(ptt)):
					dt = ptt[i]-ptt[i-1]
					[D[i], CaDi[i], Fn[i], CaFi[i]] = dstep(dt, table, D[i-1], CaDi[i-1], Fn[i-1], CaFi[i-1])
					glu[i] = (glu[i-1] + table.F * D[i] )* numpy.exp(-dt/table.ts)
					if desense:
						S[i] = table.ks/(table.ks + glu[i])
					else:
						S[i] = 1.0
					E[i] = Fn[i] * D[i] * S[i]
				if j == 0:
					ESum = E
					DSum = D
					SSum = S
					CaDiSum = CaDi
					CaFiSum = CaFi
					FnSum = Fn
					gluSum = glu
				else:
					ESum = numpy.append(ESum, E[-1])
					DSum = numpy.append(DSum, D[-1])
					SSum = numpy.append(SSum, S[-1])
					CaDiSum = numpy.append(CaDiSum, CaDi[-1])
					CaFiSum = numpy.append(CaFiSum, CaFi[-1])
					FnSum = numpy.append(FnSum, Fn[-1])
					gluSum = numpy.append(gluSum, glu[-1])
			if nt == 0:
				ER = ESum[len(pt):]
			else:
				ER = ER + ESum[len(pt):]
		ER = ER/ntrial

		xout = numpy.append(pt, numpy.amax(pt) + recovery)

		yout = ESum/ESum[0]

		if plot is True: # plot if requested

			p1.plot(xout, yout, colors[n] + 's-', markerfacecolor = mark_face, markersize = 1.5)
			p1.set_ylim((0, 2.5))
			p1.hold(True)
			#p1.plot(pt, E[0:-1]/esum[0],  colors[n] + 's-', markerfacecolor =mark_face, markersize =1.5)
			p1.set_ylabel('E(normalized)')
			if plotsvn:
				xp = range(0,len(D))
			else:
				xp = ptt

			p2.plot(xout, DSum, colors[n] + 'o-', markerfacecolor = mark_face, markersize = 1.5)
			p2.set_ylim(0,1)
			p2.set_ylabel('D')
			p2.hold(True)
			p3.plot(xout, FnSum, colors[n] + 'o-', markerfacecolor = mark_face, markersize = 1.5)
			p3.set_ylim(0,1)
			p3.set_ylabel('Fn')
			p3.hold(True)
			p4.plot(xout, CaFiSum, colors[n] + 'd-', markerfacecolor = mark_face, markersize = 1.5)
			p4.set_ylabel( 'CaD (d) CaF (x)')
			#    set(gca, 'Ylim', [0 1])
			p4.hold(True)
			p4.plot(xout, CaDiSum, colors[n] + 'x-', markerfacecolor = mark_face, markersize = 1.5)
			p5.plot(xout, SSum, colors[n] + 'd-', markerfacecolor = mark_face, markersize = 1.5)
			p5.set_ylabel('S')
			#    set(gca, 'Ylim', [0 1])
			p5.hold(True)

		xout = numpy.append(xout, max(pt) + recovery)
		yout = numpy.append(yout, ESum/ESum[0])
	if plot is True:
		MP.show()
	return(xout, yout)

# [CaDi'; CaFi']'
# this is the discreet step for solving the analytical equations as
# described in the papers. This version, including F, is from Dittman et
# al.
def dstep(dt, T, Di, CaDi, Fi, CaFi):
# calculate next D from Equations 15/16, Dittman et al. 2000
# dt is time since last step
# Di is D at last time step
# T is the table (parameters/constants).
# CaDi is the calcium term for depression
# CaFi is the calcium term for facilitiation

	CaDi = CaDi + T.dD
	CaFi = CaFi + T.dF
	CaDn = CaDi*numpy.exp(-dt/T.td) # change with next step
	CaFn = CaFi*numpy.exp(-dt/T.tf)
	# r = 1/dt # convert to hz (units are SECONDS).
	eta = (T.kd/CaDi + 1)/(T.kd/CaDi + numpy.exp(-dt/T.td))
	eta = numpy.power(eta,(-(T.kmax-T.k0)*T.td))
	Dn = 1-(1-(1-Fi)*Di)*numpy.exp(-T.k0*dt)*eta
	Fn = T.F + (1-T.F)/(1+T.kf/CaFn)
	return(Dn, CaDn, Fn, CaFn)

def test():
    XuF(plot=True, celltype = 'bushy_epsc')

if __name__ == "__main__":

    parser=OptionParser()
    parser.add_option("-b", "--bushyepsc", dest="be", action='store_true',
                      help="test bushy_epsc", default=False)
    parser.add_option("-s", "--stellateepsc", dest="se", action='store_true',
                      help="test stellate_epsc", default=False)
    parser.add_option("-B", "--bushyipsc", dest="bi", action='store_true',
                      help="test bushy_ipsc", default=False)
    parser.add_option("-S", "--stellateipsc", dest="si", action='store_true',
                      help="test bushy_epsc", default=False)

    (options, args) = parser.parse_args()
    model = 'bushy_epsc'
    if options.be:
        model = 'bushy_epsc'
    if options.bi:
        model = 'bushy_ipsc'
    if options.se:
        model = 'stellate_epsc'
    if options.si:
        model = 'stellate_ipsc'
    XuF(plot=True, celltype = model)
