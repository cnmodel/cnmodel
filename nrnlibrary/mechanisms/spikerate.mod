COMMENT
-----------------------------------------------------------------------------

  Spike Rate Controller and PRC measure
  ----------------------

  This mechanism was written to be used as a Spike Rate Controller that injects
  current into a cell 

  WARNING: this current is NOT inserted as a standard current, but as a 
  point process (same way as a synapse or current injection).

  Procedure for insertion:

	objectvar sr
	sr = new spikerate()

	access <compartment_name>
	sr.loc(0.5)

  or, in python:
	sr = h.spikerate(0.5, sec=cwcell), where cwccell is the target section of the cell.
  
	sr.rate = setrate (spikes per second)
	sr.thresh = 0 (mV) for spike threshold
	sr.nwincount = 5 - number of spikes to average over for window
	sr.on = 1 - to allow the rate controller to operate
	sr.on =  0 - don't let rate controller change rate (for PRC measure)
	sr.slope - rate to ramp
	sr.slopedur = time to ramp end
	sr.rdeltaplus = current increase factor to make spike when none happens...
	
we fit ongoing spike rate versus current 

  P. Manis UNC Chapel Hill, 8/2007

-----------------------------------------------------------------------------
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
DEFINE SSIZE 10001

NEURON {
	POINT_PROCESS spikerate
	RANGE on, rate, drate, thresh
	RANGE slope, slopedur, rdeltaplus
	RANGE nspikes, nwincount
	RANGE spike, curr, isi, cisi
	RANGE prcx, prcy, prccount, prcdur, prcamp
	ELECTRODE_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	on = 1 (1) : flag to enable controller to run
	rate	= 20 (1/s)		: target spike rate (/sec)
	drate = 0.2 : delta rate (fractional error before doing PRC analysis)
	thresh = 0 (mV)   : spike detection threshold, mV

	slope = 0.0005 (nA/ms)
	slopedur = 1000.0 (ms)
	rdeltaplus = 2

	nwincount = 5 (1) : number of spikes to average over for rate measure

	prcdur = 0.1
	prcamp = 0.1
	}


ASSIGNED {
    spike[SSIZE]    (ms)    : list of spike times
    curr[SSIZE] (nA) : list of currents for each spike time
	isi[SSIZE] (ms)
	cisi[SSIZE] (nA)
	prcx[50] (ms)
	prcy[50] (ms)
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current
	intvl (ms)			: holder for interspike interval
	mrate (1/s)			: calulated mean rate over nwincount spikes
	lastt (ms)			: time of last event... 
	lastv (mV)			: temporary voltage (previous voltage)
	nspikes (1)			: number of spikes detected so far
	sx (1)				: next set of variables are for the linear regression
	sy (1)
	sxy (1)
	sx2 (1)
	sy2 (1)
	dy (1)
	a (1)
	b (1)
	c (1)
	j (1)
	a1 (1)
	a2 (1)
	iset (1)			: flag to indicate that a current adjusment was made
	inew (nA)			: temporary hold for new current (? not needed)
	ilog (1)
	nprc (1)			: count of PRC intervals tested so far
	prcskip (1)			: number of spikes to skip between PRC measures
	prccount (1)		: 
	prcskipcnt (1)		: skip counter
	prcintvl (ms)		: recent interval
	thisdelay (ms)		: current PRC delay to test
	thisperturb (1)		: state variable: 0 = waiting for interval, 1 = perturb current injected and still on, 2 = perturb current injected and done (thus save interval with next spike
	tperturb (ms)		: time that perturbation started.
	debug (1)				: debug flag for verbose messages.
	}

INITIAL {
	debug = 0
	mrate = 0.0 : rate over last nwincount spikes 
	lastv = -100.0 : voltage at previous time step (for spike detector)
	i = 0 : default starting current
	lastt = 0.0
	nspikes = 0
	iset = 0
	nprc = 50
	prccount = 0
	prcskip = 4 : number of spikes to skip between PRC measures (allow stabiization of rate)
	prcskipcnt = 0
	prcintvl = (1000.0)/(nprc*rate)
	thisdelay = 0
	thisperturb = 0 : state variable: 0 - no perturbation done, 1, 
	tperturb = 0

}

BREAKPOINT {
	if(on == 1) { : adjust current only if enabled
		if(v >= thresh && lastv < thresh && nspikes < SSIZE) { : detect spikes
			iset = 0 : reset flag when spike detected
			nspikes = nspikes + 1 : count spike
			spike[nspikes] = t : save current spike time
			curr[nspikes] = i : save the current with this time too
			if(nspikes >= 2) {
			  if (debug) {
				VERBATIM
			   fprintf(stdout, "Spike %6.0f at %8.3f ms, isi = %8.3f\n", nspikes, lastt, spike[(int)nspikes]-spike[(int)(nspikes-1)]);
			   ENDVERBATIM
			  }
			  isi[nspikes] = spike[nspikes] - spike[nspikes-1]
			  cisi[nspikes] = (curr[nspikes] + curr[nspikes-1])/2
			  prcskipcnt = prcskipcnt + 1 : how many spikes have we skipped? 
			}
		  if(thisperturb == 2 && prccount <= nprc) {
				prcy[prccount] =  spike[nspikes] - spike[nspikes-1] : save the interval
				prcx[prccount] = thisdelay : and the delay time
				if (debug) {
				  VERBATIM
			    	fprintf(stdout, "Perturbation #: %d delay: %7.3f   isi: %7.3f\n", \
					(int) prccount, prcx[(int)prccount], prcy[(int)prccount]);
					ENDVERBATIM
				}
				prccount = prccount + 1
				thisdelay = prcintvl*prccount : compute the next interval
				prcskipcnt = 0 : reset the skip counter
				thisperturb = 0 : reset the perturbation flag
			}
	
		iset = 0 : clear the flag indicating we set current during this spike
		if(nspikes > nwincount && t > slopedur) { : wait until enough intervals measured to adjust I
: calculate the running mean firing rate over nwincount spikes
			if (debug) {
			  VERBATIM
			  fprintf(stdout, "time and spike count meet criteria\n");
			  ENDVERBATIM
			}
			intvl = 0
				FROM j = 0 TO nwincount-1 { : compute average interval over previous n spikes
					intvl = intvl + (spike[nspikes-j] - spike[nspikes-j-1])
				}
				mrate = (1000.0)*nwincount/intvl : convert to rate

				if(thisperturb == 0) { : if we are not doing PRC, allow current adjustment
					reglin()
					:ilog = (log10(rate) - b) / a
					:i = pow(10, ilog ) : predict current for target rate
					i = (rate-b ) / a
					inew = i
					if (debug) {
					  VERBATIM
				  	fprintf(stdout, "current spike rate: %9.3f, target rate: %9.1f, ilog=%f i = %8.3f\n  ", mrate, rate, ilog, i);
				  	fprintf(stdout, "reglin: slope %7.3f  intcpt %7.3f\n", a, b);
				  	ENDVERBATIM
				  }
				}
			} : end of current adjusment block for detected spike
		} : end of spike detection
: 
: ramp current - independent of spike detection during start of test
:
	  if(t < slopedur) { : start with a ramp - no dependence on spikes
		  i = slope * t
		  inew = i
		  nspikes = 0
	  } : end of ramp block
	
:
: if no spikes occur after ramp detection, current may be too low - try increasing it a bit
: 		
	  if((t >= slopedur) && (nspikes < SSIZE) && ((t-lastt) > (2000.0/rate))) { : rate too low - so increase I 
			  i = rdeltaplus * i : (rate - b)/a
			  prcskipcnt = 0 : reset the skip count... need to stabilize again
:				if (debug) {
:				  VERBATIM
:				fprintf(stdout, "Updating i at %7.3f, lastt=%7.3f, waiting = %7.3f, new i = %7.3f\n", \
:					t, lastt, 2.0*1000.0/rate, i);
:				ENDVERBATIM
:				}
			  iset = 1
:				  lastt = t : update time when we changed current last --  
	  } : end of current update block for long intervals

: code for prc - wait until after spike then kick it
: Must meet all of these conditions:
: time after most recent spike is current delay time for test pulse
: we are not already perturbing
: the recent rate is within limits of target rate
: we have not finished all the PRC intervals

	if (((t-spike[nspikes]) >= thisdelay)  && (thisperturb == 0) && 
	  ((mrate > (1.0-drate)*rate) && (mrate < (1.0+drate)*rate)) 
	  && (prccount < nprc) && (prcskipcnt >= prcskip)) { 

			VERBATIM
			fprintf(stdout, "Perturbing at t=%7.3f, last spike t = %7.3f, delay = %7.3f, skip count = %d\n", t, spike[(int) nspikes], thisdelay, (int) prcskipcnt);
			ENDVERBATIM
			thisperturb = 1 : set perturbation flag
			i = i + prcamp : step current
			tperturb = t
			prcskipcnt = 0 : reset the skip count flag
	} : end of block that initiates perturbation step 
	if(thisperturb == 1 && (t-tperturb) >= prcdur) {
		i = i - prcamp : step current back
		thisperturb = 2 : set flag - so we can save the interval when the next spike occurs
		VERBATIM
		fprintf(stdout, "Perturbation ended at t=%7.3f, tpert = %7.3f, prcdur = %7.3f\n", t, tperturb, prcdur);
		ENDVERBATIM
	} : end of block that terminates perturbation step
	}
	lastv = v : save voltage 	
  	lastt = t : save last spike time
}

PROCEDURE reglin()
{
:	compute the linear regression between two arrays x and y
:	regression is on y=ax+b
: 
: Note that for FI, data is closer to log-log in characters, so we transform before the fit
: and then return the values back in the original space.

	sx=0
	sy=0
	sx2=0
	sy2=0
	sxy=0
:VERBATIM
:fprintf(stdout, "nspikes: %d\n", (int) nspikes);
:ENDVERBATIM
		FROM j = 1 TO nspikes-1 {
			c = ((curr[j+1] + curr[j])/2) : would be log10(cur...)
			sx = sx + c
			dy = (1000.0/(spike[j+1] - spike[j]) ): convert isi (msec) to rate (s/s) : would be log10(1000/s1-s2)
:VERBATIM
:fprintf(stdout, "c: %7.3f   rate: %7.3f \n", c, dy);
:ENDVERBATIM
			sy = sy + dy
			sx2 = sx2 + c*c
			sy2 = sy2 + dy*dy
			sxy = sxy + c*dy
  	}

	a1 = (nspikes * sxy) - ( sx * sy )
	a2 = (nspikes * sx2) - ( sx * sx )
	a = (a1 / a2)

	b = ( (sy - a * sx) / nspikes)

}

