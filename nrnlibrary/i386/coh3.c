/* Created by Language version: 6.2.0 */
/* NOT VECTORIZED */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib.h"
#undef PI
 
#include "md1redef.h"
#include "section.h"
#include "md2redef.h"

#if METHOD3
extern int _method3;
#endif

#undef exp
#define exp hoc_Exp
extern double hoc_Exp();
 
#define _threadargscomma_ /**/
#define _threadargs_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define TAmp _p[0]
#define TDur _p[1]
#define dD _p[2]
#define dF _p[3]
#define F _p[4]
#define k0 _p[5]
#define kmax _p[6]
#define taud _p[7]
#define kd _p[8]
#define tauf _p[9]
#define kf _p[10]
#define taus _p[11]
#define ks _p[12]
#define glu _p[13]
#define rseed _p[14]
#define latency _p[15]
#define latstd _p[16]
#define debug _p[17]
#define Identifier _p[18]
#define nZones _p[19]
#define nRequests _p[20]
#define nReleases _p[21]
#define LatencyDist (_p + 22)
#define EventDist (_p + 10022)
#define ScopDist (_p + 11022)
#define Fn _p[12022]
#define Dn _p[12023]
#define TTotal _p[12024]
#define ev_index _p[12025]
#define sc_index _p[12026]
#define nVesicles (_p + 12027)
#define XMTR (_p + 13027)
#define TCount (_p + 14027)
#define ZoneActive (_p + 15027)
#define CaDn _p[16027]
#define CaFn _p[16028]
#define CaDi _p[16029]
#define CaFi _p[16030]
#define eta _p[16031]
#define tSpike _p[16032]
#define tstep _p[16033]
#define tRelease (_p + 16034)
#define relThisSpike (_p + 17034)
#define ZoneLatency (_p + 18034)
#define nTotal _p[19034]
#define hasReleased _p[19035]
#define tz _p[19036]
#define tspike _p[19037]
#define latzone _p[19038]
#define iZone _p[19039]
#define gindex _p[19040]
#define gammasize _p[19041]
#define scrand _p[19042]
#define DnVesicles (_p + 19043)
#define DXMTR (_p + 20043)
#define DTCount (_p + 21043)
#define DZoneActive (_p + 22043)
#define _g _p[23043]
#define _tsav _p[23044]
#define _nd_area  *_ppvar[0]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 /* declaration of user functions */
 static double _hoc_release();
 static double _hoc_update();
 static int _mechtype;
extern int nrn_get_mechtype();
 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(_ho) Object* _ho; { void* create_point_process();
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt();
 static double _hoc_loc_pnt(_vptr) void* _vptr; {double loc_point_process();
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(_vptr) void* _vptr; {double has_loc_point();
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(_vptr)void* _vptr; {
 double get_loc_point_process(); return (get_loc_point_process(_vptr));
}
 static _hoc_setdata(_vptr) void* _vptr; { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 /* connect user functions to hoc names */
 static IntFunc hoc_intfunc[] = {
 0,0
};
 static struct Member_func {
	char* _name; double (*_member)();} _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 "release", _hoc_release,
 "update", _hoc_update,
 0, 0
};
 /* declare global and static user variables */
#define n0 n0_COH3
 double n0 = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "n0_COH3", "1",
 "TAmp", "mM",
 "TDur", "ms",
 "dD", "1",
 "dF", "1",
 "F", "1",
 "k0", "/ms",
 "kmax", "/ms",
 "taud", "ms",
 "kd", "1",
 "tauf", "ms",
 "kf", "1",
 "taus", "ms",
 "ks", "1",
 "glu", "mM",
 "rseed", "1",
 "latency", "ms",
 "latstd", "ms",
 "nVesicles", "1",
 "XMTR", "mM",
 "TCount", "1",
 "ZoneActive", "1",
 "nZones", "1",
 "nRequests", "1",
 "nReleases", "1",
 "Fn", "1",
 "Dn", "1",
 0,0
};
 static double TCount0 = 0;
 static double XMTR0 = 0;
 static double ZoneActive0 = 0;
 static double delta_t = 1;
 static double nVesicles0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "n0_COH3", &n0_COH3,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(), nrn_init(), nrn_state();
 static void nrn_cur(), nrn_jacob();
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 
static int _ode_count();
 /* connect range variables in _p that hoc is supposed to know about */
 static char *_mechanism[] = {
 "6.2.0",
"COH3",
 "TAmp",
 "TDur",
 "dD",
 "dF",
 "F",
 "k0",
 "kmax",
 "taud",
 "kd",
 "tauf",
 "kf",
 "taus",
 "ks",
 "glu",
 "rseed",
 "latency",
 "latstd",
 "debug",
 "Identifier",
 0,
 "nZones",
 "nRequests",
 "nReleases",
 "LatencyDist[10000]",
 "EventDist[1000]",
 "ScopDist[1000]",
 "Fn",
 "Dn",
 "TTotal",
 "ev_index",
 "sc_index",
 0,
 "nVesicles[1000]",
 "XMTR[1000]",
 "TCount[1000]",
 "ZoneActive[1000]",
 0,
 0};
 
static void nrn_alloc(_prop)
	Prop *_prop;
{
	Prop *prop_ion, *need_memb();
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 23045, _prop);
 	/*initialize range parameters*/
 	TAmp = 1;
 	TDur = 0.5;
 	dD = 0.02;
 	dF = 0.02;
 	F = 0.5;
 	k0 = 0.0005714;
 	kmax = 0.04;
 	taud = 50;
 	kd = 0.7;
 	tauf = 10;
 	kf = 0.5;
 	taus = 1;
 	ks = 0.5;
 	glu = 1;
 	rseed = 0;
 	latency = 0;
 	latstd = 0;
 	debug = 0;
 	Identifier = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 23045;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 2, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static _initlists();
 static _net_receive();
 typedef (*_Pfrv)();
 extern _Pfrv* pnt_receive;
 extern short* pnt_receive_size;
 _coh3_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func,
	 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
  hoc_register_dparam_size(_mechtype, 2);
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 COH3 /Users/pbmanis/Desktop/Python/nrnlibrary/nrnlibrary/i386/coh3.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "Calyx of Held Version 3";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
static release();
static update();
 
static int  release (  )  {
   nTotal = 0.0 ;
   tz = 0.0 ;
   latzone = 0.0 ;
   gindex = 1.0 ;
   
/*VERBATIM*/
		if (debug == 1) {
// 			  fprintf(stderr, "** Entering release @ t=%f, nZones = %d\n", t, (int)nZones);
		}
 iZone = 0.0 ;
   {int  _li ;for ( _li = 0 ; _li <= ( ((int) nZones ) - 1 ) ; _li ++ ) {
     
/*VERBATIM*/
			if (debug == 1) {
/* 
				  fprintf(stderr, "      Checking Zones, zone= %.0f of %d, ", iZone, (int)nZones);
				  fprintf(stderr, " relthisspike = %.0f, trel: %.3f, tspike: %.3f\n", \
				  relThisSpike[(int)iZone], tRelease[(int)iZone], tSpike);
 */
			}
 scrand = scop_random ( ) ;
     if ( sc_index < 1000.0 ) {
       ScopDist [ ((int) sc_index ) ] = scrand ;
       sc_index = sc_index + 1.0 ;
       }
     
/*VERBATIM*/
			if (debug == 1) {
// 				fprintf(stderr, "        scop_random= %8.3f, fn*dn = %8.3f\n", scrand, Fn*Dn);
			}
 if ( relThisSpike [ ((int) iZone ) ]  == 0.0  && ( tRelease [ ((int) iZone ) ] < tSpike ) ) {
       
/*VERBATIM*/
			if (debug == 1) {
				fprintf(stderr, "Zone %.0f: Checking Probability\n", iZone);
			}
 if ( scrand < Fn * Dn ) {
         
/*VERBATIM*/
				if (debug == 1) {
					fprintf(stderr, "   *** Zone %d RELEASES at %7.2f.  scrnd: %8.3f   fndn: %8.3f\n", \
					(int)iZone, tSpike+latzone, scrand, Fn*Dn);
				}
 nReleases = nReleases + 1.0 ;
         tRelease [ ((int) iZone ) ] = tSpike ;
         TTotal = TTotal + 1.0 ;
         hasReleased = 0.0 ;
         relThisSpike [ ((int) iZone ) ] = 1.0 ;
         ZoneActive [ ((int) iZone ) ] = 1.0 ;
         
/*VERBATIM*/
					gindex = (int)(10000.0*scop_random()); // select a value randomly
//					fprintf(stderr, "       gindex: %d", (int)gindex);
 if ( gindex <= 0.0 ) {
           gindex = 1.0 ;
           }
         if ( gindex > 10000.0 ) {
           gindex = 10000.0 ;
           }
         latzone = LatencyDist [ ((int) gindex ) ] ;
         
/*VERBATIM*/
					if (debug == 1) {
						fprintf(stderr, "     gindex: %d zone: %d  latency: %g \n", \
						(int)gindex, (int)iZone, latzone);
					}
 if ( ev_index < 1000.0 ) {
           EventDist [ ((int) ev_index ) ] = latzone ;
           ev_index = ev_index + 1.0 ;
           }
         TCount [ ((int) iZone ) ] = ( latzone / dt ) + ( 5.0 * TDur / dt ) ;
         ZoneLatency [ ((int) iZone ) ] = latzone ;
         XMTR [ ((int) iZone ) ] = XMTR [ ((int) iZone ) ] + TAmp ;
         }
       else {
         relThisSpike [ ((int) iZone ) ] = 2.0 ;
         
/*VERBATIM*/
					if (debug == 1) {
						fprintf(stderr, "  ### Zone %d: Looked but did not release - scrnd: %8.3f   fndn: %8.3f\n",\
						 (int) iZone, scrand, Fn*Dn);
					}
 }
       }
     if ( ZoneActive [ ((int) iZone ) ] > 0.0 ) {
       if ( t >= tRelease [ ((int) iZone ) ] ) {
         TCount [ ((int) iZone ) ] = TCount [ ((int) iZone ) ] - 1.0 ;
         if ( TCount [ ((int) iZone ) ] < 0.0 ) {
           XMTR [ ((int) iZone ) ] = 0.0 ;
           ZoneActive [ ((int) iZone ) ] = 0.0 ;
           TTotal = 0.0 ;
           TCount [ ((int) iZone ) ] = 0.0 ;
           
/*VERBATIM*/
						if (debug == 1) {
						//	fprintf(stderr, " Reset XMTR\n");
							}
 }
         else {
           tz = t - ( tRelease [ ((int) iZone ) ] + ZoneLatency [ ((int) iZone ) ] ) ;
           if ( tz < 0.0 ) {
             XMTR [ ((int) iZone ) ] = 0.0 ;
             }
           else {
             
/*VERBATIM*/
						if (hasReleased == 0 && debug == 1) {
						//	fprintf(stderr,  "First latency is %8.4f\n", t);
						}
 hasReleased = 1.0 ;
             XMTR [ ((int) iZone ) ] = TAmp * ( 1.0 - exp ( - tz / ( TDur / 3.0 ) ) ) * exp ( - ( tz - ( TDur / 3.0 ) ) / TDur ) ;
             }
           
/*VERBATIM*/
						if (debug == 1) {
						//	fprintf(stderr, " Updated XMTR\n");
						}
 }
         }
       }
     iZone = iZone + 1.0 ;
     } }
    return 0; }
 
static double _hoc_release(_vptr) void* _vptr; {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 release (  ) ;
 return(_r);
}
 
static int  update (  _ltstep )  
	double _ltstep ;
 {
   
/*VERBATIM*/
		if (debug == 2 ){
		fprintf(stdout, "update start t = %f ts=%f: F=%7.2f CaDi = %g CaFi = %g\n", \
		t, tstep, F, CaDi, CaFi);
		fprintf(stdout, "	vars:  taud=%g: tauf=%g kd = %g kmax= %g\n", taud, tauf, kd, kmax);
		}
 if ( _ltstep < 0.0 ) {
     return ( _threadargscomma_ 0.0 ) ;
     }
   CaDi = CaDi + dD ;
   CaFi = CaFi + dF ;
   
/*VERBATIM*/
		if (debug >= 2 ){
			fprintf(stdout, "	 CaDi = %g CaFi = %g\n", CaDi, CaFi);
		}
 CaDn = CaDi * exp ( - _ltstep / taud ) ;
   CaFn = CaFi * exp ( - _ltstep / tauf ) ;
   
/*VERBATIM*/
		if (debug >= 2 ){
			fprintf(stdout, "	 CaDn = %g CaFn = %g\n", CaDn, CaFn);
		}
 eta = ( kd / CaDi + 1.0 ) / ( kd / CaDi + exp ( - _ltstep / taud ) ) ;
   eta = pow( eta , ( - ( kmax - k0 ) * taud ) ) ;
   
/*VERBATIM*/
		if (debug >= 2 ){
			fprintf(stdout, "eta: %g\n", eta);
		}
 Dn = 1.0 - ( 1.0 - ( 1.0 - Fn ) * Dn ) * exp ( - k0 * _ltstep ) * eta ;
   Fn = F + ( 1.0 - F ) / ( 1.0 + kf / CaFn ) ;
   CaDi = CaDn ;
   CaFi = CaFn ;
   
/*VERBATIM*/
		if (debug >= 2 ){
			fprintf(stdout, "	  Fn=%7.2f Dn: %7.2f  CaDi = %g CaFi = %g,\n", \
			Fn, Dn, CaDi, CaFi);
		}
  return 0; }
 
static double _hoc_update(_vptr) void* _vptr; {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 update (  *getarg(1) ) ;
 return(_r);
}
 
static _net_receive (_pnt, _args, _lflag) Point_process* _pnt; double* _args; double _lflag; 
{    _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t; {
   update ( _threadargscomma_ t - tSpike ) ;
   tSpike = t ;
   TTotal = 0.0 ;
   {int  _li ;for ( _li = 0 ; _li <= ((int) nZones ) - 1 ; _li ++ ) {
     relThisSpike [ _li ] = 0.0 ;
     } }
   
/*VERBATIM*/
	  if (debug == 1) {
	    fprintf(stderr, "  ---> Spike at t = %9.3f\n", tSpike);
	  }
 nRequests = nRequests + 1.0 ;
   } }
 
static int _ode_count(_type)int _type; { hoc_execerror("COH3", "cannot be used with CVODE");}

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
 for (_i=0; _i<1000; _i++) TCount[_i] = TCount0;
 for (_i=0; _i<1000; _i++) XMTR[_i] = XMTR0;
 for (_i=0; _i<1000; _i++) ZoneActive[_i] = ZoneActive0;
 for (_i=0; _i<1000; _i++) nVesicles[_i] = nVesicles0;
 {
   
/*VERBATIM*/
		if(debug == 1) {
			fprintf(stdout, "Initial for COH3, instance %d\n", (int)Identifier);
			fprintf(stdout, "Latency: %f   Var: %f\n", latency, latstd);
		}
	//	fprintf(stdout, "Debug flag: %d for instance %d\n", (int)debug, (int)Identifier);
 nTotal = 0.0 ;
   TTotal = 0.0 ;
   nRequests = 0.0 ;
   nReleases = 0.0 ;
   set_seed ( rseed ) ;
   nVesicles [ 0 ] = n0 ;
   XMTR [ 0 ] = 0.0 ;
   tRelease [ 0 ] = 0.0 ;
   hasReleased = 0.0 ;
   relThisSpike [ 0 ] = 0.0 ;
   ZoneActive [ 0 ] = 0.0 ;
   tSpike = - 1000.0 ;
   latzone = 0.0 ;
   gindex = 0.0 ;
   ev_index = 0.0 ;
   sc_index = 0.0 ;
   scrand = 0.0 ;
   CaDi = 1.0 ;
   CaFi = 0.0 ;
   CaDn = 1.0 ;
   CaFn = 0.0 ;
   Fn = F ;
   Dn = 1.0 ;
   iZone = 0.0 ;
   {int  _li ;for ( _li = 0 ; _li <= ( ((int) nZones ) - 1 ) ; _li ++ ) {
     
/*VERBATIM*/
			if (debug == 1) {
		//		fprintf(stderr, " Initializing release for zone= %.0f of %.0f\n", iZone, nZones);
			}
 nVesicles [ _li ] = n0 ;
     XMTR [ _li ] = 0.0 ;
     tRelease [ _li ] = 0.0 ;
     relThisSpike [ _li ] = 0.0 ;
     ZoneActive [ _li ] = 0.0 ;
     iZone = iZone + 1.0 ;
     } }
   update ( _threadargscomma_ t - tSpike ) ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _tsav = -1e20;
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{
} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
 double _break, _save;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _break = t + .5*dt; _save = t;
 v=_v;
{
 { {
 for (; t < _break; t += dt) {
 error =  release();
 if(error){fprintf(stderr,"at line 180 in file coh3.mod:\n	SOLVE release\n"); nrn_complain(_p); abort_run(error);}
 
}}
 t = _save;
 } {
   }
}}

}

static terminal(){}

static _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}
