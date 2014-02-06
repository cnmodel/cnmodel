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
#define on _p[0]
#define rate _p[1]
#define drate _p[2]
#define thresh _p[3]
#define slope _p[4]
#define slopedur _p[5]
#define rdeltaplus _p[6]
#define nwincount _p[7]
#define prcdur _p[8]
#define prcamp _p[9]
#define spike (_p + 10)
#define curr (_p + 10011)
#define isi (_p + 20012)
#define cisi (_p + 30013)
#define prcx (_p + 40014)
#define prcy (_p + 40064)
#define i _p[40114]
#define nspikes _p[40115]
#define prccount _p[40116]
#define intvl _p[40117]
#define mrate _p[40118]
#define lastt _p[40119]
#define lastv _p[40120]
#define sx _p[40121]
#define sy _p[40122]
#define sxy _p[40123]
#define sx2 _p[40124]
#define sy2 _p[40125]
#define dy _p[40126]
#define a _p[40127]
#define b _p[40128]
#define c _p[40129]
#define j _p[40130]
#define a1 _p[40131]
#define a2 _p[40132]
#define iset _p[40133]
#define inew _p[40134]
#define ilog _p[40135]
#define nprc _p[40136]
#define prcskip _p[40137]
#define prcskipcnt _p[40138]
#define prcintvl _p[40139]
#define thisdelay _p[40140]
#define thisperturb _p[40141]
#define tperturb _p[40142]
#define debug _p[40143]
#define _g _p[40144]
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
 static double _hoc_reglin();
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
 "reglin", _hoc_reglin,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "on", "1",
 "rate", "1/s",
 "thresh", "mV",
 "slope", "nA/ms",
 "slopedur", "ms",
 "nwincount", "1",
 "spike", "ms",
 "curr", "nA",
 "isi", "ms",
 "cisi", "nA",
 "prcx", "ms",
 "prcy", "ms",
 "i", "nA",
 "nspikes", "1",
 "prccount", "1",
 0,0
};
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
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
 /* connect range variables in _p that hoc is supposed to know about */
 static char *_mechanism[] = {
 "6.2.0",
"spikerate",
 "on",
 "rate",
 "drate",
 "thresh",
 "slope",
 "slopedur",
 "rdeltaplus",
 "nwincount",
 "prcdur",
 "prcamp",
 0,
 "spike[10001]",
 "curr[10001]",
 "isi[10001]",
 "cisi[10001]",
 "prcx[50]",
 "prcy[50]",
 "i",
 "nspikes",
 "prccount",
 0,
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
 	_p = nrn_prop_data_alloc(_mechtype, 40145, _prop);
 	/*initialize range parameters*/
 	on = 1;
 	rate = 20;
 	drate = 0.2;
 	thresh = 0;
 	slope = 0.0005;
 	slopedur = 1000;
 	rdeltaplus = 2;
 	nwincount = 5;
 	prcdur = 0.1;
 	prcamp = 0.1;
  }
 	_prop->param = _p;
 	_prop->param_size = 40145;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 2, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static _initlists();
 _spikerate_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func,
	 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
  hoc_register_dparam_size(_mechtype, 2);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 spikerate /Users/pbmanis/Desktop/Python/nrnlibrary/nrnlibrary/i386/spikerate.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
static reglin();
 
static int  reglin (  )  {
   sx = 0.0 ;
   sy = 0.0 ;
   sx2 = 0.0 ;
   sy2 = 0.0 ;
   sxy = 0.0 ;
   {int  _lj ;for ( _lj = 1 ; _lj <= ((int) nspikes ) - 1 ; _lj ++ ) {
     c = ( ( curr [ _lj + 1 ] + curr [ _lj ] ) / 2.0 ) ;
     sx = sx + c ;
     dy = ( 1000.0 / ( spike [ _lj + 1 ] - spike [ _lj ] ) ) ;
     sy = sy + dy ;
     sx2 = sx2 + c * c ;
     sy2 = sy2 + dy * dy ;
     sxy = sxy + c * dy ;
     } }
   a1 = ( nspikes * sxy ) - ( sx * sy ) ;
   a2 = ( nspikes * sx2 ) - ( sx * sx ) ;
   a = ( a1 / a2 ) ;
   b = ( ( sy - a * sx ) / nspikes ) ;
    return 0; }
 
static double _hoc_reglin(_vptr) void* _vptr; {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 reglin (  ) ;
 return(_r);
}

static void initmodel() {
  int _i; double _save;_ninits++;
{
 {
   debug = 0.0 ;
   mrate = 0.0 ;
   lastv = - 100.0 ;
   i = 0.0 ;
   lastt = 0.0 ;
   nspikes = 0.0 ;
   iset = 0.0 ;
   nprc = 50.0 ;
   prccount = 0.0 ;
   prcskip = 4.0 ;
   prcskipcnt = 0.0 ;
   prcintvl = ( 1000.0 ) / ( nprc * rate ) ;
   thisdelay = 0.0 ;
   thisperturb = 0.0 ;
   tperturb = 0.0 ;
   }

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
#if EXTRACELLULAR
 _nd = _ml->_nodelist[_iml];
 if (_nd->_extnode) {
    _v = NODEV(_nd) +_nd->_extnode->_v[0];
 }else
#endif
 {
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 }
 v = _v;
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   if ( on  == 1.0 ) {
     if ( v >= thresh  && lastv < thresh  && nspikes < 10001.0 ) {
       iset = 0.0 ;
       nspikes = nspikes + 1.0 ;
       spike [ ((int) nspikes ) ] = t ;
       curr [ ((int) nspikes ) ] = i ;
       if ( nspikes >= 2.0 ) {
         if ( debug ) {
           
/*VERBATIM*/
			   fprintf(stdout, "Spike %6.0f at %8.3f ms, isi = %8.3f\n", nspikes, lastt, spike[(int)nspikes]-spike[(int)(nspikes-1)]);
 }
         isi [ ((int) nspikes ) ] = spike [ ((int) nspikes ) ] - spike [ ((int) nspikes ) - 1 ] ;
         cisi [ ((int) nspikes ) ] = ( curr [ ((int) nspikes ) ] + curr [ ((int) nspikes ) - 1 ] ) / 2.0 ;
         prcskipcnt = prcskipcnt + 1.0 ;
         }
       if ( thisperturb  == 2.0  && prccount <= nprc ) {
         prcy [ ((int) prccount ) ] = spike [ ((int) nspikes ) ] - spike [ ((int) nspikes ) - 1 ] ;
         prcx [ ((int) prccount ) ] = thisdelay ;
         if ( debug ) {
           
/*VERBATIM*/
			    	fprintf(stdout, "Perturbation #: %d delay: %7.3f   isi: %7.3f\n", \
					(int) prccount, prcx[(int)prccount], prcy[(int)prccount]);
 }
         prccount = prccount + 1.0 ;
         thisdelay = prcintvl * prccount ;
         prcskipcnt = 0.0 ;
         thisperturb = 0.0 ;
         }
       iset = 0.0 ;
       if ( nspikes > nwincount  && t > slopedur ) {
         if ( debug ) {
           
/*VERBATIM*/
			  fprintf(stdout, "time and spike count meet criteria\n");
 }
         intvl = 0.0 ;
         {int  _lj ;for ( _lj = 0 ; _lj <= ((int) nwincount ) - 1 ; _lj ++ ) {
           intvl = intvl + ( spike [ ((int) nspikes ) - _lj ] - spike [ ((int) nspikes ) - _lj - 1 ] ) ;
           } }
         mrate = ( 1000.0 ) * nwincount / intvl ;
         if ( thisperturb  == 0.0 ) {
           reglin ( _threadargs_ ) ;
           i = ( rate - b ) / a ;
           inew = i ;
           if ( debug ) {
             
/*VERBATIM*/
				  	fprintf(stdout, "current spike rate: %9.3f, target rate: %9.1f, ilog=%f i = %8.3f\n  ", mrate, rate, ilog, i);
				  	fprintf(stdout, "reglin: slope %7.3f  intcpt %7.3f\n", a, b);
 }
           }
         }
       }
     if ( t < slopedur ) {
       i = slope * t ;
       inew = i ;
       nspikes = 0.0 ;
       }
     if ( ( t >= slopedur )  && ( nspikes < 10001.0 )  && ( ( t - lastt ) > ( 2000.0 / rate ) ) ) {
       i = rdeltaplus * i ;
       prcskipcnt = 0.0 ;
       iset = 1.0 ;
       }
     if ( ( ( t - spike [ ((int) nspikes ) ] ) >= thisdelay )  && ( thisperturb  == 0.0 )  && ( ( mrate > ( 1.0 - drate ) * rate )  && ( mrate < ( 1.0 + drate ) * rate ) )  && ( prccount < nprc )  && ( prcskipcnt >= prcskip ) ) {
       
/*VERBATIM*/
			fprintf(stdout, "Perturbing at t=%7.3f, last spike t = %7.3f, delay = %7.3f, skip count = %d\n", t, spike[(int) nspikes], thisdelay, (int) prcskipcnt);
 thisperturb = 1.0 ;
       i = i + prcamp ;
       tperturb = t ;
       prcskipcnt = 0.0 ;
       }
     if ( thisperturb  == 1.0  && ( t - tperturb ) >= prcdur ) {
       i = i - prcamp ;
       thisperturb = 2.0 ;
       
/*VERBATIM*/
		fprintf(stdout, "Perturbation ended at t=%7.3f, tpert = %7.3f, prcdur = %7.3f\n", t, tperturb, prcdur);
 }
     }
   lastv = v ;
   lastt = t ;
   }
 _current += i;

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
#if EXTRACELLULAR
 _nd = _ml->_nodelist[_iml];
 if (_nd->_extnode) {
    _v = NODEV(_nd) +_nd->_extnode->_v[0];
 }else
#endif
 {
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 }
 _g = _nrn_current(_v + .001);
 	{ _rhs = _nrn_current(_v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) += _rhs;
  }else
#endif
  {
	NODERHS(_nd) += _rhs;
  }
#if EXTRACELLULAR
 if (_nd->_extnode) {
   *_nd->_extnode->_rhs[0] += _rhs;
 }
#endif
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) -= _g;
  }else
#endif
  {
	NODED(_nd) -= _g;
  }
#if EXTRACELLULAR
 if (_nd->_extnode) {
   *_nd->_extnode->_d[0] += _g;
 }
#endif
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){

}

static terminal(){}

static _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}
