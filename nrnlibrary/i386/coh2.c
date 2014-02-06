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
#define rseed _p[0]
#define PR0 _p[1]
#define nzones _p[2]
#define ntot _p[3]
#define Ttot _p[4]
#define Camp _p[5]
#define NVes (_p + 6)
#define RScale _p[1006]
#define RO (_p + 1007)
#define CGLU (_p + 1009)
#define C _p[2009]
#define Cn _p[2010]
#define pv (_p + 2011)
#define Ccount _p[3011]
#define Cncount _p[3012]
#define Tcnt (_p + 3013)
#define tspike _p[4013]
#define trel (_p + 4014)
#define relthisspike (_p + 5014)
#define km _p[6014]
#define inf (_p + 6015)
#define tau (_p + 6017)
#define fac (_p + 6019)
#define up _p[6021]
#define ui _p[6022]
#define DNVes (_p + 6023)
#define DRScale _p[7023]
#define DRO (_p + 7024)
#define DCGLU (_p + 7026)
#define DC _p[8026]
#define DCn _p[8027]
#define Dpv (_p + 8028)
#define DCcount _p[9028]
#define DCncount _p[9029]
#define DTcnt (_p + 9030)
#define _g _p[10030]
#define _tsav _p[10031]
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
 static double _hoc_onerel();
 static double _hoc_prel();
 static double _hoc_release();
 static double _hoc_rates();
 static double _hoc_unirand();
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
 "onerel", _hoc_onerel,
 "prel", _hoc_prel,
 "release", _hoc_release,
 "rates", _hoc_rates,
 "unirand", _hoc_unirand,
 0, 0
};
#define unirand unirand_COH2
 extern double unirand();
 /* declare global and static user variables */
#define Cndur Cndur_COH2
 double Cndur = 0.5;
#define Cnres Cnres_COH2
 double Cnres = 0;
#define Cnamp Cnamp_COH2
 double Cnamp = 0.01;
#define Cdur Cdur_COH2
 double Cdur = 0.25;
#define Cres Cres_COH2
 double Cres = 0.05;
#define KO KO_COH2
 double KO[2];
#define KC KC_COH2
 double KC[2];
#define R R_COH2
 double R = 0;
#define Tdur Tdur_COH2
 double Tdur = 0.1;
#define Tamp Tamp_COH2
 double Tamp = 1;
#define kd kd_COH2
 double kd = 0.002;
#define ke ke_COH2
 double ke = 80;
#define n0 n0_COH2
 double n0 = 1;
#define usetable usetable_COH2
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_COH2", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "n0_COH2", "1",
 "ke_COH2", "/mM",
 "kd_COH2", "/ms",
 "Cres_COH2", "mM",
 "Cdur_COH2", "ms",
 "Cnamp_COH2", "mM",
 "Cnres_COH2", "mM",
 "Cndur_COH2", "ms",
 "Tamp_COH2", "mM",
 "Tdur_COH2", "ms",
 "R_COH2", "1",
 "KO_COH2", "/mM /ms",
 "KC_COH2", "/ms",
 "NVes", "1",
 "RScale", "1",
 "RO", "1",
 "CGLU", "mM",
 "C", "mM",
 "Cn", "mM",
 "pv", "1",
 "Ccount", "1",
 "Cncount", "1",
 "Tcnt", "1",
 "nzones", "1",
 "ntot", "1",
 "Ttot", "mM",
 "Camp", "mM",
 0,0
};
 static double Cncount0 = 0;
 static double Ccount0 = 0;
 static double Cn0 = 0;
 static double C0 = 0;
 static double CGLU0 = 0;
 static double NVes0 = 0;
 static double RO0 = 0;
 static double RScale0 = 0;
 static double Tcnt0 = 0;
 static double delta_t = 1;
 static double pv0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "n0_COH2", &n0_COH2,
 "ke_COH2", &ke_COH2,
 "kd_COH2", &kd_COH2,
 "Cres_COH2", &Cres_COH2,
 "Cdur_COH2", &Cdur_COH2,
 "Cnamp_COH2", &Cnamp_COH2,
 "Cnres_COH2", &Cnres_COH2,
 "Cndur_COH2", &Cndur_COH2,
 "Tamp_COH2", &Tamp_COH2,
 "Tdur_COH2", &Tdur_COH2,
 "R_COH2", &R_COH2,
 "usetable_COH2", &usetable_COH2,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 "KO_COH2", KO_COH2, 2,
 "KC_COH2", KC_COH2, 2,
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
"COH2",
 0,
 "rseed",
 "PR0",
 "nzones",
 "ntot",
 "Ttot",
 "Camp",
 0,
 "NVes[1000]",
 "RScale",
 "RO[2]",
 "CGLU[1000]",
 "C",
 "Cn",
 "pv[1000]",
 "Ccount",
 "Cncount",
 "Tcnt[1000]",
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
 	_p = nrn_prop_data_alloc(_mechtype, 10032, _prop);
 	/*initialize range parameters*/
  }
 	_prop->param = _p;
 	_prop->param_size = 10032;
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
 _coh2_reg() {
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
 	ivoc_help("help ?1 COH2 /Users/pbmanis/Desktop/Python/nrnlibrary/nrnlibrary/i386/coh2.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_inf[2];
 static double *_t_fac[2];
 static double *_t_tau[2];
static int _reset;
static char *modelname = "Calyx of Held Version 2";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
static _f_rates();
static onerel();
static prel();
static release();
static rates();
 static _n_rates();
 
static int  release (  )  {
   prel ( _threadargs_ ) ;
   ntot = 0.0 ;
   {int  _li ;for ( _li = 0 ; _li <= ((int) nzones ) - 1 ; _li ++ ) {
     if ( unirand ( _threadargs_ ) < dt * km ) {
       NVes [ _li ] = NVes [ _li ] + 1.0 ;
       }
     if ( NVes [ _li ] > 0.0  && unirand ( _threadargs_ ) < dt * kd ) {
       NVes [ _li ] = NVes [ _li ] - 1.0 ;
       }
     if ( Cn > 0.0 ) {
       if ( unirand ( _threadargs_ ) < dt * ke * Cn ) {
         NVes [ _li ] = NVes [ _li ] + 1.0 ;
         }
       }
     onerel ( _threadargscomma_ ((double) _li ) ) ;
     if ( CGLU [ _li ] > 0.0 ) {
       Tcnt [ _li ] = Tcnt [ _li ] - 1.0 ;
       if ( Tcnt [ _li ] < 0.0 ) {
         CGLU [ _li ] = 0.0 ;
         Ttot = 0.0 ;
         }
       }
     ntot = ntot + NVes [ _li ] ;
     } }
   if ( C > Cres ) {
     Ccount = Ccount - 1.0 ;
     if ( Ccount < 0.0 ) {
       C = Cres ;
       Cn = Cnamp ;
       Cncount = Cndur / dt ;
       }
     }
   if ( Cn > Cnres ) {
     Cncount = Cncount - 1.0 ;
     if ( Cncount < 0.0 ) {
       Cn = Cnres ;
       }
     }
   
/*VERBATIM*/
	return 0;
  return 0; }
 
static double _hoc_release(_vptr) void* _vptr; {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 release (  ) ;
 return(_r);
}
 
static int  onerel (  _lizone )  
	double _lizone ;
 {
   if ( relthisspike [ ((int) _lizone ) ]  == 0.0 ) {
     if ( trel [ ((int) _lizone ) ] < tspike ) {
       ui = R * pv [ ((int) _lizone ) ] / RScale ;
       if ( NVes [ ((int) _lizone ) ] > 0.0  && unirand ( _threadargs_ ) < ui ) {
         NVes [ ((int) _lizone ) ] = NVes [ ((int) _lizone ) ] - 1.0 ;
         CGLU [ ((int) _lizone ) ] = Tamp ;
         Tcnt [ ((int) _lizone ) ] = Tdur / dt ;
         trel [ ((int) _lizone ) ] = t ;
         Ttot = Ttot + 1.0 ;
         relthisspike [ ((int) _lizone ) ] = 1.0 ;
         }
       else {
         relthisspike [ ((int) _lizone ) ] = 2.0 ;
         }
       }
     }
    return 0; }
 
static double _hoc_onerel(_vptr) void* _vptr; {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 onerel (  *getarg(1) ) ;
 return(_r);
}
 
static int  prel (  )  {
   rates ( _threadargscomma_ C ) ;
   R = 1.0 ;
   {int  _li ;for ( _li = 0 ; _li <= 1 ; _li ++ ) {
     RO [ _li ] = RO [ _li ] + fac [ _li ] * ( inf [ _li ] - RO [ _li ] ) ;
     R = R * RO [ _li ] ;
     } }
    return 0; }
 
static double _hoc_prel(_vptr) void* _vptr; {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 prel (  ) ;
 return(_r);
}
 static double _mfac_rates, _tmin_rates;
 static _check_rates();
 static _check_rates() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_dt;
  if (!usetable) {return;}
  if (_sav_dt != dt) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_rates =  0.0 ;
   _tmax =  0.2 ;
   _dx = (_tmax - _tmin_rates)/200.; _mfac_rates = 1./_dx;
   for (_i=0, _x=_tmin_rates; _i < 201; _x += _dx, _i++) {
    _f_rates(_x);
    for (_j = 0; _j < 2; _j++) { _t_inf[_j][_i] = inf[_j];
}    for (_j = 0; _j < 2; _j++) { _t_fac[_j][_i] = fac[_j];
}    for (_j = 0; _j < 2; _j++) { _t_tau[_j][_i] = tau[_j];
}   }
   _sav_dt = dt;
  }
 }

 static rates(double _lC){ _check_rates();
 _n_rates(_lC);
 return;
 }

 static _n_rates(double _lC){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_rates(_lC); return; 
}
 _xi = _mfac_rates * (_lC - _tmin_rates);
 _i = (int) _xi;
 if (_xi <= 0.) {
 for (_j = 0; _j < 2; _j++) { inf[_j] = _t_inf[_j][0];
} for (_j = 0; _j < 2; _j++) { fac[_j] = _t_fac[_j][0];
} for (_j = 0; _j < 2; _j++) { tau[_j] = _t_tau[_j][0];
} return; }
 if (_i >= 200) {
 for (_j = 0; _j < 2; _j++) { inf[_j] = _t_inf[_j][200];
} for (_j = 0; _j < 2; _j++) { fac[_j] = _t_fac[_j][200];
} for (_j = 0; _j < 2; _j++) { tau[_j] = _t_tau[_j][200];
} return; }
 _theta = _xi - (double)_i;
 for (_j = 0; _j < 2; _j++) {double *_t = _t_inf[_j]; inf[_j] = _t[_i] + _theta*(_t[_i+1] - _t[_i]);}
 for (_j = 0; _j < 2; _j++) {double *_t = _t_fac[_j]; fac[_j] = _t[_i] + _theta*(_t[_i+1] - _t[_i]);}
 for (_j = 0; _j < 2; _j++) {double *_t = _t_tau[_j]; tau[_j] = _t[_i] + _theta*(_t[_i+1] - _t[_i]);}
 }

 
static int  _f_rates (  _lC )  
	double _lC ;
 {
   double _la , _lb ;
 {int  _lj ;for ( _lj = 0 ; _lj <= 1 ; _lj ++ ) {
     _la = KO [ _lj ] * _lC ;
     _lb = KC [ _lj ] ;
     tau [ _lj ] = 1.0 / ( _la + _lb ) ;
     inf [ _lj ] = _la / ( _la + _lb ) ;
     fac [ _lj ] = ( 1.0 - exp ( - dt / tau [ _lj ] ) ) ;
     } }
    return 0; }
 
static double _hoc_rates(_vptr) void* _vptr; {
 double _r;
    _hoc_setdata(_vptr);
  _r = 1.;
 rates (  *getarg(1) ) ;
 return(_r);
}
 
double unirand (  )  {
   double _lunirand;
 return ( _threadargscomma_ scop_random ( ) ) ;
   
return _lunirand;
 }
 
static double _hoc_unirand(_vptr) void* _vptr; {
 double _r;
    _hoc_setdata(_vptr);
 _r =  unirand (  ) ;
 return(_r);
}
 
static _net_receive (_pnt, _args, _lflag) Point_process* _pnt; double* _args; double _lflag; 
{    _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t; {
   if ( C < Camp ) {
     C = Camp ;
     tspike = t ;
     Ccount = Cdur / dt ;
     Ttot = 0.0 ;
     {int  _li ;for ( _li = 0 ; _li <= ((int) nzones ) - 1 ; _li ++ ) {
       relthisspike [ _li ] = 0.0 ;
       } }
     }
   } }
 
static int _ode_count(_type)int _type; { hoc_execerror("COH2", "cannot be used with CVODE");}

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  Cncount = Cncount0;
  Ccount = Ccount0;
  Cn = Cn0;
  C = C0;
 for (_i=0; _i<1000; _i++) CGLU[_i] = CGLU0;
 for (_i=0; _i<1000; _i++) NVes[_i] = NVes0;
 for (_i=0; _i<2; _i++) RO[_i] = RO0;
  RScale = RScale0;
 for (_i=0; _i<1000; _i++) Tcnt[_i] = Tcnt0;
 for (_i=0; _i<1000; _i++) pv[_i] = pv0;
 {
   NVes [ 0 ] = n0 ;
   CGLU [ 0 ] = 0.0 ;
   trel [ 0 ] = 0.0 ;
   relthisspike [ 0 ] = 0.0 ;
   tspike = 0.0 ;
   pv [ 0 ] = PR0 ;
   {int  _li ;for ( _li = 0 ; _li <= ((int) nzones ) - 1 ; _li ++ ) {
     NVes [ _li ] = n0 ;
     CGLU [ _li ] = 0.0 ;
     trel [ _li ] = 0.0 ;
     relthisspike [ _li ] = 0.0 ;
     pv [ _li ] = PR0 ;
     } }
   R = 0.0 ;
   C = Cres ;
   Cn = Cnres ;
   ntot = 0.0 ;
   Ttot = 0.0 ;
   set_seed ( rseed ) ;
   km = n0 * kd ;
   KO [ 0 ] = 150.0 ;
   KO [ 1 ] = 1.0 ;
   KC [ 0 ] = 30.0 ;
   KC [ 1 ] = 0.1 ;
   C = Camp ;
   rates ( _threadargscomma_ C ) ;
   {int  _li ;for ( _li = 0 ; _li <= 1 ; _li ++ ) {
     RO [ _li ] = inf [ _li ] ;
     } }
   prel ( _threadargs_ ) ;
   RScale = sqrt ( R ) ;
   C = Cres ;
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
 if(error){fprintf(stderr,"at line 144 in file coh2.mod:\n	SOLVE release\n"); nrn_complain(_p); abort_run(error);}
 
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
  for (_i=0; _i < 2; _i++) {  _t_inf[_i] = makevector(201*sizeof(double)); }
  for (_i=0; _i < 2; _i++) {  _t_fac[_i] = makevector(201*sizeof(double)); }
  for (_i=0; _i < 2; _i++) {  _t_tau[_i] = makevector(201*sizeof(double)); }
_first = 0;
}
