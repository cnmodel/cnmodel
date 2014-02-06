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
#define spike (_p + 0)
#define rseed _p[10001]
#define PR0 _p[10002]
#define nzones _p[10003]
#define ntot _p[10004]
#define Ttot _p[10005]
#define Camp _p[10006]
#define n (_p + 10007)
#define Rmax _p[11007]
#define RO (_p + 11008)
#define CGLU (_p + 11010)
#define C _p[12010]
#define Cn _p[12011]
#define pv (_p + 12012)
#define Ccount _p[13012]
#define Cncount _p[13013]
#define Tcnt (_p + 13014)
#define index _p[14014]
#define tspike _p[14015]
#define trel (_p + 14016)
#define relthisspike (_p + 15016)
#define km _p[16016]
#define inf (_p + 16017)
#define tau (_p + 16019)
#define fac (_p + 16021)
#define up _p[16023]
#define ui _p[16024]
#define Dn (_p + 16025)
#define DRmax _p[17025]
#define DRO (_p + 17026)
#define DCGLU (_p + 17028)
#define DC _p[18028]
#define DCn _p[18029]
#define Dpv (_p + 18030)
#define DCcount _p[19030]
#define DCncount _p[19031]
#define DTcnt (_p + 19032)
#define _g _p[20032]
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
 static double _hoc_rates();
 static double _hoc_release();
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
 "rates", _hoc_rates,
 "release", _hoc_release,
 "unirand", _hoc_unirand,
 0, 0
};
#define unirand unirand_COH
 extern double unirand();
 /* declare global and static user variables */
#define Cndur Cndur_COH
 double Cndur = 0.5;
#define Cnres Cnres_COH
 double Cnres = 0;
#define Cnamp Cnamp_COH
 double Cnamp = 0.01;
#define Cdur Cdur_COH
 double Cdur = 0.25;
#define Cres Cres_COH
 double Cres = 0;
#define KO KO_COH
 double KO[2];
#define KC KC_COH
 double KC[2];
#define R R_COH
 double R = 0;
#define Tdur Tdur_COH
 double Tdur = 0.1;
#define Tamp Tamp_COH
 double Tamp = 1;
#define kd kd_COH
 double kd = 0.002;
#define ke ke_COH
 double ke = 80;
#define n0 n0_COH
 double n0 = 1;
#define usetable usetable_COH
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_COH", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "n0_COH", "1",
 "ke_COH", "/mM",
 "kd_COH", "/ms",
 "Cres_COH", "mM",
 "Cdur_COH", "ms",
 "Cnamp_COH", "mM",
 "Cnres_COH", "mM",
 "Cndur_COH", "ms",
 "Tamp_COH", "mM",
 "Tdur_COH", "ms",
 "R_COH", "1",
 "KO_COH", "/mM /ms",
 "KC_COH", "/ms",
 "n", "1",
 "Rmax", "1",
 "RO", "1",
 "CGLU", "mM",
 "C", "mM",
 "Cn", "mM",
 "pv", "1",
 "Ccount", "1",
 "Cncount", "1",
 "Tcnt", "1",
 "spike", "ms",
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
 static double RO0 = 0;
 static double Rmax0 = 0;
 static double Tcnt0 = 0;
 static double delta_t = 1;
 static double pv0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "n0_COH", &n0_COH,
 "ke_COH", &ke_COH,
 "kd_COH", &kd_COH,
 "Cres_COH", &Cres_COH,
 "Cdur_COH", &Cdur_COH,
 "Cnamp_COH", &Cnamp_COH,
 "Cnres_COH", &Cnres_COH,
 "Cndur_COH", &Cndur_COH,
 "Tamp_COH", &Tamp_COH,
 "Tdur_COH", &Tdur_COH,
 "R_COH", &R_COH,
 "usetable_COH", &usetable_COH,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 "KO_COH", KO_COH, 2,
 "KC_COH", KC_COH, 2,
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
"COH",
 0,
 "spike[10001]",
 "rseed",
 "PR0",
 "nzones",
 "ntot",
 "Ttot",
 "Camp",
 0,
 "n[1000]",
 "Rmax",
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
 	_p = nrn_prop_data_alloc(_mechtype, 20033, _prop);
 	/*initialize range parameters*/
  }
 	_prop->param = _p;
 	_prop->param_size = 20033;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 2, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static _initlists();
 _coh_reg() {
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
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 COH /Users/pbmanis/Desktop/Python/nrnlibrary/nrnlibrary/i386/coh.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_inf[2];
 static double *_t_fac[2];
 static double *_t_tau[2];
static int _reset;
static char *modelname = "Calyx of Held";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
static _f_rates();
static onerel();
static prel();
static rates();
static release();
 static _n_rates();
 
static int  release (  )  {
   if ( ( index < 10001.0 )  && ( ( t - spike [ ((int) index ) ] ) >= 0.0 )  && ( ( t - spike [ ((int) index ) ] ) < dt )  && ( C < Camp ) ) {
     C = Camp ;
     index = index + 1.0 ;
     tspike = t ;
     Ccount = Cdur / dt ;
     Ttot = 0.0 ;
     {int  _li ;for ( _li = 0 ; _li <= ((int) nzones ) - 1 ; _li ++ ) {
       relthisspike [ _li ] = 0.0 ;
       } }
     }
   prel ( _threadargs_ ) ;
   ntot = 0.0 ;
   {int  _li ;for ( _li = 0 ; _li <= ((int) nzones ) - 1 ; _li ++ ) {
     if ( unirand ( _threadargs_ ) < dt * km ) {
       n [ _li ] = n [ _li ] + 1.0 ;
       }
     if ( n [ _li ] > 0.0  && unirand ( _threadargs_ ) < dt * kd ) {
       n [ _li ] = n [ _li ] - 1.0 ;
       }
     if ( Cn > 0.0 ) {
       if ( unirand ( _threadargs_ ) < dt * ke * Cn ) {
         n [ _li ] = n [ _li ] + 1.0 ;
         }
       }
     onerel ( _threadargscomma_ ((double) _li ) ) ;
     if ( CGLU [ _li ] > 0.0 ) {
       Tcnt [ _li ] = Tcnt [ _li ] - 1.0 ;
       if ( Tcnt [ _li ] < 0.0 ) {
         CGLU [ _li ] = 0.0 ;
         }
       }
     ntot = ntot + n [ _li ] ;
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
 
static int  onerel (  _li )  
	double _li ;
 {
   if ( relthisspike [ ((int) _li ) ]  == 0.0 ) {
     if ( trel [ ((int) _li ) ] < tspike ) {
       ui = R * pv [ ((int) _li ) ] * n [ ((int) _li ) ] / Rmax ;
       if ( n [ ((int) _li ) ] > 0.0  && unirand ( _threadargs_ ) < ui * dt * 40.0 ) {
         n [ ((int) _li ) ] = n [ ((int) _li ) ] - 1.0 ;
         CGLU [ ((int) _li ) ] = Tamp ;
         Tcnt [ ((int) _li ) ] = Tdur / dt ;
         trel [ ((int) _li ) ] = t ;
         Ttot = Ttot + 1.0 ;
         relthisspike [ ((int) _li ) ] = 1.0 ;
         }
       else {
         relthisspike [ ((int) _li ) ] = 2.0 ;
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
 
static int _ode_count(_type)int _type; { hoc_execerror("COH", "cannot be used with CVODE");}

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
 for (_i=0; _i<2; _i++) RO[_i] = RO0;
  Rmax = Rmax0;
 for (_i=0; _i<1000; _i++) Tcnt[_i] = Tcnt0;
 for (_i=0; _i<1000; _i++) n[_i] = n0;
 for (_i=0; _i<1000; _i++) pv[_i] = pv0;
 {
   index = 0.0 ;
   n [ 0 ] = n0 ;
   CGLU [ 0 ] = 0.0 ;
   trel [ 0 ] = 0.0 ;
   relthisspike [ 0 ] = 0.0 ;
   pv [ 0 ] = PR0 ;
   {int  _li ;for ( _li = 0 ; _li <= ((int) nzones ) - 1 ; _li ++ ) {
     n [ _li ] = n0 ;
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
   prel ( _threadargs_ ) ;
   Rmax = 0.000275953 ;
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
 if(error){fprintf(stderr,"at line 143 in file coh.mod:\n    SOLVE release\n"); nrn_complain(_p); abort_run(error);}
 
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
