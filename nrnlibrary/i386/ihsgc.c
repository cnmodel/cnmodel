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
#define ghbar _p[0]
#define eh _p[1]
#define vh _p[2]
#define k _p[3]
#define vshift _p[4]
#define gh _p[5]
#define i _p[6]
#define ih _p[7]
#define rinf _p[8]
#define rtau _p[9]
#define sinf _p[10]
#define stau _p[11]
#define r _p[12]
#define s _p[13]
#define Dr _p[14]
#define Ds _p[15]
#define _g _p[16]
 
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
 extern double celsius;
 /* declaration of user functions */
 static int _hoc_rates();
 static int _hoc_states();
 static int _hoc_trates();
 static int _mechtype;
extern int nrn_get_mechtype();
 static _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range();
 _prop = hoc_getdata_range(_mechtype);
 _p = _prop->param; _ppvar = _prop->dparam;
 ret(1.);
}
 /* connect user functions to hoc names */
 static IntFunc hoc_intfunc[] = {
 "setdata_ihsgc", _hoc_setdata,
 "rates_ihsgc", _hoc_rates,
 "states_ihsgc", _hoc_states,
 "trates_ihsgc", _hoc_trates,
 0, 0
};
 /* declare global and static user variables */
#define aslow aslow_ihsgc
 double aslow = 0.465;
#define afast afast_ihsgc
 double afast = 0.535;
#define kfac2 kfac2_ihsgc
 double kfac2 = 11.2766;
#define kfac1 kfac1_ihsgc
 double kfac1 = 26.278;
#define skfac2 skfac2_ihsgc
 double skfac2 = 19.3753;
#define svtau2 svtau2_ihsgc
 double svtau2 = 13.9087;
#define stausc2 stausc2_ihsgc
 double stausc2 = 1.12963e-05;
#define skfac1 skfac1_ihsgc
 double skfac1 = 37.3538;
#define svtau1 svtau1_ihsgc
 double svtau1 = 30.7426;
#define stausc1 stausc1_ihsgc
 double stausc1 = 0.00792101;
#define staumin staumin_ihsgc
 double staumin = 40.7778;
#define staufac staufac_ihsgc
 double staufac = 1;
#define sba2 sba2_ihsgc
 double sba2 = 0.295536;
#define sba1 sba1_ihsgc
 double sba1 = 0.704464;
#define svshift svshift_ihsgc
 double svshift = 0;
#define sk2 sk2_ihsgc
 double sk2 = -5.281;
#define svh2 svh2_ihsgc
 double svh2 = -116.164;
#define sk1 sk1_ihsgc
 double sk1 = 2.343;
#define svh1 svh1_ihsgc
 double svh1 = -87.981;
#define tausc2 tausc2_ihsgc
 double tausc2 = 2.65845e-05;
#define tausc1 tausc1_ihsgc
 double tausc1 = 0.0174863;
#define taumin taumin_ihsgc
 double taumin = 29.3334;
#define taufac taufac_ihsgc
 double taufac = 1;
#define usetable usetable_ihsgc
 double usetable = 1;
#define vtau2 vtau2_ihsgc
 double vtau2 = 30.7333;
#define vtau1 vtau1_ihsgc
 double vtau1 = 61.2494;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "ghbar_ihsgc", 0, 1e+09,
 "usetable_ihsgc", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "taufac_ihsgc", "1",
 "taumin_ihsgc", "ms",
 "tausc1_ihsgc", "ms",
 "vtau1_ihsgc", "mV",
 "kfac1_ihsgc", "mV",
 "tausc2_ihsgc", "ms",
 "vtau2_ihsgc", "mV",
 "kfac2_ihsgc", "mV",
 "svh1_ihsgc", "mV",
 "sk1_ihsgc", "mV",
 "svh2_ihsgc", "mV",
 "sk2_ihsgc", "mV",
 "svshift_ihsgc", "mV",
 "staufac_ihsgc", "1",
 "staumin_ihsgc", "ms",
 "stausc1_ihsgc", "ms",
 "svtau1_ihsgc", "mV",
 "skfac1_ihsgc", "mV",
 "stausc2_ihsgc", "ms",
 "svtau2_ihsgc", "mV",
 "skfac2_ihsgc", "mV",
 "ghbar_ihsgc", "mho/cm2",
 "eh_ihsgc", "mV",
 "vh_ihsgc", "mV",
 "k_ihsgc", "mV",
 "vshift_ihsgc", "mV",
 "gh_ihsgc", "mho/cm2",
 "i_ihsgc", "mA/cm2",
 "ih_ihsgc", "mA/cm2",
 "rtau_ihsgc", "ms",
 "stau_ihsgc", "ms",
 0,0
};
 static double delta_t = 1;
 static double r0 = 0;
 static double s0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "afast_ihsgc", &afast_ihsgc,
 "taufac_ihsgc", &taufac_ihsgc,
 "taumin_ihsgc", &taumin_ihsgc,
 "tausc1_ihsgc", &tausc1_ihsgc,
 "vtau1_ihsgc", &vtau1_ihsgc,
 "kfac1_ihsgc", &kfac1_ihsgc,
 "tausc2_ihsgc", &tausc2_ihsgc,
 "vtau2_ihsgc", &vtau2_ihsgc,
 "kfac2_ihsgc", &kfac2_ihsgc,
 "svh1_ihsgc", &svh1_ihsgc,
 "sk1_ihsgc", &sk1_ihsgc,
 "svh2_ihsgc", &svh2_ihsgc,
 "sk2_ihsgc", &sk2_ihsgc,
 "svshift_ihsgc", &svshift_ihsgc,
 "sba1_ihsgc", &sba1_ihsgc,
 "sba2_ihsgc", &sba2_ihsgc,
 "aslow_ihsgc", &aslow_ihsgc,
 "staufac_ihsgc", &staufac_ihsgc,
 "staumin_ihsgc", &staumin_ihsgc,
 "stausc1_ihsgc", &stausc1_ihsgc,
 "svtau1_ihsgc", &svtau1_ihsgc,
 "skfac1_ihsgc", &skfac1_ihsgc,
 "stausc2_ihsgc", &stausc2_ihsgc,
 "svtau2_ihsgc", &svtau2_ihsgc,
 "skfac2_ihsgc", &skfac2_ihsgc,
 "usetable_ihsgc", &usetable_ihsgc,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(), nrn_init(), nrn_state();
 static void nrn_cur(), nrn_jacob();
 
static int _ode_count();
 /* connect range variables in _p that hoc is supposed to know about */
 static char *_mechanism[] = {
 "6.2.0",
"ihsgc",
 "ghbar_ihsgc",
 "eh_ihsgc",
 "vh_ihsgc",
 "k_ihsgc",
 "vshift_ihsgc",
 0,
 "gh_ihsgc",
 "i_ihsgc",
 "ih_ihsgc",
 "rinf_ihsgc",
 "rtau_ihsgc",
 "sinf_ihsgc",
 "stau_ihsgc",
 0,
 "r_ihsgc",
 "s_ihsgc",
 0,
 0};
 
static void nrn_alloc(_prop)
	Prop *_prop;
{
	Prop *prop_ion, *need_memb();
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 17, _prop);
 	/*initialize range parameters*/
 	ghbar = 0.00318;
 	eh = -43;
 	vh = -103.035;
 	k = 11.223;
 	vshift = 0;
 	_prop->param = _p;
 	_prop->param_size = 17;
 
}
 static _initlists();
 _ihsgc_reg() {
	int _vectorized = 0;
  _initlists();
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
  hoc_register_dparam_size(_mechtype, 0);
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 ihsgc /Users/pbmanis/Desktop/Python/nrnlibrary/nrnlibrary/i386/ihsgc.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double _zrexp , _zsexp ;
 static double _zq10 , _zst ;
 static double *_t_rinf;
 static double *_t__zrexp;
 static double *_t_sinf;
 static double *_t__zsexp;
static int _reset;
static char *modelname = "ihsgc.mod - Spiral Ganglion Cell Ih current";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
static _f_trates();
static rates();
static states();
static trates();
 static _n_trates();
 
static int  states (  )  {
   trates ( _threadargscomma_ v ) ;
   r = r + _zrexp * ( rinf - r ) ;
   s = s + _zsexp * ( sinf - s ) ;
   
/*VERBATIM*/
	return 0;
  return 0; }
 
static int _hoc_states() {
  double _r;
   _r = 1.;
 states (  ) ;
 ret(_r);
}
 
static int  rates (  _lv )  
	double _lv ;
 {
   _zq10 = pow( 3.0 , ( ( celsius - 22.0 ) / 10.0 ) ) ;
   rinf = 1.0 / pow( ( 1.0 + exp ( ( _lv - vh + vshift ) / k ) ) , 0.5 ) ;
   rtau = tausc1 * exp ( ( _lv + vtau1 + vshift ) / kfac1 ) + tausc2 * exp ( - ( _lv + vtau2 + vshift ) / kfac2 ) ;
   rtau = ( taumin + taufac / rtau ) ;
   sinf = 1. / ( 1.0 + exp ( ( _lv - svh1 + vshift ) / sk1 ) ) ;
   _zst = 1. / ( 1.0 + exp ( ( _lv - svh2 + vshift ) / sk2 ) ) ;
   sinf = sinf + sba2 * ( _zst - 1.0 ) ;
   stau = staufac / ( stausc1 * exp ( ( _lv + svtau1 + vshift ) / skfac1 ) + stausc2 * exp ( - ( _lv + svtau2 + vshift ) / skfac2 ) ) ;
   stau = ( stau + staumin ) ;
    return 0; }
 
static int _hoc_rates() {
  double _r;
   _r = 1.;
 rates (  *getarg(1) ) ;
 ret(_r);
}
 static double _mfac_trates, _tmin_trates;
 static _check_trates();
 static _check_trates() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_dt;
  static double _sav_celsius;
  if (!usetable) {return;}
  if (_sav_dt != dt) { _maktable = 1;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_trates =  - 200.0 ;
   _tmax =  150.0 ;
   _dx = (_tmax - _tmin_trates)/350.; _mfac_trates = 1./_dx;
   for (_i=0, _x=_tmin_trates; _i < 351; _x += _dx, _i++) {
    _f_trates(_x);
    _t_rinf[_i] = rinf;
    _t__zrexp[_i] = _zrexp;
    _t_sinf[_i] = sinf;
    _t__zsexp[_i] = _zsexp;
   }
   _sav_dt = dt;
   _sav_celsius = celsius;
  }
 }

 static trates(double _lv){ _check_trates();
 _n_trates(_lv);
 return;
 }

 static _n_trates(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_trates(_lv); return; 
}
 _xi = _mfac_trates * (_lv - _tmin_trates);
 _i = (int) _xi;
 if (_xi <= 0.) {
 rinf = _t_rinf[0];
 _zrexp = _t__zrexp[0];
 sinf = _t_sinf[0];
 _zsexp = _t__zsexp[0];
 return; }
 if (_i >= 350) {
 rinf = _t_rinf[350];
 _zrexp = _t__zrexp[350];
 sinf = _t_sinf[350];
 _zsexp = _t__zsexp[350];
 return; }
 _theta = _xi - (double)_i;
 rinf = _t_rinf[_i] + _theta*(_t_rinf[_i+1] - _t_rinf[_i]);
 _zrexp = _t__zrexp[_i] + _theta*(_t__zrexp[_i+1] - _t__zrexp[_i]);
 sinf = _t_sinf[_i] + _theta*(_t_sinf[_i+1] - _t_sinf[_i]);
 _zsexp = _t__zsexp[_i] + _theta*(_t__zsexp[_i+1] - _t__zsexp[_i]);
 }

 
static int  _f_trates (  _lv )  
	double _lv ;
 {
   double _ltinc ;
 rates ( _threadargscomma_ _lv ) ;
   _ltinc = - dt * _zq10 ;
   _zrexp = 1.0 - exp ( _ltinc / rtau ) ;
   _zsexp = 1.0 - exp ( _ltinc / stau ) ;
    return 0; }
 
static int _hoc_trates() {
  double _r;
    _r = 1.;
 trates (  *getarg(1) ) ;
 ret(_r);
}
 
static int _ode_count(_type)int _type; { hoc_execerror("ihsgc", "cannot be used with CVODE");}

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  r = r0;
  s = s0;
 {
   trates ( _threadargscomma_ v ) ;
   r = rinf ;
   s = sinf ;
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

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   gh = ghbar * ( afast * ( pow( r , 2.0 ) ) + aslow * s ) ;
   ih = gh * ( v - eh ) ;
   i = ih ;
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
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_v + .001);
 	{ _rhs = _nrn_current(_v);
 	}
 _g = (_g - _rhs)/.001;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
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
 error =  states();
 if(error){fprintf(stderr,"at line 117 in file ihsgc.mod:\n	SOLVE states\n"); nrn_complain(_p); abort_run(error);}
 
}}
 t = _save;
 }}}

}

static terminal(){}

static _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
   _t_rinf = makevector(351*sizeof(double));
   _t__zrexp = makevector(351*sizeof(double));
   _t_sinf = makevector(351*sizeof(double));
   _t__zsexp = makevector(351*sizeof(double));
_first = 0;
}
