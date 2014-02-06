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
#define gkhtbar _p[0]
#define ik _p[1]
#define gkht _p[2]
#define n _p[3]
#define p _p[4]
#define ek _p[5]
#define Dn _p[6]
#define Dp _p[7]
#define _g _p[8]
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
 
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
 static int _hoc_vtrap();
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
 "setdata_kht", _hoc_setdata,
 "rates_kht", _hoc_rates,
 "states_kht", _hoc_states,
 "trates_kht", _hoc_trates,
 "vtrap_kht", _hoc_vtrap,
 0, 0
};
#define vtrap vtrap_kht
 extern double vtrap();
 /* declare global and static user variables */
#define nf nf_kht
 double nf = 0.85;
#define ntau ntau_kht
 double ntau = 0;
#define ninf ninf_kht
 double ninf = 0;
#define ptau ptau_kht
 double ptau = 0;
#define pinf pinf_kht
 double pinf = 0;
#define usetable usetable_kht
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "gkhtbar_kht", 0, 1e+09,
 "nf_kht", 0, 1,
 "usetable_kht", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "ptau_kht", "ms",
 "ntau_kht", "ms",
 "gkhtbar_kht", "mho/cm2",
 "ik_kht", "mA/cm",
 "gkht_kht", "mho/cm2",
 0,0
};
 static double delta_t = 1;
 static double n0 = 0;
 static double p0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "nf_kht", &nf_kht,
 "pinf_kht", &pinf_kht,
 "ninf_kht", &ninf_kht,
 "ptau_kht", &ptau_kht,
 "ntau_kht", &ntau_kht,
 "usetable_kht", &usetable_kht,
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
"kht",
 "gkhtbar_kht",
 0,
 "ik_kht",
 "gkht_kht",
 0,
 "n_kht",
 "p_kht",
 0,
 0};
 static Symbol* _k_sym;
 
static void nrn_alloc(_prop)
	Prop *_prop;
{
	Prop *prop_ion, *need_memb();
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 9, _prop);
 	/*initialize range parameters*/
 	gkhtbar = 0.01592;
 	_prop->param = _p;
 	_prop->param_size = 9;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static _initlists();
 static void _update_ion_pointer(Datum*);
 _kht_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("k", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_dparam_size(_mechtype, 3);
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 kht /Users/pbmanis/Desktop/Python/nrnlibrary/nrnlibrary/i386/kht.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double _znexp , _zpexp ;
 static double _zq10 ;
 static double *_t_ninf;
 static double *_t__znexp;
 static double *_t_pinf;
 static double *_t__zpexp;
static int _reset;
static char *modelname = "kht.mod  The high threshold conductance of cochlear nucleus neurons";

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
   n = n + _znexp * ( ninf - n ) ;
   p = p + _zpexp * ( pinf - p ) ;
   
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
   ninf = pow( ( 1.0 + exp ( - ( _lv + 15.0 ) / 5.0 ) ) , - 0.5 ) ;
   pinf = 1.0 / ( 1.0 + exp ( - ( _lv + 23.0 ) / 6.0 ) ) ;
   ntau = ( 100.0 / ( 11.0 * exp ( ( _lv + 60.0 ) / 24.0 ) + 21.0 * exp ( - ( _lv + 60.0 ) / 23.0 ) ) ) + 0.7 ;
   ptau = ( 100.0 / ( 4.0 * exp ( ( _lv + 60.0 ) / 32.0 ) + 5.0 * exp ( - ( _lv + 60.0 ) / 22.0 ) ) ) + 5.0 ;
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
   _tmin_trates =  - 150.0 ;
   _tmax =  150.0 ;
   _dx = (_tmax - _tmin_trates)/300.; _mfac_trates = 1./_dx;
   for (_i=0, _x=_tmin_trates; _i < 301; _x += _dx, _i++) {
    _f_trates(_x);
    _t_ninf[_i] = ninf;
    _t__znexp[_i] = _znexp;
    _t_pinf[_i] = pinf;
    _t__zpexp[_i] = _zpexp;
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
 ninf = _t_ninf[0];
 _znexp = _t__znexp[0];
 pinf = _t_pinf[0];
 _zpexp = _t__zpexp[0];
 return; }
 if (_i >= 300) {
 ninf = _t_ninf[300];
 _znexp = _t__znexp[300];
 pinf = _t_pinf[300];
 _zpexp = _t__zpexp[300];
 return; }
 _theta = _xi - (double)_i;
 ninf = _t_ninf[_i] + _theta*(_t_ninf[_i+1] - _t_ninf[_i]);
 _znexp = _t__znexp[_i] + _theta*(_t__znexp[_i+1] - _t__znexp[_i]);
 pinf = _t_pinf[_i] + _theta*(_t_pinf[_i+1] - _t_pinf[_i]);
 _zpexp = _t__zpexp[_i] + _theta*(_t__zpexp[_i+1] - _t__zpexp[_i]);
 }

 
static int  _f_trates (  _lv )  
	double _lv ;
 {
   double _ltinc ;
 rates ( _threadargscomma_ _lv ) ;
   _ltinc = - dt * _zq10 ;
   _znexp = 1.0 - exp ( _ltinc / ntau ) ;
   _zpexp = 1.0 - exp ( _ltinc / ptau ) ;
    return 0; }
 
static int _hoc_trates() {
  double _r;
    _r = 1.;
 trates (  *getarg(1) ) ;
 ret(_r);
}
 
double vtrap (  _lx , _ly )  
	double _lx , _ly ;
 {
   double _lvtrap;
 if ( fabs ( _lx / _ly ) < 1e-6 ) {
     _lvtrap = _ly * ( 1.0 - _lx / _ly / 2.0 ) ;
     }
   else {
     _lvtrap = _lx / ( exp ( _lx / _ly ) - 1.0 ) ;
     }
   
return _lvtrap;
 }
 
static int _hoc_vtrap() {
  double _r;
   _r =  vtrap (  *getarg(1) , *getarg(2) ) ;
 ret(_r);
}
 
static int _ode_count(_type)int _type; { hoc_execerror("kht", "cannot be used with CVODE");}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  n = n0;
  p = p0;
 {
   trates ( _threadargscomma_ v ) ;
   p = pinf ;
   n = ninf ;
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
  ek = _ion_ek;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   gkht = gkhtbar * ( nf * ( pow( n , 2.0 ) ) + ( 1.0 - nf ) * p ) ;
   ik = gkht * ( v - ek ) ;
   }
 _current += ik;

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
  ek = _ion_ek;
 _g = _nrn_current(_v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
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
  ek = _ion_ek;
 { {
 for (; t < _break; t += dt) {
 error =  states();
 if(error){fprintf(stderr,"at line 74 in file kht.mod:\n    \n"); nrn_complain(_p); abort_run(error);}
 
}}
 t = _save;
 } }}

}

static terminal(){}

static _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
   _t_ninf = makevector(301*sizeof(double));
   _t__znexp = makevector(301*sizeof(double));
   _t_pinf = makevector(301*sizeof(double));
   _t__zpexp = makevector(301*sizeof(double));
_first = 0;
}
