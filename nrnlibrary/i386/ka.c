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
#define gkabar _p[0]
#define ik _p[1]
#define gka _p[2]
#define a _p[3]
#define b _p[4]
#define c _p[5]
#define ek _p[6]
#define Da _p[7]
#define Db _p[8]
#define Dc _p[9]
#define _g _p[10]
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
 "setdata_ka", _hoc_setdata,
 "rates_ka", _hoc_rates,
 "states_ka", _hoc_states,
 "trates_ka", _hoc_trates,
 "vtrap_ka", _hoc_vtrap,
 0, 0
};
#define vtrap vtrap_ka
 extern double vtrap();
 /* declare global and static user variables */
#define atau atau_ka
 double atau = 0;
#define ainf ainf_ka
 double ainf = 0;
#define btau btau_ka
 double btau = 0;
#define binf binf_ka
 double binf = 0;
#define ctau ctau_ka
 double ctau = 0;
#define cinf cinf_ka
 double cinf = 0;
#define usetable usetable_ka
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "gkabar_ka", 0, 1e+09,
 "usetable_ka", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "atau_ka", "ms",
 "btau_ka", "ms",
 "ctau_ka", "ms",
 "gkabar_ka", "mho/cm2",
 "ik_ka", "mA/cm2",
 "gka_ka", "mho/cm2",
 0,0
};
 static double a0 = 0;
 static double b0 = 0;
 static double c0 = 0;
 static double delta_t = 1;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "ainf_ka", &ainf_ka,
 "binf_ka", &binf_ka,
 "cinf_ka", &cinf_ka,
 "atau_ka", &atau_ka,
 "btau_ka", &btau_ka,
 "ctau_ka", &ctau_ka,
 "usetable_ka", &usetable_ka,
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
"ka",
 "gkabar_ka",
 0,
 "ik_ka",
 "gka_ka",
 0,
 "a_ka",
 "b_ka",
 "c_ka",
 0,
 0};
 static Symbol* _k_sym;
 
static void nrn_alloc(_prop)
	Prop *_prop;
{
	Prop *prop_ion, *need_memb();
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 11, _prop);
 	/*initialize range parameters*/
 	gkabar = 0.00477;
 	_prop->param = _p;
 	_prop->param_size = 11;
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
 _ka_reg() {
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
 	ivoc_help("help ?1 ka /Users/pbmanis/Desktop/Python/nrnlibrary/nrnlibrary/i386/ka.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double _zaexp , _zbexp , _zcexp ;
 static double _zq10 ;
 static double *_t_ainf;
 static double *_t__zaexp;
 static double *_t_binf;
 static double *_t__zbexp;
 static double *_t_cinf;
 static double *_t__zcexp;
static int _reset;
static char *modelname = "klt.mod  The low threshold conductance of cochlear nucleus neurons";

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
   a = a + _zaexp * ( ainf - a ) ;
   b = b + _zbexp * ( binf - b ) ;
   c = c + _zcexp * ( cinf - c ) ;
   
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
   ainf = pow( ( 1.0 / ( 1.0 + exp ( - 1.0 * ( _lv + 31.0 ) / 6.0 ) ) ) , 0.25 ) ;
   binf = 1.0 / pow( ( 1.0 + exp ( ( _lv + 66.0 ) / 7.0 ) ) , 0.5 ) ;
   cinf = 1.0 / pow( ( 1.0 + exp ( ( _lv + 66.0 ) / 7.0 ) ) , 0.5 ) ;
   atau = ( 100.0 / ( 7.0 * exp ( ( _lv + 60.0 ) / 14.0 ) + 29.0 * exp ( - ( _lv + 60.0 ) / 24.0 ) ) ) + 0.1 ;
   btau = ( 1000.0 / ( 14.0 * exp ( ( _lv + 60.0 ) / 27.0 ) + 29.0 * exp ( - ( _lv + 60.0 ) / 24.0 ) ) ) + 1.0 ;
   ctau = ( 90.0 / ( 1.0 + exp ( ( - 66.0 - _lv ) / 17.0 ) ) ) + 10.0 ;
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
    _t_ainf[_i] = ainf;
    _t__zaexp[_i] = _zaexp;
    _t_binf[_i] = binf;
    _t__zbexp[_i] = _zbexp;
    _t_cinf[_i] = cinf;
    _t__zcexp[_i] = _zcexp;
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
 ainf = _t_ainf[0];
 _zaexp = _t__zaexp[0];
 binf = _t_binf[0];
 _zbexp = _t__zbexp[0];
 cinf = _t_cinf[0];
 _zcexp = _t__zcexp[0];
 return; }
 if (_i >= 300) {
 ainf = _t_ainf[300];
 _zaexp = _t__zaexp[300];
 binf = _t_binf[300];
 _zbexp = _t__zbexp[300];
 cinf = _t_cinf[300];
 _zcexp = _t__zcexp[300];
 return; }
 _theta = _xi - (double)_i;
 ainf = _t_ainf[_i] + _theta*(_t_ainf[_i+1] - _t_ainf[_i]);
 _zaexp = _t__zaexp[_i] + _theta*(_t__zaexp[_i+1] - _t__zaexp[_i]);
 binf = _t_binf[_i] + _theta*(_t_binf[_i+1] - _t_binf[_i]);
 _zbexp = _t__zbexp[_i] + _theta*(_t__zbexp[_i+1] - _t__zbexp[_i]);
 cinf = _t_cinf[_i] + _theta*(_t_cinf[_i+1] - _t_cinf[_i]);
 _zcexp = _t__zcexp[_i] + _theta*(_t__zcexp[_i+1] - _t__zcexp[_i]);
 }

 
static int  _f_trates (  _lv )  
	double _lv ;
 {
   double _ltinc ;
 rates ( _threadargscomma_ _lv ) ;
   _ltinc = - dt * _zq10 ;
   _zaexp = 1.0 - exp ( _ltinc / atau ) ;
   _zbexp = 1.0 - exp ( _ltinc / btau ) ;
   _zcexp = 1.0 - exp ( _ltinc / ctau ) ;
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
 
static int _ode_count(_type)int _type; { hoc_execerror("ka", "cannot be used with CVODE");}
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
  a = a0;
  b = b0;
  c = c0;
 {
   trates ( _threadargscomma_ v ) ;
   a = ainf ;
   b = binf ;
   c = cinf ;
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
   gka = gkabar * ( pow( a , 4.0 ) ) * b * c ;
   ik = gka * ( v - ek ) ;
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
 if(error){fprintf(stderr,"at line 61 in file ka.mod:\n    \n"); nrn_complain(_p); abort_run(error);}
 
}}
 t = _save;
 } }}

}

static terminal(){}

static _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
   _t_ainf = makevector(301*sizeof(double));
   _t__zaexp = makevector(301*sizeof(double));
   _t_binf = makevector(301*sizeof(double));
   _t__zbexp = makevector(301*sizeof(double));
   _t_cinf = makevector(301*sizeof(double));
   _t__zcexp = makevector(301*sizeof(double));
_first = 0;
}
