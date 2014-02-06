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
#define vshift _p[1]
#define tau_s _p[2]
#define tau_f _p[3]
#define tau_c _p[4]
#define s1 _p[5]
#define s2 _p[6]
#define f1 _p[7]
#define f2 _p[8]
#define cai _p[9]
#define eh _p[10]
#define Ds1 _p[11]
#define Ds2 _p[12]
#define Df1 _p[13]
#define Df2 _p[14]
#define ih _p[15]
#define h_inf _p[16]
#define alpha_s _p[17]
#define alpha_f _p[18]
#define beta_s _p[19]
#define beta_f _p[20]
#define C _p[21]
#define k2 _p[22]
#define tadj _p[23]
#define s0 _p[24]
#define f0 _p[25]
#define _g _p[26]
#define _ion_eh	*_ppvar[0]._pval
#define _ion_ih	*_ppvar[1]._pval
#define _ion_dihdv	*_ppvar[2]._pval
#define _ion_cai	*_ppvar[3]._pval
 
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
 static int _hoc_evaluate_fct();
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
 "setdata_iH_std", _hoc_setdata,
 "evaluate_fct_iH_std", _hoc_evaluate_fct,
 0, 0
};
 /* declare global and static user variables */
#define cac cac_iH_std
 double cac = 0.0005;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "cac_iH_std", "mM",
 "ghbar_iH_std", "mho/cm2",
 "vshift_iH_std", "mV",
 "tau_s_iH_std", "ms",
 "tau_f_iH_std", "ms",
 "tau_c_iH_std", "ms",
 0,0
};
 static double delta_t = 1;
 static double f20 = 0;
 static double f10 = 0;
 static double s20 = 0;
 static double s10 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "cac_iH_std", &cac_iH_std,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(), nrn_init(), nrn_state();
 static void nrn_cur(), nrn_jacob();
 
static int _ode_count(), _ode_map(), _ode_spec(), _ode_matsol();
 
#define _cvode_ieq _ppvar[4]._i
 /* connect range variables in _p that hoc is supposed to know about */
 static char *_mechanism[] = {
 "6.2.0",
"iH_std",
 "ghbar_iH_std",
 "vshift_iH_std",
 0,
 "tau_s_iH_std",
 "tau_f_iH_std",
 "tau_c_iH_std",
 0,
 "s1_iH_std",
 "s2_iH_std",
 "f1_iH_std",
 "f2_iH_std",
 0,
 0};
 static Symbol* _h_sym;
 static Symbol* _ca_sym;
 
static void nrn_alloc(_prop)
	Prop *_prop;
{
	Prop *prop_ion, *need_memb();
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 27, _prop);
 	/*initialize range parameters*/
 	ghbar = 4e-05;
 	vshift = 0;
 	_prop->param = _p;
 	_prop->param_size = 27;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_h_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* eh */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ih */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dihdv */
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[3]._pval = &prop_ion->param[1]; /* cai */
 
}
 static _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 _Ih_std_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("h", 1.0);
 	ion_reg("ca", -10000.);
 	_h_sym = hoc_lookup("h_ion");
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_dparam_size(_mechtype, 5);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 iH_std /Users/pbmanis/Desktop/Python/nrnlibrary/nrnlibrary/i386/Ih_std.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "hyperpolarization-activated current (H-current)";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
static evaluate_fct();
 
static int _ode_spec1(), _ode_matsol1();
 static double *_temp1;
 static int _slist1[4], _dlist1[4];
 static int states();
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   evaluate_fct ( _threadargscomma_ v , cai ) ;
   Ds1 = alpha_s * s0 - beta_s * s1 + k2 * ( s2 - C * s1 ) ;
   Df1 = alpha_f * f0 - beta_f * f1 + k2 * ( f2 - C * f1 ) ;
   Ds2 = - k2 * ( s2 - C * s1 ) ;
   Df2 = - k2 * ( f2 - C * f1 ) ;
   s0 = 1.0 - s1 - s2 ;
   f0 = 1.0 - f1 - f2 ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 evaluate_fct ( _threadargscomma_ v , cai ) ;
 Ds1 = Ds1  / (1. - dt*( ( - (beta_s)*(1.0) ) + (k2)*(( ( - (C)*(1.0) ) )) )) ;
 Df1 = Df1  / (1. - dt*( ( - (beta_f)*(1.0) ) + (k2)*(( ( - (C)*(1.0) ) )) )) ;
 Ds2 = Ds2  / (1. - dt*( (- k2)*(( 1.0 )) )) ;
 Df2 = Df2  / (1. - dt*( (- k2)*(( 1.0 )) )) ;
}
 /*END CVODE*/
 
static int states () {_reset=0;
 {
   evaluate_fct ( _threadargscomma_ v , cai ) ;
   Ds1 = alpha_s * s0 - beta_s * s1 + k2 * ( s2 - C * s1 ) ;
   Df1 = alpha_f * f0 - beta_f * f1 + k2 * ( f2 - C * f1 ) ;
   Ds2 = - k2 * ( s2 - C * s1 ) ;
   Df2 = - k2 * ( f2 - C * f1 ) ;
   s0 = 1.0 - s1 - s2 ;
   f0 = 1.0 - f1 - f2 ;
   }
 return _reset;}
 
static int  evaluate_fct (  _lv , _lcai )  
	double _lv , _lcai ;
 {
   h_inf = 1.0 / ( 1.0 + exp ( ( _lv + 88.9 + vshift ) / 7.8 ) ) ;
   tau_s = 0.25 * exp ( ( _lv + 183.6 + vshift ) / 15.24 ) / tadj ;
   tau_f = 0.25 * exp ( ( _lv + 158.6 + vshift ) / 11.2 ) / ( 1.0 + exp ( ( _lv + 75.0 + vshift ) / 5.5 ) ) / tadj ;
   alpha_s = h_inf / tau_s ;
   alpha_f = h_inf / tau_f ;
   beta_s = ( 1.0 - h_inf ) / tau_s ;
   beta_f = ( 1.0 - h_inf ) / tau_f ;
   C = _lcai * _lcai / ( cac * cac ) ;
   k2 = 4e-4 * tadj ;
    return 0; }
 
static int _hoc_evaluate_fct() {
  double _r;
   _r = 1.;
 evaluate_fct (  *getarg(1) , *getarg(2) ) ;
 ret(_r);
}
 
static int _ode_count(_type) int _type;{ return 4;}
 
static int _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  eh = _ion_eh;
  cai = _ion_cai;
     _ode_spec1 ();
  }}
 
static int _ode_map(_ieq, _pv, _pvdot, _pp, _ppd, _atol, _type) int _ieq, _type; double** _pv, **_pvdot, *_pp, *_atol; Datum* _ppd; { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 4; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static int _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  eh = _ion_eh;
  cai = _ion_cai;
 _ode_matsol1 ();
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_h_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_h_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_h_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 3, 1);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  f2 = f20;
  f1 = f10;
  s2 = s20;
  s1 = s10;
 {
   tadj = pow( 3.0 , ( ( celsius - 35.5 ) / 10.0 ) ) ;
   evaluate_fct ( _threadargscomma_ v , cai ) ;
   s1 = alpha_s / ( beta_s + alpha_s * ( 1.0 + C ) ) ;
   s2 = alpha_s * C / ( beta_s + alpha_s * ( 1.0 + C ) ) ;
   s0 = 1.0 - s1 - s2 ;
   f1 = alpha_f / ( beta_f + alpha_f * ( 1.0 + C ) ) ;
   f2 = alpha_f * C / ( beta_f + alpha_f * ( 1.0 + C ) ) ;
   f0 = 1.0 - f1 - f2 ;
   tau_c = 1.0 / ( 1.0 + C ) / k2 ;
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
  eh = _ion_eh;
  cai = _ion_cai;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   ih = ghbar * ( s1 + s2 ) * ( f1 + f2 ) * ( v - eh ) ;
   }
 _current += ih;

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
  eh = _ion_eh;
  cai = _ion_cai;
 _g = _nrn_current(_v + .001);
 	{ double _dih;
  _dih = ih;
 _rhs = _nrn_current(_v);
  _ion_dihdv += (_dih - ih)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ih += ih ;
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
  eh = _ion_eh;
  cai = _ion_cai;
 { {
 for (; t < _break; t += dt) {
 error =  euler(_ninits, 4, _slist1, _dlist1, _p, &t, dt, states, &_temp1);
 if(error){fprintf(stderr,"at line 86 in file Ih_std.mod:\n	SOLVE states METHOD euler\n"); nrn_complain(_p); abort_run(error);}
 
}}
 t = _save;
  states();
 } }}

}

static terminal(){}

static _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(s1) - _p;  _dlist1[0] = &(Ds1) - _p;
 _slist1[1] = &(f1) - _p;  _dlist1[1] = &(Df1) - _p;
 _slist1[2] = &(s2) - _p;  _dlist1[2] = &(Ds2) - _p;
 _slist1[3] = &(f2) - _p;  _dlist1[3] = &(Df2) - _p;
_first = 0;
}
