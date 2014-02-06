/* Created by Language version: 6.2.0 */
/* VECTORIZED */
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
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define etaalpha _p[0]
#define etabeta _p[1]
#define kalpha _p[2]
#define kbeta _p[3]
#define kfs _p[4]
#define kbs _p[5]
#define gkna _p[6]
#define ninf _p[7]
#define sinf _p[8]
#define alpha _p[9]
#define beta _p[10]
#define stau _p[11]
#define n _p[12]
#define s _p[13]
#define nai _p[14]
#define ek _p[15]
#define ina _p[16]
#define ik _p[17]
#define Dn _p[18]
#define Ds _p[19]
#define v _p[20]
#define _g _p[21]
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
#define _ion_nai	*_ppvar[3]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static int _hoc_Kna_s();
 static int _hoc_Kna_n();
 static int _hoc_knabeta();
 static int _hoc_knalpha();
 static int _hoc_rates();
 static int _hoc_tau();
 static int _mechtype;
extern int nrn_get_mechtype();
 static _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range();
 _prop = hoc_getdata_range(_mechtype);
 _extcall_prop = _prop;
 ret(1.);
}
 /* connect user functions to hoc names */
 static IntFunc hoc_intfunc[] = {
 "setdata_kna", _hoc_setdata,
 "Kna_s_kna", _hoc_Kna_s,
 "Kna_n_kna", _hoc_Kna_n,
 "knabeta_kna", _hoc_knabeta,
 "knalpha_kna", _hoc_knalpha,
 "rates_kna", _hoc_rates,
 "tau_kna", _hoc_tau,
 0, 0
};
#define Kna_s Kna_s_kna
#define Kna_n Kna_n_kna
#define knabeta knabeta_kna
#define knalpha knalpha_kna
#define tau tau_kna
 extern double Kna_s();
 extern double Kna_n();
 extern double knabeta();
 extern double knalpha();
 extern double tau();
 /* declare global and static user variables */
#define gknabar gknabar_kna
 double gknabar = 1e-06;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gknabar_kna", "mho/cm2",
 "etaalpha_kna", "/mV",
 "etabeta_kna", "/mV",
 "kalpha_kna", "/ms",
 "kbeta_kna", "/ms",
 "kfs_kna", "/mM",
 "kbs_kna", "/mM",
 "gkna_kna", "mho/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double n0 = 0;
 static double s0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "gknabar_kna", &gknabar_kna,
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
"kna",
 "etaalpha_kna",
 "etabeta_kna",
 "kalpha_kna",
 "kbeta_kna",
 "kfs_kna",
 "kbs_kna",
 0,
 "gkna_kna",
 "ninf_kna",
 "sinf_kna",
 "alpha_kna",
 "beta_kna",
 "stau_kna",
 0,
 "n_kna",
 "s_kna",
 0,
 0};
 static Symbol* _k_sym;
 static Symbol* _na_sym;
 
static void nrn_alloc(_prop)
	Prop *_prop;
{
	Prop *prop_ion, *need_memb();
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 22, _prop);
 	/*initialize range parameters*/
 	etaalpha = 0.0558;
 	etabeta = 0.0308;
 	kalpha = 0.141;
 	kbeta = 0.0587;
 	kfs = 0.14;
 	kbs = 2;
 	_prop->param = _p;
 	_prop->param_size = 22;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[3]._pval = &prop_ion->param[1]; /* nai */
 
}
 static _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 _KNa_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("k", -10000.);
 	ion_reg("na", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_dparam_size(_mechtype, 5);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 kna /Users/pbmanis/Desktop/Python/nrnlibrary/nrnlibrary/i386/KNa.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "kna.mod";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
static rates();
 
static int _ode_spec1(), _ode_matsol1();
 static int _slist1[2], _dlist1[2];
 static int states();
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   rates ( _threadargscomma_ v , nai ) ;
   Dn = ( alpha * ( 1.0 - ninf ) ) - beta ;
   Ds = ( sinf - s ) / stau ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 rates ( _threadargscomma_ v , nai ) ;
 Dn = Dn  / (1. - dt*( 0.0 )) ;
 Ds = Ds  / (1. - dt*( ( ( ( - 1.0 ) ) ) / stau )) ;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   rates ( _threadargscomma_ v , nai ) ;
    n = n - dt*(- ( ( (alpha)*(( 1.0 - ninf )) ) - beta ) ) ;
    s = s + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / stau)))*(- ( ( ( sinf ) ) / stau ) / ( ( ( ( - 1.0) ) ) / stau ) - s) ;
   }
  return 0;
}
 
static int  rates ( _p, _ppvar, _thread, _nt, _lv , _lnai ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lv , _lnai ;
 {
   alpha = knalpha ( _threadargscomma_ _lv ) ;
   beta = knabeta ( _threadargscomma_ _lv ) ;
   ninf = Kna_n ( _threadargscomma_ alpha , beta ) ;
   sinf = Kna_s ( _threadargscomma_ _lnai ) ;
   stau = tau ( _threadargscomma_ _lnai ) ;
    return 0; }
 
static int _hoc_rates() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) ) ;
 ret(_r);
}
 
double Kna_s ( _p, _ppvar, _thread, _nt, _lnai ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lnai ;
 {
   double _lKna_s;
 _lKna_s = ( kfs * _lnai ) / ( kbs + ( kfs * _lnai ) ) ;
   
return _lKna_s;
 }
 
static int _hoc_Kna_s() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  Kna_s ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double tau ( _p, _ppvar, _thread, _nt, _lnai ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lnai ;
 {
   double _ltau;
 _ltau = 1.0 / ( kbs + ( kfs * _lnai ) ) ;
   
return _ltau;
 }
 
static int _hoc_tau() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  tau ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double Kna_n ( _p, _ppvar, _thread, _nt, _lalpha , _lbeta ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lalpha , _lbeta ;
 {
   double _lKna_n;
 _lKna_n = _lalpha / ( _lalpha + _lbeta ) ;
   
return _lKna_n;
 }
 
static int _hoc_Kna_n() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  Kna_n ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) ) ;
 ret(_r);
}
 
double knalpha ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lknalpha;
 _lknalpha = kalpha * ( exp ( etaalpha * _lx ) ) ;
   
return _lknalpha;
 }
 
static int _hoc_knalpha() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  knalpha ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double knabeta ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lknabeta;
 _lknabeta = kbeta * ( exp ( etabeta * _lx ) ) ;
   
return _lknabeta;
 }
 
static int _hoc_knabeta() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  knabeta ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
static int _ode_count(_type) int _type;{ return 2;}
 
static int _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
  nai = _ion_nai;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static int _ode_map(_ieq, _pv, _pvdot, _pp, _ppd, _atol, _type) int _ieq, _type; double** _pv, **_pvdot, *_pp, *_atol; Datum* _ppd; { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static int _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
  nai = _ion_nai;
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_na_sym, _ppvar, 3, 1);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  n = n0;
  s = s0;
 {
   rates ( _threadargscomma_ v , nai ) ;
   s = sinf ;
   n = ninf ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
  nai = _ion_nai;
 initmodel(_p, _ppvar, _thread, _nt);
 }}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   gkna = gknabar * n * n * s * s * s * s ;
   ik = gkna * ( v + 80.0 ) ;
   }
 _current += ik;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
  nai = _ion_nai;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
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

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
 double _break, _save;
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
  nai = _ion_nai;
 { {
 for (; t < _break; t += dt) {
   states(_p, _ppvar, _thread, _nt);
  
}}
 t = _save;
 } }}

}

static terminal(){}

static _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(n) - _p;  _dlist1[0] = &(Dn) - _p;
 _slist1[1] = &(s) - _p;  _dlist1[1] = &(Ds) - _p;
_first = 0;
}
