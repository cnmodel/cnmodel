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
#define gkbar _p[0]
#define ik _p[1]
#define m _p[2]
#define z _p[3]
#define h _p[4]
#define ek _p[5]
#define cai _p[6]
#define Dm _p[7]
#define Dz _p[8]
#define Dh _p[9]
#define _g _p[10]
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
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
 /* declaration of user functions */
 static int _hoc_rates();
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
 "setdata_bkpkj", _hoc_setdata,
 "rates_bkpkj", _hoc_rates,
 0, 0
};
 /* declare global and static user variables */
#define htau_k2 htau_k2_bkpkj
 double htau_k2 = 5.2;
#define htau_k1 htau_k1_bkpkj
 double htau_k1 = -12.9;
#define htau_vh2 htau_vh2_bkpkj
 double htau_vh2 = 48.5;
#define htau_vh1 htau_vh1_bkpkj
 double htau_vh1 = -54.2;
#define htau_y0 htau_y0_bkpkj
 double htau_y0 = 0.0019;
#define h_k h_k_bkpkj
 double h_k = 5.8;
#define h_vh h_vh_bkpkj
 double h_vh = -32;
#define h_y0 h_y0_bkpkj
 double h_y0 = 0.085;
#define htau htau_bkpkj
 double htau = 0;
#define hinf hinf_bkpkj
 double hinf = 0;
#define m1 m1_bkpkj
 double m1 = 0;
#define mtau_k2 mtau_k2_bkpkj
 double mtau_k2 = 10.1;
#define mtau_k1 mtau_k1_bkpkj
 double mtau_k1 = -10;
#define mtau_vh2 mtau_vh2_bkpkj
 double mtau_vh2 = 86.4;
#define mtau_vh1 mtau_vh1_bkpkj
 double mtau_vh1 = -33.3;
#define mtau_y0 mtau_y0_bkpkj
 double mtau_y0 = 0.000505;
#define m_k m_k_bkpkj
 double m_k = 6.2;
#define m_vh m_vh_bkpkj
 double m_vh = -28.9;
#define mtau mtau_bkpkj
 double mtau = 0;
#define minf minf_bkpkj
 double minf = 0;
#define z_coef z_coef_bkpkj
 double z_coef = 0.001;
#define ztau ztau_bkpkj
 double ztau = 1;
#define zinf zinf_bkpkj
 double zinf = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "m1_bkpkj", "mV",
 "m_vh_bkpkj", "mV",
 "m_k_bkpkj", "mV",
 "mtau_y0_bkpkj", "s",
 "mtau_vh1_bkpkj", "mV",
 "mtau_k1_bkpkj", "mV",
 "mtau_vh2_bkpkj", "mV",
 "mtau_k2_bkpkj", "mV",
 "z_coef_bkpkj", "mM",
 "ztau_bkpkj", "ms",
 "h_vh_bkpkj", "mV",
 "h_k_bkpkj", "mV",
 "htau_y0_bkpkj", "s",
 "htau_vh1_bkpkj", "mV",
 "htau_k1_bkpkj", "mV",
 "htau_vh2_bkpkj", "mV",
 "htau_k2_bkpkj", "mV",
 "mtau_bkpkj", "ms",
 "htau_bkpkj", "ms",
 "gkbar_bkpkj", "mho/cm2",
 "ik_bkpkj", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 static double v = 0;
 static double z0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "m1_bkpkj", &m1_bkpkj,
 "m_vh_bkpkj", &m_vh_bkpkj,
 "m_k_bkpkj", &m_k_bkpkj,
 "mtau_y0_bkpkj", &mtau_y0_bkpkj,
 "mtau_vh1_bkpkj", &mtau_vh1_bkpkj,
 "mtau_k1_bkpkj", &mtau_k1_bkpkj,
 "mtau_vh2_bkpkj", &mtau_vh2_bkpkj,
 "mtau_k2_bkpkj", &mtau_k2_bkpkj,
 "z_coef_bkpkj", &z_coef_bkpkj,
 "ztau_bkpkj", &ztau_bkpkj,
 "h_y0_bkpkj", &h_y0_bkpkj,
 "h_vh_bkpkj", &h_vh_bkpkj,
 "h_k_bkpkj", &h_k_bkpkj,
 "htau_y0_bkpkj", &htau_y0_bkpkj,
 "htau_vh1_bkpkj", &htau_vh1_bkpkj,
 "htau_k1_bkpkj", &htau_k1_bkpkj,
 "htau_vh2_bkpkj", &htau_vh2_bkpkj,
 "htau_k2_bkpkj", &htau_k2_bkpkj,
 "minf_bkpkj", &minf_bkpkj,
 "mtau_bkpkj", &mtau_bkpkj,
 "hinf_bkpkj", &hinf_bkpkj,
 "htau_bkpkj", &htau_bkpkj,
 "zinf_bkpkj", &zinf_bkpkj,
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
"bkpkj",
 "gkbar_bkpkj",
 0,
 "ik_bkpkj",
 0,
 "m_bkpkj",
 "z_bkpkj",
 "h_bkpkj",
 0,
 0};
 static Symbol* _k_sym;
 static Symbol* _ca_sym;
 
static void nrn_alloc(_prop)
	Prop *_prop;
{
	Prop *prop_ion, *need_memb();
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 11, _prop);
 	/*initialize range parameters*/
 	gkbar = 0.007;
 	_prop->param = _p;
 	_prop->param_size = 11;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
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
 _bkpkj_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("k", -10000.);
 	ion_reg("ca", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_dparam_size(_mechtype, 5);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 bkpkj /Users/pbmanis/Desktop/Python/nrnlibrary/nrnlibrary/i386/bkpkj.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
static rates();
 
static int _ode_spec1(), _ode_matsol1();
 static int _slist1[3], _dlist1[3];
 static int states();
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / mtau ;
   Dh = ( hinf - h ) / htau ;
   Dz = ( zinf - z ) / ztau ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau )) ;
 Dz = Dz  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ztau )) ;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / mtau)))*(- ( ( ( minf ) ) / mtau ) / ( ( ( ( - 1.0) ) ) / mtau ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / htau)))*(- ( ( ( hinf ) ) / htau ) / ( ( ( ( - 1.0) ) ) / htau ) - h) ;
    z = z + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ztau)))*(- ( ( ( zinf ) ) / ztau ) / ( ( ( ( - 1.0) ) ) / ztau ) - z) ;
   }
  return 0;
}
 
static int  rates (  _lVm )  
	double _lVm ;
 {
   double _lv ;
 _lv = _lVm + 5.0 ;
   minf = 1.0 / ( 1.0 + exp ( - ( _lv - ( m_vh ) ) / m_k ) ) ;
   m1 = mtau_y0 + ( 1.0 ) / ( exp ( ( _lv + mtau_vh1 ) / mtau_k1 ) ) ;
   mtau = ( 1e0 ) * ( m1 + exp ( ( _lv + mtau_vh2 ) / mtau_k2 ) ) ;
   zinf = 1.0 / ( 1.0 + z_coef / cai ) ;
   hinf = h_y0 + ( 1.0 - h_y0 ) / ( 1.0 + exp ( ( _lv - h_vh ) / h_k ) ) ;
   htau = ( 1e3 ) * ( htau_y0 + 1.0 / ( exp ( ( _lv + htau_vh1 ) / htau_k1 ) + exp ( ( _lv + htau_vh2 ) / htau_k2 ) ) ) ;
    return 0; }
 
static int _hoc_rates() {
  double _r;
   _r = 1.;
 rates (  *getarg(1) ) ;
 ret(_r);
}
 
static int _ode_count(_type) int _type;{ return 3;}
 
static int _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
  cai = _ion_cai;
     _ode_spec1 ();
  }}
 
static int _ode_map(_ieq, _pv, _pvdot, _pp, _ppd, _atol, _type) int _ieq, _type; double** _pv, **_pvdot, *_pp, *_atol; Datum* _ppd; { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
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
  ek = _ion_ek;
  cai = _ion_cai;
 _ode_matsol1 ();
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 3, 1);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
  z = z0;
 {
   rates ( _threadargscomma_ v ) ;
   m = minf ;
   z = zinf ;
   h = hinf ;
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
  cai = _ion_cai;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   ik = gkbar * m * m * m * z * z * h * ( v - ek ) ;
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
  cai = _ion_cai;
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
  cai = _ion_cai;
 { {
 for (; t < _break; t += dt) {
 error =  states();
 if(error){fprintf(stderr,"at line 67 in file bkpkj.mod:\n           SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 
}}
 t = _save;
 } }}

}

static terminal(){}

static _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
 _slist1[2] = &(z) - _p;  _dlist1[2] = &(Dz) - _p;
_first = 0;
}
