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
#define inakmax _p[0]
#define Nai_inf _p[1]
#define Ki_inf _p[2]
#define ATPi _p[3]
#define inak _p[4]
#define y _p[5]
#define nao _p[6]
#define ko _p[7]
#define ina _p[8]
#define ik _p[9]
#define nai _p[10]
#define ki _p[11]
#define pE1Na _p[12]
#define pE1K _p[13]
#define pE2Na _p[14]
#define pE2K _p[15]
#define alfay _p[16]
#define betay _p[17]
#define Naeff _p[18]
#define k1 _p[19]
#define drugblock _p[20]
#define kelvin _p[21]
#define Dy _p[22]
#define v _p[23]
#define _g _p[24]
#define _ion_nai	*_ppvar[0]._pval
#define _ion_nao	*_ppvar[1]._pval
#define _ion_ina	*_ppvar[2]._pval
#define _ion_dinadv	*_ppvar[3]._pval
#define _ion_ki	*_ppvar[4]._pval
#define _ion_ko	*_ppvar[5]._pval
#define _ion_ik	*_ppvar[6]._pval
#define _ion_dikdv	*_ppvar[7]._pval
 
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
 extern double celsius;
 /* declaration of user functions */
 static int _hoc_settables();
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
 "setdata_nakpump", _hoc_setdata,
 "settables_nakpump", _hoc_settables,
 0, 0
};
 /* declare global and static user variables */
#define F F_nakpump
 double F = 96.4853;
#define KdKi KdKi_nakpump
 double KdKi = 32.88;
#define KdNai KdNai_nakpump
 double KdNai = 4.05;
#define KdKo KdKo_nakpump
 double KdKo = 0.258;
#define KdNao KdNao_nakpump
 double KdNao = 69.8;
#define R R_nakpump
 double R = 8.31447;
#define k4 k4_nakpump
 double k4 = 0.165;
#define k3 k3_nakpump
 double k3 = 0.01;
#define k2 k2_nakpump
 double k2 = 0.04;
#define ouabain ouabain_nakpump
 double ouabain = 1;
#define y0 y0_nakpump
 double y0 = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "R_nakpump", "j/(K.molar",
 "F_nakpump", "C/mol",
 "ouabain_nakpump", "mM",
 "inakmax_nakpump", "mA/cm2",
 "Nai_inf_nakpump", "mM",
 "Ki_inf_nakpump", "mM",
 "ATPi_nakpump", "mM",
 "inak_nakpump", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "R_nakpump", &R_nakpump,
 "F_nakpump", &F_nakpump,
 "ouabain_nakpump", &ouabain_nakpump,
 "KdNao_nakpump", &KdNao_nakpump,
 "KdKo_nakpump", &KdKo_nakpump,
 "KdNai_nakpump", &KdNai_nakpump,
 "KdKi_nakpump", &KdKi_nakpump,
 "k2_nakpump", &k2_nakpump,
 "k3_nakpump", &k3_nakpump,
 "k4_nakpump", &k4_nakpump,
 "y0_nakpump", &y0_nakpump,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(), nrn_init(), nrn_state();
 static void nrn_cur(), nrn_jacob();
 
static int _ode_count(), _ode_map(), _ode_spec(), _ode_matsol();
 
#define _cvode_ieq _ppvar[8]._i
 /* connect range variables in _p that hoc is supposed to know about */
 static char *_mechanism[] = {
 "6.2.0",
"nakpump",
 "inakmax_nakpump",
 "Nai_inf_nakpump",
 "Ki_inf_nakpump",
 "ATPi_nakpump",
 0,
 "inak_nakpump",
 0,
 "y_nakpump",
 0,
 0};
 static Symbol* _na_sym;
 static Symbol* _k_sym;
 
static void nrn_alloc(_prop)
	Prop *_prop;
{
	Prop *prop_ion, *need_memb();
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 25, _prop);
 	/*initialize range parameters*/
 	inakmax = 1;
 	Nai_inf = 5;
 	Ki_inf = 14;
 	ATPi = 6.97;
 	_prop->param = _p;
 	_prop->param_size = 25;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 9, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* nai */
 	_ppvar[1]._pval = &prop_ion->param[2]; /* nao */
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[4]._pval = &prop_ion->param[1]; /* ki */
 	_ppvar[5]._pval = &prop_ion->param[2]; /* ko */
 	_ppvar[6]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[7]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 _nak_pump_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("na", -10000.);
 	ion_reg("k", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_dparam_size(_mechtype, 9);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 nakpump /Users/pbmanis/Desktop/Python/nrnlibrary/nrnlibrary/i386/nak_pump.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "NaKPump adapted from Kyoto's model";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
static settables();
 
static int _ode_spec1(), _ode_matsol1();
 static int _slist1[1], _dlist1[1];
 static int states();
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   settables ( _threadargscomma_ nai , ki , v , celsius , ATPi , ouabain ) ;
   Dy = alfay * ( 1.0 - y ) - betay * y ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 settables ( _threadargscomma_ nai , ki , v , celsius , ATPi , ouabain ) ;
 Dy = Dy  / (1. - dt*( (alfay)*(( ( - 1.0 ) )) - (betay)*(1.0) )) ;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   settables ( _threadargscomma_ nai , ki , v , celsius , ATPi , ouabain ) ;
    y = y + (1. - exp(dt*((alfay)*(( ( - 1.0 ) )) - (betay)*(1.0))))*(- ( (alfay)*(( 1.0 )) ) / ( (alfay)*(( ( - 1.0) )) - (betay)*(1.0) ) - y) ;
   }
  return 0;
}
 
static int  settables ( _p, _ppvar, _thread, _nt, _lnai , _lki , _lv , _lcelsius , _lATPi ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lnai , _lki , _lv , _lcelsius , _lATPi ;
 {
   drugblock = 1.0 / ( 1.0 + ouabain / 0.0006 ) ;
   k1 = 0.37 * ( 1.0 / ( 1.0 + 0.094 / _lATPi ) ) ;
   kelvin = _lcelsius + 273.16 ;
   Naeff = nao * exp ( - 0.82 * F * _lv / R / kelvin ) ;
   alfay = k2 * pE2Na + k4 * pE2K ;
   betay = k1 * pE1Na + k3 * pE1K ;
   pE1Na = 1.0 / ( 1.0 + pow( ( KdNai / _lnai ) , 1.06 ) * ( 1.0 + pow( ( _lki / KdKi ) , 1.12 ) ) ) ;
   pE1K = 1.0 / ( 1.0 + pow( ( KdKi / _lki ) , 1.12 ) * ( 1.0 + pow( ( _lnai / KdNai ) , 1.06 ) ) ) ;
   pE2Na = 1.0 / ( 1.0 + pow( ( KdNao / Naeff ) , 1.06 ) * ( 1.0 + pow( ( ko / KdKo ) , 1.12 ) ) ) ;
   pE2K = 1.0 / ( 1.0 + pow( ( KdKo / ko ) , 1.12 ) * ( 1.0 + pow( ( Naeff / KdNao ) , 1.06 ) ) ) ;
    return 0; }
 
static int _hoc_settables() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 settables ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) , *getarg(5) ) ;
 ret(_r);
}
 
static int _ode_count(_type) int _type;{ return 1;}
 
static int _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  nai = _ion_nai;
  nao = _ion_nao;
  ki = _ion_ki;
  ko = _ion_ko;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
   }}
 
static int _ode_map(_ieq, _pv, _pvdot, _pp, _ppd, _atol, _type) int _ieq, _type; double** _pv, **_pvdot, *_pp, *_atol; Datum* _ppd; { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 1; ++_i) {
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
  nai = _ion_nai;
  nao = _ion_nao;
  ki = _ion_ki;
  ko = _ion_ko;
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 2);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 3, 4);
   nrn_update_ion_pointer(_k_sym, _ppvar, 4, 1);
   nrn_update_ion_pointer(_k_sym, _ppvar, 5, 2);
   nrn_update_ion_pointer(_k_sym, _ppvar, 6, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 7, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  y = y0;
 {
   settables ( _threadargscomma_ Nai_inf , Ki_inf , v , celsius , ATPi , ouabain ) ;
   y = .6224 ;
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
  nai = _ion_nai;
  nao = _ion_nao;
  ki = _ion_ki;
  ko = _ion_ko;
 initmodel(_p, _ppvar, _thread, _nt);
  }}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   settables ( _threadargscomma_ nai , ki , v , celsius , ATPi , ouabain ) ;
   inak = inakmax * 0.001 * 11.5 * ( k1 * pE1Na * y - k2 * pE2Na * ( 1.0 - y ) ) * drugblock ;
   ina = 3.0 * inak ;
   ik = - 2.0 * inak ;
   }
 _current += ina;
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
  nai = _ion_nai;
  nao = _ion_nao;
  ki = _ion_ki;
  ko = _ion_ko;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dik;
 double _dina;
  _dina = ina;
  _dik = ik;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dinadv += (_dina - ina)/.001 ;
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
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
  nai = _ion_nai;
  nao = _ion_nao;
  ki = _ion_ki;
  ko = _ion_ko;
 { {
 for (; t < _break; t += dt) {
   states(_p, _ppvar, _thread, _nt);
  
}}
 t = _save;
 }  }}

}

static terminal(){}

static _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(y) - _p;  _dlist1[0] = &(Dy) - _p;
_first = 0;
}
