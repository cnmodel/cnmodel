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
#define gbar _p[0]
#define i _p[1]
#define h1 _p[2]
#define h2 _p[3]
#define thegna _p[4]
#define Dh1 _p[5]
#define Dh2 _p[6]
#define _g _p[7]
 
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
 static int _hoc_alp2();
 static int _hoc_alp1();
 static int _hoc_bet2();
 static int _hoc_bet1();
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
 "setdata_hcno", _hoc_setdata,
 "alp2_hcno", _hoc_alp2,
 "alp1_hcno", _hoc_alp1,
 "bet2_hcno", _hoc_bet2,
 "bet1_hcno", _hoc_bet1,
 "trates_hcno", _hoc_trates,
 0, 0
};
#define alp2 alp2_hcno
#define alp1 alp1_hcno
#define bet2 bet2_hcno
#define bet1 bet1_hcno
 extern double alp2();
 extern double alp1();
 extern double bet2();
 extern double bet1();
 /* declare global and static user variables */
#define a02 a02_hcno
 double a02 = 0.0029;
#define a01 a01_hcno
 double a01 = 0.008;
#define eh eh_hcno
 double eh = 0;
#define frac frac_hcno
 double frac = 0;
#define gm2 gm2_hcno
 double gm2 = 0.6;
#define gm1 gm1_hcno
 double gm1 = 0.3;
#define hinf hinf_hcno
 double hinf = 0;
#define q10 q10_hcno
 double q10 = 4.5;
#define qinf qinf_hcno
 double qinf = 7;
#define thinf thinf_hcno
 double thinf = -66;
#define tau2 tau2_hcno
 double tau2 = 0;
#define tau1 tau1_hcno
 double tau1 = 0;
#define vhalf2 vhalf2_hcno
 double vhalf2 = -84;
#define vhalf1 vhalf1_hcno
 double vhalf1 = -50;
#define zeta2 zeta2_hcno
 double zeta2 = 3;
#define zeta1 zeta1_hcno
 double zeta1 = 3;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "vhalf1_hcno", "mV",
 "vhalf2_hcno", "mV",
 "gm1_hcno", "mV",
 "gm2_hcno", "mV",
 "zeta1_hcno", "/ms",
 "zeta2_hcno", "/ms",
 "thinf_hcno", "mV",
 "qinf_hcno", "mV",
 "eh_hcno", "mV",
 "gbar_hcno", "mho/cm2",
 "i_hcno", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double h20 = 0;
 static double h10 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "vhalf1_hcno", &vhalf1_hcno,
 "vhalf2_hcno", &vhalf2_hcno,
 "gm1_hcno", &gm1_hcno,
 "gm2_hcno", &gm2_hcno,
 "zeta1_hcno", &zeta1_hcno,
 "zeta2_hcno", &zeta2_hcno,
 "a01_hcno", &a01_hcno,
 "a02_hcno", &a02_hcno,
 "frac_hcno", &frac_hcno,
 "thinf_hcno", &thinf_hcno,
 "qinf_hcno", &qinf_hcno,
 "q10_hcno", &q10_hcno,
 "eh_hcno", &eh_hcno,
 "hinf_hcno", &hinf_hcno,
 "tau1_hcno", &tau1_hcno,
 "tau2_hcno", &tau2_hcno,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(), nrn_init(), nrn_state();
 static void nrn_cur(), nrn_jacob();
 
static int _ode_count(), _ode_map(), _ode_spec(), _ode_matsol();
 
#define _cvode_ieq _ppvar[0]._i
 /* connect range variables in _p that hoc is supposed to know about */
 static char *_mechanism[] = {
 "6.2.0",
"hcno",
 "gbar_hcno",
 0,
 "i_hcno",
 0,
 "h1_hcno",
 "h2_hcno",
 0,
 0};
 
static void nrn_alloc(_prop)
	Prop *_prop;
{
	Prop *prop_ion, *need_memb();
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 8, _prop);
 	/*initialize range parameters*/
 	gbar = 0.0005;
 	_prop->param = _p;
 	_prop->param_size = 8;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 1, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 _hcno_reg() {
	int _vectorized = 0;
  _initlists();
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
  hoc_register_dparam_size(_mechtype, 1);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 hcno /Users/pbmanis/Desktop/Python/nrnlibrary/nrnlibrary/i386/hcno.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "h current for Octopus cells of Cochlear Nucleus";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
static trates();
 static int _deriv1_advance = 0;
 
static int _ode_spec1(), _ode_matsol1();
 static int _slist2[2]; static double _dlist2[2];
 static double _savstate1[2], *_temp1 = _savstate1;
 static int _slist1[2], _dlist1[2];
 static int states();
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   trates ( _threadargscomma_ v ) ;
   Dh1 = ( hinf - h1 ) / tau1 ;
   Dh2 = ( hinf - h2 ) / tau2 ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 trates ( _threadargscomma_ v ) ;
 Dh1 = Dh1  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau1 )) ;
 Dh2 = Dh2  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau2 )) ;
}
 /*END CVODE*/
 
static int states () {_reset=0;
 { static int _recurse = 0;
 int _counte = -1;
 if (!_recurse) {
 _recurse = 1;
 {int _id; for(_id=0; _id < 2; _id++) { _savstate1[_id] = _p[_slist1[_id]];}}
 error = newton(2,_slist2, _p, states, _dlist2);
 _recurse = 0; if(error) {abort_run(error);}}
 {
   trates ( _threadargscomma_ v ) ;
   Dh1 = ( hinf - h1 ) / tau1 ;
   Dh2 = ( hinf - h2 ) / tau2 ;
   {int _id; for(_id=0; _id < 2; _id++) {
if (_deriv1_advance) {
 _dlist2[++_counte] = _p[_dlist1[_id]] - (_p[_slist1[_id]] - _savstate1[_id])/dt;
 }else{
_dlist2[++_counte] = _p[_slist1[_id]] - _savstate1[_id];}}}
 } }
 return _reset;}
 
static int  trates (  _lv )  
	double _lv ;
 {
   double _lqt ;
 _lqt = pow( q10 , ( ( celsius - 33.0 ) / 10.0 ) ) ;
   tau1 = bet1 ( _threadargscomma_ _lv ) / ( _lqt * a01 * ( 1.0 + alp1 ( _threadargscomma_ _lv ) ) ) ;
   tau2 = bet2 ( _threadargscomma_ _lv ) / ( _lqt * a02 * ( 1.0 + alp2 ( _threadargscomma_ _lv ) ) ) ;
   hinf = 1.0 / ( 1.0 + exp ( ( _lv - thinf ) / qinf ) ) ;
    return 0; }
 
static int _hoc_trates() {
  double _r;
   _r = 1.;
 trates (  *getarg(1) ) ;
 ret(_r);
}
 
double alp1 (  _lv )  
	double _lv ;
 {
   double _lalp1;
 _lalp1 = exp ( 1.e-3 * zeta1 * ( _lv - vhalf1 ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lalp1;
 }
 
static int _hoc_alp1() {
  double _r;
   _r =  alp1 (  *getarg(1) ) ;
 ret(_r);
}
 
double bet1 (  _lv )  
	double _lv ;
 {
   double _lbet1;
 _lbet1 = exp ( 1.e-3 * zeta1 * gm1 * ( _lv - vhalf1 ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lbet1;
 }
 
static int _hoc_bet1() {
  double _r;
   _r =  bet1 (  *getarg(1) ) ;
 ret(_r);
}
 
double alp2 (  _lv )  
	double _lv ;
 {
   double _lalp2;
 _lalp2 = exp ( 1.e-3 * zeta2 * ( _lv - vhalf2 ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lalp2;
 }
 
static int _hoc_alp2() {
  double _r;
   _r =  alp2 (  *getarg(1) ) ;
 ret(_r);
}
 
double bet2 (  _lv )  
	double _lv ;
 {
   double _lbet2;
 _lbet2 = exp ( 1.e-3 * zeta2 * gm2 * ( _lv - vhalf2 ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lbet2;
 }
 
static int _hoc_bet2() {
  double _r;
   _r =  bet2 (  *getarg(1) ) ;
 ret(_r);
}
 
static int _ode_count(_type) int _type;{ return 2;}
 
static int _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 ();
 }}
 
static int _ode_map(_ieq, _pv, _pvdot, _pp, _ppd, _atol, _type) int _ieq, _type; double** _pv, **_pvdot, *_pp, *_atol; Datum* _ppd; { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
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
 _ode_matsol1 ();
 }}

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h2 = h20;
  h1 = h10;
 {
   trates ( _threadargscomma_ v ) ;
   h1 = hinf ;
   h2 = hinf ;
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
   thegna = gbar * ( h1 * frac + h2 * ( 1.0 - frac ) ) ;
   i = thegna * ( v - eh ) ;
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
 error = 0; _deriv1_advance = 1;
 derivimplicit(_ninits, 2, _slist1, _dlist1, _p, &t, dt, states, &_temp1);
_deriv1_advance = 0;
 if(error){fprintf(stderr,"at line 54 in file hcno.mod:\n        SOLVE states METHOD derivimplicit\n"); nrn_complain(_p); abort_run(error);}
 
}}
 t = _save;
 }}}

}

static terminal(){}

static _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(h1) - _p;  _dlist1[0] = &(Dh1) - _p;
 _slist1[1] = &(h2) - _p;  _dlist1[1] = &(Dh2) - _p;
 _slist2[0] = &(h2) - _p;
 _slist2[1] = &(h1) - _p;
_first = 0;
}
