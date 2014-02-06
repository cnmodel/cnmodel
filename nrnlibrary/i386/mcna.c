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
#define lp _p[0]
#define ml _p[1]
#define nm _p[2]
#define porate _p[3]
#define ina _p[4]
#define gna _p[5]
#define P _p[6]
#define L _p[7]
#define M _p[8]
#define N _p[9]
#define O _p[10]
#define ena _p[11]
#define DP _p[12]
#define DL _p[13]
#define DM _p[14]
#define DN _p[15]
#define DO _p[16]
#define v _p[17]
#define _g _p[18]
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
 
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
 static int _hoc_alp();
 static int _hoc_bet();
 static int _hoc_expM1();
 static int _hoc_rate();
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
 "setdata_MCna", _hoc_setdata,
 "alp_MCna", _hoc_alp,
 "bet_MCna", _hoc_bet,
 "expM1_MCna", _hoc_expM1,
 "rate_MCna", _hoc_rate,
 0, 0
};
#define alp alp_MCna
#define bet bet_MCna
#define expM1 expM1_MCna
 extern double alp();
 extern double bet();
 extern double expM1();
 
static void _check_rate(double*, Datum*, Datum*, _NrnThread*); 
static void _check_table_thread(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, int _type) {
   _check_rate(_p, _ppvar, _thread, _nt);
 }
 #define _zam _thread[2]._pval[0]
 #define _zbm _thread[2]._pval[1]
 #define _zah _thread[2]._pval[2]
 #define _zbh _thread[2]._pval[3]
 /* declare global and static user variables */
#define gnabar gnabar_MCna
 double gnabar = 0.12;
#define usetable usetable_MCna
 double usetable = 1;
#define vrest vrest_MCna
 double vrest = -55;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_MCna", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gnabar_MCna", "mho/cm2",
 "vrest_MCna", "mV",
 "ina_MCna", "mA/cm2",
 "gna_MCna", "mho/cm2",
 0,0
};
 static double L0 = 0;
 static double M0 = 0;
 static double N0 = 0;
 static double O0 = 0;
 static double P0 = 0;
 static double delta_t = 1;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "gnabar_MCna", &gnabar_MCna,
 "vrest_MCna", &vrest_MCna,
 "usetable_MCna", &usetable_MCna,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(), nrn_init(), nrn_state();
 static void nrn_cur(), nrn_jacob();
 
static int _ode_count(), _ode_map(), _ode_spec(), _ode_matsol();
 
#define _cvode_ieq _ppvar[3]._i
 /* connect range variables in _p that hoc is supposed to know about */
 static char *_mechanism[] = {
 "6.2.0",
"MCna",
 "lp_MCna",
 "ml_MCna",
 "nm_MCna",
 "porate_MCna",
 0,
 "ina_MCna",
 "gna_MCna",
 0,
 "P_MCna",
 "L_MCna",
 "M_MCna",
 "N_MCna",
 "O_MCna",
 0,
 0};
 static Symbol* _na_sym;
 
static void nrn_alloc(_prop)
	Prop *_prop;
{
	Prop *prop_ion, *need_memb();
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 19, _prop);
 	/*initialize range parameters*/
 	lp = 1.9;
 	ml = 0.75;
 	nm = 0.3;
 	porate = 1;
 	_prop->param = _p;
 	_prop->param_size = 19;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 
}
 static _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _thread_mem_init(Datum*);
 static void _thread_cleanup(Datum*);
 static void _update_ion_pointer(Datum*);
 _mcna_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("na", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 4);
  _extcall_thread = (Datum*)ecalloc(3, sizeof(Datum));
  _thread_mem_init(_extcall_thread);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_thread_reg(_mechtype, 1, _thread_mem_init);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
     _nrn_thread_table_reg(_mechtype, _check_table_thread);
  hoc_register_dparam_size(_mechtype, 4);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 MCna /Users/pbmanis/Desktop/Python/nrnlibrary/nrnlibrary/i386/mcna.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 /*Top LOCAL _zam , _zbm , _zah , _zbh */
 static double *_t__zam;
 static double *_t__zah;
 static double *_t__zbm;
 static double *_t__zbh;
static int _reset;
static char *modelname = "Moore-Cox sodium channel";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
static _f_rate();
static rate();
 
#define _MATELM1(_row,_col) *(_nrn_thread_getelm(_so, _row + 1, _col + 1))
 
#define _RHS1(_arg) _rhs[_arg+1]
  static int _cvspth1 = 1;
 
static int _ode_spec1(), _ode_matsol1();
 extern double *_nrn_thread_getelm();
 
#define _MATELM1(_row,_col) *(_nrn_thread_getelm(_so, _row + 1, _col + 1))
 
#define _RHS1(_arg) _rhs[_arg+1]
  
#define _linmat1  1
 static int _spth1 = 0;
 static _n_rate();
 static int _slist1[5], _dlist1[5]; static double *_temp1;
 static int states();
 
static int states (void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt)
 {int _reset=0;
 {
   double b_flux, f_flux, _term; int _i;
 {int _i; double _dt1 = 1.0/dt;
for(_i=1;_i<5;_i++){
  	_RHS1(_i) = -_dt1*(_p[_slist1[_i]] - _p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 rate ( _threadargscomma_ v * 1.0 ) ;
    /* P + L + M + N + O = 1.0 */
 _RHS1(0) =  1.0;
 _MATELM1(0, 0) = 1;
 _RHS1(0) -= O ;
 _MATELM1(0, 3) = 1;
 _RHS1(0) -= N ;
 _MATELM1(0, 2) = 1;
 _RHS1(0) -= M ;
 _MATELM1(0, 1) = 1;
 _RHS1(0) -= L ;
 _MATELM1(0, 4) = 1;
 _RHS1(0) -= P ;
 /*CONSERVATION*/
 /* ~ P <-> L ( _zam , lp * _zbm )*/
 f_flux =  _zam * P ;
 b_flux =  lp * _zbm * L ;
 _RHS1( 4) -= (f_flux - b_flux);
 _RHS1( 1) += (f_flux - b_flux);
 
 _term =  _zam ;
 _MATELM1( 4 ,4)  += _term;
 _MATELM1( 1 ,4)  -= _term;
 _term =  lp * _zbm ;
 _MATELM1( 4 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ L <-> M ( 2.0 * _zam , ml * _zbm )*/
 f_flux =  2.0 * _zam * L ;
 b_flux =  ml * _zbm * M ;
 _RHS1( 1) -= (f_flux - b_flux);
 _RHS1( 2) += (f_flux - b_flux);
 
 _term =  2.0 * _zam ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 2 ,1)  -= _term;
 _term =  ml * _zbm ;
 _MATELM1( 1 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  /* ~ M <-> N ( 3.0 * _zam , nm * _zbm )*/
 f_flux =  3.0 * _zam * M ;
 b_flux =  nm * _zbm * N ;
 _RHS1( 2) -= (f_flux - b_flux);
 _RHS1( 3) += (f_flux - b_flux);
 
 _term =  3.0 * _zam ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 3 ,2)  -= _term;
 _term =  nm * _zbm ;
 _MATELM1( 2 ,3)  -= _term;
 _MATELM1( 3 ,3)  += _term;
 /*REACTION*/
  /* ~ N <-> O ( 1.1 * _zbh , 0.0 )*/
 f_flux =  1.1 * _zbh * N ;
 b_flux =  0.0 * O ;
 _RHS1( 3) -= (f_flux - b_flux);
 
 _term =  1.1 * _zbh ;
 _MATELM1( 3 ,3)  += _term;
 _term =  0.0 ;
 _MATELM1( 3 ,0)  -= _term;
 /*REACTION*/
  /* ~ N <-> P ( 3.0 * _zbm , 0.0 )*/
 f_flux =  3.0 * _zbm * N ;
 b_flux =  0.0 * P ;
 _RHS1( 3) -= (f_flux - b_flux);
 _RHS1( 4) += (f_flux - b_flux);
 
 _term =  3.0 * _zbm ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 4 ,3)  -= _term;
 _term =  0.0 ;
 _MATELM1( 3 ,4)  -= _term;
 _MATELM1( 4 ,4)  += _term;
 /*REACTION*/
  /* ~ P <-> O ( _zbh , _zah )*/
 f_flux =  _zbh * P ;
 b_flux =  _zah * O ;
 _RHS1( 4) -= (f_flux - b_flux);
 
 _term =  _zbh ;
 _MATELM1( 4 ,4)  += _term;
 _term =  _zah ;
 _MATELM1( 4 ,0)  -= _term;
 /*REACTION*/
    } return _reset;
 }
 
double alp ( _p, _ppvar, _thread, _nt, _lv , _li ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lv , _li ;
 {
   double _lalp;
 double _la , _lb , _lc , _lq10 ;
 _lv = - _lv + vrest ;
   _lq10 = pow( 3.0 , ( ( celsius - 6.3 ) / 10.0 ) ) ;
   if ( _li  == 0.0 ) {
     _lalp = _lq10 * 0.1 * expM1 ( _threadargscomma_ _lv + 25.0 , 10.0 ) ;
     }
   else if ( _li  == 1.0 ) {
     _lalp = _lq10 * 0.07 * exp ( _lv / 20.0 ) ;
     }
   
return _lalp;
 }
 
static int _hoc_alp() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  alp ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) ) ;
 ret(_r);
}
 
double bet ( _p, _ppvar, _thread, _nt, _lv , _li ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lv , _li ;
 {
   double _lbet;
 double _la , _lb , _lc , _lq10 ;
 _lv = - _lv + vrest ;
   _lq10 = pow( 3.0 , ( ( celsius - 6.3 ) / 10.0 ) ) ;
   if ( _li  == 0.0 ) {
     _lbet = _lq10 * 4.0 * exp ( _lv / 18.0 ) ;
     }
   else if ( _li  == 1.0 ) {
     _lbet = _lq10 * 1.0 / ( exp ( 0.1 * _lv + 3.0 ) + 1.0 ) ;
     }
   
return _lbet;
 }
 
static int _hoc_bet() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  bet ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) ) ;
 ret(_r);
}
 
double expM1 ( _p, _ppvar, _thread, _nt, _lx , _ly ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx , _ly ;
 {
   double _lexpM1;
 if ( fabs ( _lx / _ly ) < 1e-6 ) {
     _lexpM1 = _ly * ( 1.0 - _lx / _ly / 2.0 ) ;
     }
   else {
     _lexpM1 = _lx / ( exp ( _lx / _ly ) - 1.0 ) ;
     }
   
return _lexpM1;
 }
 
static int _hoc_expM1() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  expM1 ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) ) ;
 ret(_r);
}
 static double _mfac_rate, _tmin_rate;
  static void _check_rate(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_dt;
  static double _sav_celsius;
  if (!usetable) {return;}
  if (_sav_dt != dt) { _maktable = 1;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_rate =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_rate)/200.; _mfac_rate = 1./_dx;
   for (_i=0, _x=_tmin_rate; _i < 201; _x += _dx, _i++) {
    _f_rate(_p, _ppvar, _thread, _nt, _x);
    _t__zam[_i] = _zam;
    _t__zah[_i] = _zah;
    _t__zbm[_i] = _zbm;
    _t__zbh[_i] = _zbh;
   }
   _sav_dt = dt;
   _sav_celsius = celsius;
  }
 }

 static rate(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv) { 
#if 0
_check_rate(_p, _ppvar, _thread, _nt);
#endif
 _n_rate(_p, _ppvar, _thread, _nt, _lv);
 return;
 }

 static _n_rate(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_rate(_p, _ppvar, _thread, _nt, _lv); return; 
}
 _xi = _mfac_rate * (_lv - _tmin_rate);
 _i = (int) _xi;
 if (_xi <= 0.) {
 _zam = _t__zam[0];
 _zah = _t__zah[0];
 _zbm = _t__zbm[0];
 _zbh = _t__zbh[0];
 return; }
 if (_i >= 200) {
 _zam = _t__zam[200];
 _zah = _t__zah[200];
 _zbm = _t__zbm[200];
 _zbh = _t__zbh[200];
 return; }
 _theta = _xi - (double)_i;
 _zam = _t__zam[_i] + _theta*(_t__zam[_i+1] - _t__zam[_i]);
 _zah = _t__zah[_i] + _theta*(_t__zah[_i+1] - _t__zah[_i]);
 _zbm = _t__zbm[_i] + _theta*(_t__zbm[_i+1] - _t__zbm[_i]);
 _zbh = _t__zbh[_i] + _theta*(_t__zbh[_i+1] - _t__zbh[_i]);
 }

 
static int  _f_rate ( _p, _ppvar, _thread, _nt, _lv ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lv ;
 {
   double _la , _lb , _ltau ;
 _zam = alp ( _threadargscomma_ _lv , 0.0 ) ;
   _zah = alp ( _threadargscomma_ _lv , 1.0 ) ;
   _zbm = bet ( _threadargscomma_ _lv , 0.0 ) ;
   _zbh = bet ( _threadargscomma_ _lv , 1.0 ) ;
    return 0; }
 
static int _hoc_rate() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 
#if 1
 _check_rate(_p, _ppvar, _thread, _nt);
#endif
 _r = 1.;
 rate ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
/*CVODE ode begin*/
 static int _ode_spec1(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
 {int _i; for(_i=0;_i<5;_i++) _p[_dlist1[_i]] = 0.0;}
 rate ( _threadargscomma_ v * 1.0 ) ;
  /* P + L + M + N + O = 1.0 */
 /*CONSERVATION*/
 /* ~ P <-> L ( _zam , lp * _zbm )*/
 f_flux =  _zam * P ;
 b_flux =  lp * _zbm * L ;
 DP -= (f_flux - b_flux);
 DL += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ L <-> M ( 2.0 * _zam , ml * _zbm )*/
 f_flux =  2.0 * _zam * L ;
 b_flux =  ml * _zbm * M ;
 DL -= (f_flux - b_flux);
 DM += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ M <-> N ( 3.0 * _zam , nm * _zbm )*/
 f_flux =  3.0 * _zam * M ;
 b_flux =  nm * _zbm * N ;
 DM -= (f_flux - b_flux);
 DN += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ N <-> O ( 1.1 * _zbh , 0.0 )*/
 f_flux =  1.1 * _zbh * N ;
 b_flux =  0.0 * O ;
 DN -= (f_flux - b_flux);
 DO += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ N <-> P ( 3.0 * _zbm , 0.0 )*/
 f_flux =  3.0 * _zbm * N ;
 b_flux =  0.0 * P ;
 DN -= (f_flux - b_flux);
 DP += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ P <-> O ( _zbh , _zah )*/
 f_flux =  _zbh * P ;
 b_flux =  _zah * O ;
 DP -= (f_flux - b_flux);
 DO += (f_flux - b_flux);
 
 /*REACTION*/
    } return _reset;
 }
 
/*CVODE matsol*/
 static int _ode_matsol1(void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
   b_flux = f_flux = 0.;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<5;_i++){
  	_RHS1(_i) = _dt1*(_p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 rate ( _threadargscomma_ v * 1.0 ) ;
  /* P + L + M + N + O = 1.0 */
 /*CONSERVATION*/
 /* ~ P <-> L ( _zam , lp * _zbm )*/
 _term =  _zam ;
 _MATELM1( 4 ,4)  += _term;
 _MATELM1( 1 ,4)  -= _term;
 _term =  lp * _zbm ;
 _MATELM1( 4 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ L <-> M ( 2.0 * _zam , ml * _zbm )*/
 _term =  2.0 * _zam ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 2 ,1)  -= _term;
 _term =  ml * _zbm ;
 _MATELM1( 1 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  /* ~ M <-> N ( 3.0 * _zam , nm * _zbm )*/
 _term =  3.0 * _zam ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 3 ,2)  -= _term;
 _term =  nm * _zbm ;
 _MATELM1( 2 ,3)  -= _term;
 _MATELM1( 3 ,3)  += _term;
 /*REACTION*/
  /* ~ N <-> O ( 1.1 * _zbh , 0.0 )*/
 _term =  1.1 * _zbh ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 0 ,3)  -= _term;
 _term =  0.0 ;
 _MATELM1( 3 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
  /* ~ N <-> P ( 3.0 * _zbm , 0.0 )*/
 _term =  3.0 * _zbm ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 4 ,3)  -= _term;
 _term =  0.0 ;
 _MATELM1( 3 ,4)  -= _term;
 _MATELM1( 4 ,4)  += _term;
 /*REACTION*/
  /* ~ P <-> O ( _zbh , _zah )*/
 _term =  _zbh ;
 _MATELM1( 4 ,4)  += _term;
 _MATELM1( 0 ,4)  -= _term;
 _term =  _zah ;
 _MATELM1( 4 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
    } return _reset;
 }
 
/*CVODE end*/
 
static int _ode_count(_type) int _type;{ return 5;}
 
static int _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static int _ode_map(_ieq, _pv, _pvdot, _pp, _ppd, _atol, _type) int _ieq, _type; double** _pv, **_pvdot, *_pp, *_atol; Datum* _ppd; { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 5; ++_i) {
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
  ena = _ion_ena;
 _cvode_sparse_thread(&_thread[_cvspth1]._pvoid, 5, _dlist1, _p, _ode_matsol1, _ppvar, _thread, _nt);
 }}
 
static void _thread_mem_init(Datum* _thread) {
   _thread[2]._pval = (double*)ecalloc(4, sizeof(double));
 }
 
static void _thread_cleanup(Datum* _thread) {
   _nrn_destroy_sparseobj_thread(_thread[_spth1]._pvoid);
   _nrn_destroy_sparseobj_thread(_thread[_cvspth1]._pvoid);
   free((void*)(_thread[2]._pval));
 }
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  L = L0;
  M = M0;
  N = N0;
  O = O0;
  P = P0;
 {
   P = 1.0 ;
   rate ( _threadargscomma_ v * 1.0 ) ;
    _ss_sparse_thread(&_thread[_spth1]._pvoid, 5, _slist1, _dlist1, _p, &t, dt, states, _linmat1, _ppvar, _thread, _nt);
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

#if 0
 _check_rate(_p, _ppvar, _thread, _nt);
#endif
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
  ena = _ion_ena;
 initmodel(_p, _ppvar, _thread, _nt);
 }}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   ina = gnabar * N * ( v - ena ) ;
   }
 _current += ina;

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
  ena = _ion_ena;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
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
  ena = _ion_ena;
 { {
 for (; t < _break; t += dt) {
  sparse_thread(&_thread[_spth1]._pvoid, 5, _slist1, _dlist1, _p, &t, dt, states, _linmat1, _ppvar, _thread, _nt);
  
}}
 t = _save;
 } }}

}

static terminal(){}

static _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
   _t__zam = makevector(201*sizeof(double));
   _t__zah = makevector(201*sizeof(double));
   _t__zbm = makevector(201*sizeof(double));
   _t__zbh = makevector(201*sizeof(double));
 _slist1[0] = &(O) - _p;  _dlist1[0] = &(DO) - _p;
 _slist1[1] = &(L) - _p;  _dlist1[1] = &(DL) - _p;
 _slist1[2] = &(M) - _p;  _dlist1[2] = &(DM) - _p;
 _slist1[3] = &(N) - _p;  _dlist1[3] = &(DN) - _p;
 _slist1[4] = &(P) - _p;  _dlist1[4] = &(DP) - _p;
_first = 0;
}
