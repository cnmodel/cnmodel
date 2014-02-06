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
#define Erev _p[0]
#define gmax _p[1]
#define k1 _p[2]
#define km1 _p[3]
#define a1 _p[4]
#define b1 _p[5]
#define r1 _p[6]
#define d1 _p[7]
#define r2 _p[8]
#define d2 _p[9]
#define r3 _p[10]
#define d3 _p[11]
#define rd _p[12]
#define dd _p[13]
#define i _p[14]
#define g _p[15]
#define f1 _p[16]
#define f2 _p[17]
#define Open _p[18]
#define C0 _p[19]
#define C1 _p[20]
#define C2 _p[21]
#define D1 _p[22]
#define O1 _p[23]
#define D2 _p[24]
#define D3 _p[25]
#define DC0 _p[26]
#define DC1 _p[27]
#define DC2 _p[28]
#define DD1 _p[29]
#define DO1 _p[30]
#define DD2 _p[31]
#define DD3 _p[32]
#define _g _p[33]
#define _nd_area  *_ppvar[0]._pval
#define XMTR	*_ppvar[2]._pval
#define _p_XMTR	_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 static int hoc_nrnpointerindex =  2;
 /* external NEURON variables */
 /* declaration of user functions */
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
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "Erev", "mV",
 "gmax", "pS",
 "k1", "/uM",
 "km1", "/ms",
 "a1", "/ms",
 "b1", "/ms",
 "r1", "/ms",
 "d1", "/ms",
 "r2", "/ms",
 "d2", "/ms",
 "r3", "/ms",
 "d3", "/ms",
 "rd", "/ms",
 "dd", "/ms",
 "i", "nA",
 "g", "pS",
 "f1", "/ms",
 "f2", "/ms",
 "Open", "1",
 "XMTR", "mM",
 0,0
};
 static double C20 = 0;
 static double C10 = 0;
 static double C00 = 0;
 static double D30 = 0;
 static double D20 = 0;
 static double D10 = 0;
 static double O10 = 0;
 static double delta_t = 1;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(), nrn_init(), nrn_state();
 static void nrn_cur(), nrn_jacob();
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 
static int _ode_count(), _ode_map(), _ode_spec(), _ode_matsol();
 
#define _cvode_ieq _ppvar[3]._i
 /* connect range variables in _p that hoc is supposed to know about */
 static char *_mechanism[] = {
 "6.2.0",
"GLYaGC",
 "Erev",
 "gmax",
 "k1",
 "km1",
 "a1",
 "b1",
 "r1",
 "d1",
 "r2",
 "d2",
 "r3",
 "d3",
 "rd",
 "dd",
 0,
 "i",
 "g",
 "f1",
 "f2",
 "Open",
 0,
 "C0",
 "C1",
 "C2",
 "D1",
 "O1",
 "D2",
 "D3",
 0,
 "XMTR",
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
 	_p = nrn_prop_data_alloc(_mechtype, 34, _prop);
 	/*initialize range parameters*/
 	Erev = -70;
 	gmax = 500;
 	k1 = 12.81;
 	km1 = 0.0087;
 	a1 = 0.0194;
 	b1 = 1.138;
 	r1 = 5.19;
 	d1 = 0.000462;
 	r2 = 0.731;
 	d2 = 1.641;
 	r3 = 3.817;
 	d3 = 1.806;
 	rd = 1;
 	dd = 1;
  }
 	_prop->param = _p;
 	_prop->param_size = 34;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 _Gly5GC_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func,
	 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
  hoc_register_dparam_size(_mechtype, 4);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 GLYaGC /Users/pbmanis/Desktop/Python/nrnlibrary/nrnlibrary/i386/Gly5GC.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "detailed model of Glycine receptors";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
 extern double *_getelm();
 
#define _MATELM1(_row,_col)	*(_getelm(_row + 1, _col + 1))
 
#define _RHS1(_arg) _coef1[_arg + 1]
 static double *_coef1;
 
#define _linmat1  1
 static void* _sparseobj1;
 static void* _cvsparseobj1;
 
static int _ode_spec1(), _ode_matsol1();
 static int _slist1[7], _dlist1[7]; static double *_temp1;
 static int kstates();
 
static int kstates ()
 {_reset=0;
 {
   double b_flux, f_flux, _term; int _i;
 {int _i; double _dt1 = 1.0/dt;
for(_i=1;_i<7;_i++){
  	_RHS1(_i) = -_dt1*(_p[_slist1[_i]] - _p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 f1 = k1 * ( 1e3 ) * XMTR ;
   f2 = k1 * ( 1e3 ) * XMTR ;
   /* ~ C0 <-> C1 ( f1 , km1 )*/
 f_flux =  f1 * C0 ;
 b_flux =  km1 * C1 ;
 _RHS1( 3) -= (f_flux - b_flux);
 _RHS1( 2) += (f_flux - b_flux);
 
 _term =  f1 ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 2 ,3)  -= _term;
 _term =  km1 ;
 _MATELM1( 3 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  /* ~ C1 <-> C2 ( f2 , 12.5 * km1 )*/
 f_flux =  f2 * C1 ;
 b_flux =  12.5 * km1 * C2 ;
 _RHS1( 2) -= (f_flux - b_flux);
 _RHS1( 1) += (f_flux - b_flux);
 
 _term =  f2 ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 1 ,2)  -= _term;
 _term =  12.5 * km1 ;
 _MATELM1( 2 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ C2 <-> O1 ( a1 , b1 )*/
 f_flux =  a1 * C2 ;
 b_flux =  b1 * O1 ;
 _RHS1( 1) -= (f_flux - b_flux);
 
 _term =  a1 ;
 _MATELM1( 1 ,1)  += _term;
 _term =  b1 ;
 _MATELM1( 1 ,0)  -= _term;
 /*REACTION*/
  /* ~ C1 <-> D1 ( r1 , d1 )*/
 f_flux =  r1 * C1 ;
 b_flux =  d1 * D1 ;
 _RHS1( 2) -= (f_flux - b_flux);
 _RHS1( 6) += (f_flux - b_flux);
 
 _term =  r1 ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 6 ,2)  -= _term;
 _term =  d1 ;
 _MATELM1( 2 ,6)  -= _term;
 _MATELM1( 6 ,6)  += _term;
 /*REACTION*/
  /* ~ C2 <-> D2 ( r2 , d2 )*/
 f_flux =  r2 * C2 ;
 b_flux =  d2 * D2 ;
 _RHS1( 1) -= (f_flux - b_flux);
 _RHS1( 5) += (f_flux - b_flux);
 
 _term =  r2 ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 5 ,1)  -= _term;
 _term =  d2 ;
 _MATELM1( 1 ,5)  -= _term;
 _MATELM1( 5 ,5)  += _term;
 /*REACTION*/
  /* ~ D1 <-> D2 ( rd , dd )*/
 f_flux =  rd * D1 ;
 b_flux =  dd * D2 ;
 _RHS1( 6) -= (f_flux - b_flux);
 _RHS1( 5) += (f_flux - b_flux);
 
 _term =  rd ;
 _MATELM1( 6 ,6)  += _term;
 _MATELM1( 5 ,6)  -= _term;
 _term =  dd ;
 _MATELM1( 6 ,5)  -= _term;
 _MATELM1( 5 ,5)  += _term;
 /*REACTION*/
  /* ~ D2 <-> D3 ( r3 , d3 )*/
 f_flux =  r3 * D2 ;
 b_flux =  d3 * D3 ;
 _RHS1( 5) -= (f_flux - b_flux);
 _RHS1( 4) += (f_flux - b_flux);
 
 _term =  r3 ;
 _MATELM1( 5 ,5)  += _term;
 _MATELM1( 4 ,5)  -= _term;
 _term =  d3 ;
 _MATELM1( 5 ,4)  -= _term;
 _MATELM1( 4 ,4)  += _term;
 /*REACTION*/
   /* C0 + C1 + C2 + D1 + D2 + D3 + O1 = 1.0 */
 _RHS1(0) =  1.0;
 _MATELM1(0, 0) = 1;
 _RHS1(0) -= O1 ;
 _MATELM1(0, 4) = 1;
 _RHS1(0) -= D3 ;
 _MATELM1(0, 5) = 1;
 _RHS1(0) -= D2 ;
 _MATELM1(0, 6) = 1;
 _RHS1(0) -= D1 ;
 _MATELM1(0, 1) = 1;
 _RHS1(0) -= C2 ;
 _MATELM1(0, 2) = 1;
 _RHS1(0) -= C1 ;
 _MATELM1(0, 3) = 1;
 _RHS1(0) -= C0 ;
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE ode begin*/
 static int _ode_spec1() {_reset=0;{
 double b_flux, f_flux, _term; int _i;
 {int _i; for(_i=0;_i<7;_i++) _p[_dlist1[_i]] = 0.0;}
 f1 = k1 * ( 1e3 ) * XMTR ;
 f2 = k1 * ( 1e3 ) * XMTR ;
 /* ~ C0 <-> C1 ( f1 , km1 )*/
 f_flux =  f1 * C0 ;
 b_flux =  km1 * C1 ;
 DC0 -= (f_flux - b_flux);
 DC1 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ C1 <-> C2 ( f2 , 12.5 * km1 )*/
 f_flux =  f2 * C1 ;
 b_flux =  12.5 * km1 * C2 ;
 DC1 -= (f_flux - b_flux);
 DC2 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ C2 <-> O1 ( a1 , b1 )*/
 f_flux =  a1 * C2 ;
 b_flux =  b1 * O1 ;
 DC2 -= (f_flux - b_flux);
 DO1 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ C1 <-> D1 ( r1 , d1 )*/
 f_flux =  r1 * C1 ;
 b_flux =  d1 * D1 ;
 DC1 -= (f_flux - b_flux);
 DD1 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ C2 <-> D2 ( r2 , d2 )*/
 f_flux =  r2 * C2 ;
 b_flux =  d2 * D2 ;
 DC2 -= (f_flux - b_flux);
 DD2 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ D1 <-> D2 ( rd , dd )*/
 f_flux =  rd * D1 ;
 b_flux =  dd * D2 ;
 DD1 -= (f_flux - b_flux);
 DD2 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ D2 <-> D3 ( r3 , d3 )*/
 f_flux =  r3 * D2 ;
 b_flux =  d3 * D3 ;
 DD2 -= (f_flux - b_flux);
 DD3 += (f_flux - b_flux);
 
 /*REACTION*/
   /* C0 + C1 + C2 + D1 + D2 + D3 + O1 = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE matsol*/
 static int _ode_matsol1() {_reset=0;{
 double b_flux, f_flux, _term; int _i;
   b_flux = f_flux = 0.;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<7;_i++){
  	_RHS1(_i) = _dt1*(_p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 f1 = k1 * ( 1e3 ) * XMTR ;
 f2 = k1 * ( 1e3 ) * XMTR ;
 /* ~ C0 <-> C1 ( f1 , km1 )*/
 _term =  f1 ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 2 ,3)  -= _term;
 _term =  km1 ;
 _MATELM1( 3 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  /* ~ C1 <-> C2 ( f2 , 12.5 * km1 )*/
 _term =  f2 ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 1 ,2)  -= _term;
 _term =  12.5 * km1 ;
 _MATELM1( 2 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ C2 <-> O1 ( a1 , b1 )*/
 _term =  a1 ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 0 ,1)  -= _term;
 _term =  b1 ;
 _MATELM1( 1 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
  /* ~ C1 <-> D1 ( r1 , d1 )*/
 _term =  r1 ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 6 ,2)  -= _term;
 _term =  d1 ;
 _MATELM1( 2 ,6)  -= _term;
 _MATELM1( 6 ,6)  += _term;
 /*REACTION*/
  /* ~ C2 <-> D2 ( r2 , d2 )*/
 _term =  r2 ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 5 ,1)  -= _term;
 _term =  d2 ;
 _MATELM1( 1 ,5)  -= _term;
 _MATELM1( 5 ,5)  += _term;
 /*REACTION*/
  /* ~ D1 <-> D2 ( rd , dd )*/
 _term =  rd ;
 _MATELM1( 6 ,6)  += _term;
 _MATELM1( 5 ,6)  -= _term;
 _term =  dd ;
 _MATELM1( 6 ,5)  -= _term;
 _MATELM1( 5 ,5)  += _term;
 /*REACTION*/
  /* ~ D2 <-> D3 ( r3 , d3 )*/
 _term =  r3 ;
 _MATELM1( 5 ,5)  += _term;
 _MATELM1( 4 ,5)  -= _term;
 _term =  d3 ;
 _MATELM1( 5 ,4)  -= _term;
 _MATELM1( 4 ,4)  += _term;
 /*REACTION*/
   /* C0 + C1 + C2 + D1 + D2 + D3 + O1 = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE end*/
 
static int _ode_count(_type) int _type;{ return 7;}
 
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
	for (_i=0; _i < 7; ++_i) {
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
 _cvode_sparse(&_cvsparseobj1, 7, _dlist1, _p, _ode_matsol1, &_coef1);
 }}

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  C2 = C20;
  C1 = C10;
  C0 = C00;
  D3 = D30;
  D2 = D20;
  D1 = D10;
  O1 = O10;
 {
   XMTR = 0.0 ;
   C0 = 1.0 ;
   C1 = 0.0 ;
   C2 = 0.0 ;
   O1 = 0.0 ;
   D1 = 0.0 ;
   D2 = 0.0 ;
   D3 = 0.0 ;
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
   Open = ( O1 ) ;
   g = gmax * Open ;
   i = ( 1e-6 ) * g * ( v - Erev ) ;
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
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
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
 error = sparse(&_sparseobj1, 7, _slist1, _dlist1, _p, &t, dt, kstates,&_coef1, _linmat1);
 if(error){fprintf(stderr,"at line 114 in file Gly5GC.mod:\n	SOLVE kstates METHOD sparse\n"); nrn_complain(_p); abort_run(error);}
 
}}
 t = _save;
 }}}

}

static terminal(){}

static _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(O1) - _p;  _dlist1[0] = &(DO1) - _p;
 _slist1[1] = &(C2) - _p;  _dlist1[1] = &(DC2) - _p;
 _slist1[2] = &(C1) - _p;  _dlist1[2] = &(DC1) - _p;
 _slist1[3] = &(C0) - _p;  _dlist1[3] = &(DC0) - _p;
 _slist1[4] = &(D3) - _p;  _dlist1[4] = &(DD3) - _p;
 _slist1[5] = &(D2) - _p;  _dlist1[5] = &(DD2) - _p;
 _slist1[6] = &(D1) - _p;  _dlist1[6] = &(DD1) - _p;
_first = 0;
}
