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
#define Q10 _p[2]
#define Mode _p[3]
#define zd _p[4]
#define PA _p[5]
#define Kd0 _p[6]
#define i _p[7]
#define g _p[8]
#define rb _p[9]
#define C0 _p[10]
#define C1 _p[11]
#define C2 _p[12]
#define D _p[13]
#define O1 _p[14]
#define O2 _p[15]
#define gvdep _p[16]
#define qfac _p[17]
#define DC0 _p[18]
#define DC1 _p[19]
#define DC2 _p[20]
#define DD _p[21]
#define DO1 _p[22]
#define DO2 _p[23]
#define _g _p[24]
#define _nd_area  *_ppvar[0]._pval
#define CGLU	*_ppvar[2]._pval
#define _p_CGLU	_ppvar[2]._pval
 
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
 extern double celsius;
 /* declaration of user functions */
 static double _hoc_gvdepcalc();
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
 "gvdepcalc", _hoc_gvdepcalc,
 0, 0
};
 /* declare global and static user variables */
#define Rc2 Rc2_AMPA_M
 double Rc2 = 0.25;
#define Ro2 Ro2_AMPA_M
 double Ro2 = 2;
#define Rc1 Rc1_AMPA_M
 double Rc1 = 2;
#define Ro1 Ro1_AMPA_M
 double Ro1 = 100;
#define Rr Rr_AMPA_M
 double Rr = 0.02;
#define Rd Rd_AMPA_M
 double Rd = 30;
#define Ru2 Ru2_AMPA_M
 double Ru2 = 200;
#define Ru1 Ru1_AMPA_M
 double Ru1 = 0.3;
#define Rb Rb_AMPA_M
 double Rb = 13;
#define aflag aflag_AMPA_M
 double aflag = 1;
#define usetable usetable_AMPA_M
 double usetable = 1;
#define vmax vmax_AMPA_M
 double vmax = 100;
#define vmin vmin_AMPA_M
 double vmin = -120;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_AMPA_M", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "vmin_AMPA_M", "mV",
 "vmax_AMPA_M", "mV",
 "Rb_AMPA_M", "/mM",
 "Ru1_AMPA_M", "/ms",
 "Ru2_AMPA_M", "/ms",
 "Rd_AMPA_M", "/ms",
 "Rr_AMPA_M", "/ms",
 "Ro1_AMPA_M", "/ms",
 "Rc1_AMPA_M", "/ms",
 "Ro2_AMPA_M", "/ms",
 "Rc2_AMPA_M", "/ms",
 "Erev", "mV",
 "gmax", "pS",
 "i", "nA",
 "g", "pS",
 "rb", "/ms",
 "CGLU", "mM",
 0,0
};
 static double C20 = 0;
 static double C10 = 0;
 static double C00 = 0;
 static double D0 = 0;
 static double O20 = 0;
 static double O10 = 0;
 static double delta_t = 1;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "vmin_AMPA_M", &vmin_AMPA_M,
 "vmax_AMPA_M", &vmax_AMPA_M,
 "Rb_AMPA_M", &Rb_AMPA_M,
 "Ru1_AMPA_M", &Ru1_AMPA_M,
 "Ru2_AMPA_M", &Ru2_AMPA_M,
 "Rd_AMPA_M", &Rd_AMPA_M,
 "Rr_AMPA_M", &Rr_AMPA_M,
 "Ro1_AMPA_M", &Ro1_AMPA_M,
 "Rc1_AMPA_M", &Rc1_AMPA_M,
 "Ro2_AMPA_M", &Ro2_AMPA_M,
 "Rc2_AMPA_M", &Rc2_AMPA_M,
 "aflag_AMPA_M", &aflag_AMPA_M,
 "usetable_AMPA_M", &usetable_AMPA_M,
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
"AMPA_M",
 "Erev",
 "gmax",
 "Q10",
 "Mode",
 "zd",
 "PA",
 "Kd0",
 0,
 "i",
 "g",
 "rb",
 0,
 "C0",
 "C1",
 "C2",
 "D",
 "O1",
 "O2",
 0,
 "CGLU",
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
 	_p = nrn_prop_data_alloc(_mechtype, 25, _prop);
 	/*initialize range parameters*/
 	Erev = 7;
 	gmax = 500;
 	Q10 = 1.5;
 	Mode = 1;
 	zd = 1.032;
 	PA = 45;
 	Kd0 = 3.1e-05;
  }
 	_prop->param = _p;
 	_prop->param_size = 25;
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
 _ampa_m_reg() {
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
 	ivoc_help("help ?1 AMPA_M /Users/pbmanis/Desktop/Python/nrnlibrary/nrnlibrary/i386/ampa_m.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double _zg0 ;
 static double *_t_gvdep;
static int _reset;
static char *modelname = "Model of AMPA receptors";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
static _f_gvdepcalc();
static gvdepcalc();
 extern double *_getelm();
 
#define _MATELM1(_row,_col)	*(_getelm(_row + 1, _col + 1))
 
#define _RHS1(_arg) _coef1[_arg + 1]
 static double *_coef1;
 
#define _linmat1  1
 static void* _sparseobj1;
 static void* _cvsparseobj1;
 
static int _ode_spec1(), _ode_matsol1();
 static _n_gvdepcalc();
 static int _slist1[6], _dlist1[6]; static double *_temp1;
 static int kstates();
 
static int kstates ()
 {_reset=0;
 {
   double b_flux, f_flux, _term; int _i;
 {int _i; double _dt1 = 1.0/dt;
for(_i=1;_i<6;_i++){
  	_RHS1(_i) = -_dt1*(_p[_slist1[_i]] - _p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 rb = Rb * CGLU ;
   /* ~ C0 <-> C1 ( rb * qfac , Ru1 * qfac )*/
 f_flux =  rb * qfac * C0 ;
 b_flux =  Ru1 * qfac * C1 ;
 _RHS1( 3) -= (f_flux - b_flux);
 _RHS1( 2) += (f_flux - b_flux);
 
 _term =  rb * qfac ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 2 ,3)  -= _term;
 _term =  Ru1 * qfac ;
 _MATELM1( 3 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  /* ~ C1 <-> C2 ( rb * qfac , Ru2 * qfac )*/
 f_flux =  rb * qfac * C1 ;
 b_flux =  Ru2 * qfac * C2 ;
 _RHS1( 2) -= (f_flux - b_flux);
 _RHS1( 1) += (f_flux - b_flux);
 
 _term =  rb * qfac ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 1 ,2)  -= _term;
 _term =  Ru2 * qfac ;
 _MATELM1( 2 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ C2 <-> D ( Rd * qfac , Rr * qfac )*/
 f_flux =  Rd * qfac * C2 ;
 b_flux =  Rr * qfac * D ;
 _RHS1( 1) -= (f_flux - b_flux);
 _RHS1( 4) += (f_flux - b_flux);
 
 _term =  Rd * qfac ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 4 ,1)  -= _term;
 _term =  Rr * qfac ;
 _MATELM1( 1 ,4)  -= _term;
 _MATELM1( 4 ,4)  += _term;
 /*REACTION*/
  /* ~ C2 <-> O1 ( Ro1 * qfac , Rc1 * qfac )*/
 f_flux =  Ro1 * qfac * C2 ;
 b_flux =  Rc1 * qfac * O1 ;
 _RHS1( 1) -= (f_flux - b_flux);
 _RHS1( 5) += (f_flux - b_flux);
 
 _term =  Ro1 * qfac ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 5 ,1)  -= _term;
 _term =  Rc1 * qfac ;
 _MATELM1( 1 ,5)  -= _term;
 _MATELM1( 5 ,5)  += _term;
 /*REACTION*/
  /* ~ C2 <-> O2 ( Ro2 * qfac , Rc2 * qfac )*/
 f_flux =  Ro2 * qfac * C2 ;
 b_flux =  Rc2 * qfac * O2 ;
 _RHS1( 1) -= (f_flux - b_flux);
 
 _term =  Ro2 * qfac ;
 _MATELM1( 1 ,1)  += _term;
 _term =  Rc2 * qfac ;
 _MATELM1( 1 ,0)  -= _term;
 /*REACTION*/
   /* C0 + C1 + C2 + D + O1 + O2 = 1.0 */
 _RHS1(0) =  1.0;
 _MATELM1(0, 0) = 1;
 _RHS1(0) -= O2 ;
 _MATELM1(0, 5) = 1;
 _RHS1(0) -= O1 ;
 _MATELM1(0, 4) = 1;
 _RHS1(0) -= D ;
 _MATELM1(0, 1) = 1;
 _RHS1(0) -= C2 ;
 _MATELM1(0, 2) = 1;
 _RHS1(0) -= C1 ;
 _MATELM1(0, 3) = 1;
 _RHS1(0) -= C0 ;
 /*CONSERVATION*/
   } return _reset;
 }
 static double _mfac_gvdepcalc, _tmin_gvdepcalc;
 static _check_gvdepcalc();
 static _check_gvdepcalc() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_PA;
  static double _sav_Kd0;
  static double _sav_zd;
  if (!usetable) {return;}
  if (_sav_PA != PA) { _maktable = 1;}
  if (_sav_Kd0 != Kd0) { _maktable = 1;}
  if (_sav_zd != zd) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_gvdepcalc =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_gvdepcalc)/200.; _mfac_gvdepcalc = 1./_dx;
   for (_i=0, _x=_tmin_gvdepcalc; _i < 201; _x += _dx, _i++) {
    _f_gvdepcalc(_x);
    _t_gvdep[_i] = gvdep;
   }
   _sav_PA = PA;
   _sav_Kd0 = Kd0;
   _sav_zd = zd;
  }
 }

 static gvdepcalc(double _lv){ _check_gvdepcalc();
 _n_gvdepcalc(_lv);
 return;
 }

 static _n_gvdepcalc(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_gvdepcalc(_lv); return; 
}
 _xi = _mfac_gvdepcalc * (_lv - _tmin_gvdepcalc);
 _i = (int) _xi;
 if (_xi <= 0.) {
 gvdep = _t_gvdep[0];
 return; }
 if (_i >= 200) {
 gvdep = _t_gvdep[200];
 return; }
 _theta = _xi - (double)_i;
 gvdep = _t_gvdep[_i] + _theta*(_t_gvdep[_i+1] - _t_gvdep[_i]);
 }

 
static int  _f_gvdepcalc (  _lv )  
	double _lv ;
 {
   _zg0 = 1.0 + 0.6 * exp ( ( _lv - 50.0 ) / 40.0 ) ;
   gvdep = _zg0 * ( 1.0 / ( 1.0 + PA / ( Kd0 * exp ( - zd * _lv / 25.3 ) ) ) ) ;
    return 0; }
 
static double _hoc_gvdepcalc(_vptr) void* _vptr; {
 double _r;
    _hoc_setdata(_vptr);
  _r = 1.;
 gvdepcalc (  *getarg(1) ) ;
 return(_r);
}
 
/*CVODE ode begin*/
 static int _ode_spec1() {_reset=0;{
 double b_flux, f_flux, _term; int _i;
 {int _i; for(_i=0;_i<6;_i++) _p[_dlist1[_i]] = 0.0;}
 rb = Rb * CGLU ;
 /* ~ C0 <-> C1 ( rb * qfac , Ru1 * qfac )*/
 f_flux =  rb * qfac * C0 ;
 b_flux =  Ru1 * qfac * C1 ;
 DC0 -= (f_flux - b_flux);
 DC1 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ C1 <-> C2 ( rb * qfac , Ru2 * qfac )*/
 f_flux =  rb * qfac * C1 ;
 b_flux =  Ru2 * qfac * C2 ;
 DC1 -= (f_flux - b_flux);
 DC2 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ C2 <-> D ( Rd * qfac , Rr * qfac )*/
 f_flux =  Rd * qfac * C2 ;
 b_flux =  Rr * qfac * D ;
 DC2 -= (f_flux - b_flux);
 DD += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ C2 <-> O1 ( Ro1 * qfac , Rc1 * qfac )*/
 f_flux =  Ro1 * qfac * C2 ;
 b_flux =  Rc1 * qfac * O1 ;
 DC2 -= (f_flux - b_flux);
 DO1 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ C2 <-> O2 ( Ro2 * qfac , Rc2 * qfac )*/
 f_flux =  Ro2 * qfac * C2 ;
 b_flux =  Rc2 * qfac * O2 ;
 DC2 -= (f_flux - b_flux);
 DO2 += (f_flux - b_flux);
 
 /*REACTION*/
   /* C0 + C1 + C2 + D + O1 + O2 = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE matsol*/
 static int _ode_matsol1() {_reset=0;{
 double b_flux, f_flux, _term; int _i;
   b_flux = f_flux = 0.;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<6;_i++){
  	_RHS1(_i) = _dt1*(_p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 rb = Rb * CGLU ;
 /* ~ C0 <-> C1 ( rb * qfac , Ru1 * qfac )*/
 _term =  rb * qfac ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 2 ,3)  -= _term;
 _term =  Ru1 * qfac ;
 _MATELM1( 3 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  /* ~ C1 <-> C2 ( rb * qfac , Ru2 * qfac )*/
 _term =  rb * qfac ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 1 ,2)  -= _term;
 _term =  Ru2 * qfac ;
 _MATELM1( 2 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ C2 <-> D ( Rd * qfac , Rr * qfac )*/
 _term =  Rd * qfac ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 4 ,1)  -= _term;
 _term =  Rr * qfac ;
 _MATELM1( 1 ,4)  -= _term;
 _MATELM1( 4 ,4)  += _term;
 /*REACTION*/
  /* ~ C2 <-> O1 ( Ro1 * qfac , Rc1 * qfac )*/
 _term =  Ro1 * qfac ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 5 ,1)  -= _term;
 _term =  Rc1 * qfac ;
 _MATELM1( 1 ,5)  -= _term;
 _MATELM1( 5 ,5)  += _term;
 /*REACTION*/
  /* ~ C2 <-> O2 ( Ro2 * qfac , Rc2 * qfac )*/
 _term =  Ro2 * qfac ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 0 ,1)  -= _term;
 _term =  Rc2 * qfac ;
 _MATELM1( 1 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
   /* C0 + C1 + C2 + D + O1 + O2 = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE end*/
 
static int _ode_count(_type) int _type;{ return 6;}
 
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
	for (_i=0; _i < 6; ++_i) {
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
 _cvode_sparse(&_cvsparseobj1, 6, _dlist1, _p, _ode_matsol1, &_coef1);
 }}

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  C2 = C20;
  C1 = C10;
  C0 = C00;
  D = D0;
  O2 = O20;
  O1 = O10;
 {
   usetable = 0.0 ;
   C0 = 1.0 ;
   C1 = 0.0 ;
   C2 = 0.0 ;
   D = 0.0 ;
   O1 = 0.0 ;
   O2 = 0.0 ;
   qfac = pow( Q10 , ( ( celsius - 22.0 ) / 10.0 ) ) ;
   gvdepcalc ( _threadargscomma_ v ) ;
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
   gvdepcalc ( _threadargscomma_ v ) ;
   g = gmax * ( O1 + O2 ) ;
   if ( Mode  == 1.0 ) {
     i = ( 1e-6 ) * g * gvdep * ( v - Erev ) ;
     }
   else {
     i = ( 1e-6 ) * g * ( v - Erev ) ;
     }
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
 error = sparse(&_sparseobj1, 6, _slist1, _dlist1, _p, &t, dt, kstates,&_coef1, _linmat1);
 if(error){fprintf(stderr,"at line 162 in file ampa_m.mod:\n    SOLVE kstates METHOD sparse\n"); nrn_complain(_p); abort_run(error);}
 
}}
 t = _save;
 }}}

}

static terminal(){}

static _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
   _t_gvdep = makevector(201*sizeof(double));
 _slist1[0] = &(O2) - _p;  _dlist1[0] = &(DO2) - _p;
 _slist1[1] = &(C2) - _p;  _dlist1[1] = &(DC2) - _p;
 _slist1[2] = &(C1) - _p;  _dlist1[2] = &(DC1) - _p;
 _slist1[3] = &(C0) - _p;  _dlist1[3] = &(DC0) - _p;
 _slist1[4] = &(D) - _p;  _dlist1[4] = &(DD) - _p;
 _slist1[5] = &(O1) - _p;  _dlist1[5] = &(DO1) - _p;
_first = 0;
}
