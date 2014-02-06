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
#define gNAR _p[2]
#define vshift _p[3]
#define i _p[4]
#define g _p[5]
#define rb _p[6]
#define rmb _p[7]
#define rmu _p[8]
#define rbMg _p[9]
#define rmc1b _p[10]
#define rmc1u _p[11]
#define rmc2b _p[12]
#define rmc2u _p[13]
#define U _p[14]
#define Cl _p[15]
#define D1 _p[16]
#define D2 _p[17]
#define Open _p[18]
#define UMg _p[19]
#define ClMg _p[20]
#define D1Mg _p[21]
#define D2Mg _p[22]
#define OMg _p[23]
#define qfac _p[24]
#define DU _p[25]
#define DCl _p[26]
#define DD1 _p[27]
#define DD2 _p[28]
#define DOpen _p[29]
#define DUMg _p[30]
#define DClMg _p[31]
#define DD1Mg _p[32]
#define DD2Mg _p[33]
#define DOMg _p[34]
#define _g _p[35]
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
 extern double celsius;
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
#define Erev Erev_NMDA_Kampa
 double Erev = 5;
#define Q10 Q10_NMDA_Kampa
 double Q10 = 2;
#define Rmc2u Rmc2u_NMDA_Kampa
 double Rmc2u = 0.00504192;
#define Rmc2b Rmc2b_NMDA_Kampa
 double Rmc2b = 5e-08;
#define Rmc1u Rmc1u_NMDA_Kampa
 double Rmc1u = 0.00243831;
#define Rmc1b Rmc1b_NMDA_Kampa
 double Rmc1b = 5e-08;
#define Rmd2u Rmd2u_NMDA_Kampa
 double Rmd2u = 0.00295341;
#define Rmd2b Rmd2b_NMDA_Kampa
 double Rmd2b = 5e-08;
#define Rmd1u Rmd1u_NMDA_Kampa
 double Rmd1u = 0.00298874;
#define Rmd1b Rmd1b_NMDA_Kampa
 double Rmd1b = 5e-08;
#define RcMg RcMg_NMDA_Kampa
 double RcMg = 0.548;
#define RoMg RoMg_NMDA_Kampa
 double RoMg = 0.01;
#define Rr2Mg Rr2Mg_NMDA_Kampa
 double Rr2Mg = 0.00042;
#define Rd2Mg Rd2Mg_NMDA_Kampa
 double Rd2Mg = 0.00026;
#define Rr1Mg Rr1Mg_NMDA_Kampa
 double Rr1Mg = 0.00087;
#define Rd1Mg Rd1Mg_NMDA_Kampa
 double Rd1Mg = 0.0021;
#define RuMg RuMg_NMDA_Kampa
 double RuMg = 0.0171;
#define RbMg RbMg_NMDA_Kampa
 double RbMg = 0.01;
#define Rmu Rmu_NMDA_Kampa
 double Rmu = 12.8;
#define Rmb Rmb_NMDA_Kampa
 double Rmb = 5e-05;
#define Rc Rc_NMDA_Kampa
 double Rc = 0.273;
#define Ro Ro_NMDA_Kampa
 double Ro = 0.01;
#define Rr2 Rr2_NMDA_Kampa
 double Rr2 = 0.0005;
#define Rd2 Rd2_NMDA_Kampa
 double Rd2 = 0.0001;
#define Rr1 Rr1_NMDA_Kampa
 double Rr1 = 0.0016;
#define Rd1 Rd1_NMDA_Kampa
 double Rd1 = 0.1;
#define Ru Ru_NMDA_Kampa
 double Ru = 0.0056;
#define Rb Rb_NMDA_Kampa
 double Rb = 0.01;
#define memb_fraction memb_fraction_NMDA_Kampa
 double memb_fraction = 0.8;
#define mg mg_NMDA_Kampa
 double mg = 1;
#define rmd2u rmd2u_NMDA_Kampa
 double rmd2u = 0;
#define rmd2b rmd2b_NMDA_Kampa
 double rmd2b = 0;
#define rmd1u rmd1u_NMDA_Kampa
 double rmd1u = 0;
#define rmd1b rmd1b_NMDA_Kampa
 double rmd1b = 0;
#define valence valence_NMDA_Kampa
 double valence = -2;
#define vmax vmax_NMDA_Kampa
 double vmax = 100;
#define vmin vmin_NMDA_Kampa
 double vmin = -120;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "Erev_NMDA_Kampa", "mV",
 "mg_NMDA_Kampa", "mM",
 "vmin_NMDA_Kampa", "mV",
 "vmax_NMDA_Kampa", "mV",
 "Rb_NMDA_Kampa", "/uM",
 "Ru_NMDA_Kampa", "/ms",
 "Ro_NMDA_Kampa", "/ms",
 "Rc_NMDA_Kampa", "/ms",
 "Rd1_NMDA_Kampa", "/ms",
 "Rr1_NMDA_Kampa", "/ms",
 "Rd2_NMDA_Kampa", "/ms",
 "Rr2_NMDA_Kampa", "/ms",
 "Rmb_NMDA_Kampa", "/uM",
 "Rmu_NMDA_Kampa", "/ms",
 "Rmc1b_NMDA_Kampa", "/uM",
 "Rmc1u_NMDA_Kampa", "/ms",
 "Rmc2b_NMDA_Kampa", "/uM",
 "Rmc2u_NMDA_Kampa", "/ms",
 "Rmd1b_NMDA_Kampa", "/uM",
 "Rmd1u_NMDA_Kampa", "/ms",
 "Rmd2b_NMDA_Kampa", "/uM",
 "Rmd2u_NMDA_Kampa", "/ms",
 "RbMg_NMDA_Kampa", "/uM",
 "RuMg_NMDA_Kampa", "/ms",
 "RoMg_NMDA_Kampa", "/ms",
 "RcMg_NMDA_Kampa", "/ms",
 "Rd1Mg_NMDA_Kampa", "/ms",
 "Rr1Mg_NMDA_Kampa", "/ms",
 "Rd2Mg_NMDA_Kampa", "/ms",
 "Rr2Mg_NMDA_Kampa", "/ms",
 "rmd1b_NMDA_Kampa", "/ms",
 "rmd1u_NMDA_Kampa", "/ms",
 "rmd2b_NMDA_Kampa", "/ms",
 "rmd2u_NMDA_Kampa", "/ms",
 "Erev", "mV",
 "gmax", "pS",
 "i", "nA",
 "g", "pS",
 "rb", "/ms",
 "rmb", "/ms",
 "rmu", "/ms",
 "rbMg", "/ms",
 "rmc1b", "/ms",
 "rmc1u", "/ms",
 "rmc2b", "/ms",
 "rmc2u", "/ms",
 "XMTR", "mM",
 0,0
};
 static double ClMg0 = 0;
 static double Cl0 = 0;
 static double D2Mg0 = 0;
 static double D1Mg0 = 0;
 static double D20 = 0;
 static double D10 = 0;
 static double OMg0 = 0;
 static double Open0 = 0;
 static double UMg0 = 0;
 static double U0 = 0;
 static double delta_t = 1;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "Erev_NMDA_Kampa", &Erev_NMDA_Kampa,
 "mg_NMDA_Kampa", &mg_NMDA_Kampa,
 "vmin_NMDA_Kampa", &vmin_NMDA_Kampa,
 "vmax_NMDA_Kampa", &vmax_NMDA_Kampa,
 "valence_NMDA_Kampa", &valence_NMDA_Kampa,
 "memb_fraction_NMDA_Kampa", &memb_fraction_NMDA_Kampa,
 "Q10_NMDA_Kampa", &Q10_NMDA_Kampa,
 "Rb_NMDA_Kampa", &Rb_NMDA_Kampa,
 "Ru_NMDA_Kampa", &Ru_NMDA_Kampa,
 "Ro_NMDA_Kampa", &Ro_NMDA_Kampa,
 "Rc_NMDA_Kampa", &Rc_NMDA_Kampa,
 "Rd1_NMDA_Kampa", &Rd1_NMDA_Kampa,
 "Rr1_NMDA_Kampa", &Rr1_NMDA_Kampa,
 "Rd2_NMDA_Kampa", &Rd2_NMDA_Kampa,
 "Rr2_NMDA_Kampa", &Rr2_NMDA_Kampa,
 "Rmb_NMDA_Kampa", &Rmb_NMDA_Kampa,
 "Rmu_NMDA_Kampa", &Rmu_NMDA_Kampa,
 "Rmc1b_NMDA_Kampa", &Rmc1b_NMDA_Kampa,
 "Rmc1u_NMDA_Kampa", &Rmc1u_NMDA_Kampa,
 "Rmc2b_NMDA_Kampa", &Rmc2b_NMDA_Kampa,
 "Rmc2u_NMDA_Kampa", &Rmc2u_NMDA_Kampa,
 "Rmd1b_NMDA_Kampa", &Rmd1b_NMDA_Kampa,
 "Rmd1u_NMDA_Kampa", &Rmd1u_NMDA_Kampa,
 "Rmd2b_NMDA_Kampa", &Rmd2b_NMDA_Kampa,
 "Rmd2u_NMDA_Kampa", &Rmd2u_NMDA_Kampa,
 "RbMg_NMDA_Kampa", &RbMg_NMDA_Kampa,
 "RuMg_NMDA_Kampa", &RuMg_NMDA_Kampa,
 "RoMg_NMDA_Kampa", &RoMg_NMDA_Kampa,
 "RcMg_NMDA_Kampa", &RcMg_NMDA_Kampa,
 "Rd1Mg_NMDA_Kampa", &Rd1Mg_NMDA_Kampa,
 "Rr1Mg_NMDA_Kampa", &Rr1Mg_NMDA_Kampa,
 "Rd2Mg_NMDA_Kampa", &Rd2Mg_NMDA_Kampa,
 "Rr2Mg_NMDA_Kampa", &Rr2Mg_NMDA_Kampa,
 "rmd1b_NMDA_Kampa", &rmd1b_NMDA_Kampa,
 "rmd1u_NMDA_Kampa", &rmd1u_NMDA_Kampa,
 "rmd2b_NMDA_Kampa", &rmd2b_NMDA_Kampa,
 "rmd2u_NMDA_Kampa", &rmd2u_NMDA_Kampa,
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
"NMDA_Kampa",
 "Erev",
 "gmax",
 "gNAR",
 "vshift",
 0,
 "i",
 "g",
 "rb",
 "rmb",
 "rmu",
 "rbMg",
 "rmc1b",
 "rmc1u",
 "rmc2b",
 "rmc2u",
 0,
 "U",
 "Cl",
 "D1",
 "D2",
 "Open",
 "UMg",
 "ClMg",
 "D1Mg",
 "D2Mg",
 "OMg",
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
 	_p = nrn_prop_data_alloc(_mechtype, 36, _prop);
 	/*initialize range parameters*/
 	Erev = 5;
 	gmax = 500;
 	gNAR = 0.5;
 	vshift = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 36;
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
 _NMDA_Kampa_reg() {
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
 	ivoc_help("help ?1 NMDA_Kampa /Users/pbmanis/Desktop/Python/nrnlibrary/nrnlibrary/i386/NMDA_Kampa.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "kinetic NMDA receptor model";

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
 static int _slist1[10], _dlist1[10]; static double *_temp1;
 static int kstates();
 
static int kstates ()
 {_reset=0;
 {
   double b_flux, f_flux, _term; int _i;
 {int _i; double _dt1 = 1.0/dt;
for(_i=1;_i<10;_i++){
  	_RHS1(_i) = -_dt1*(_p[_slist1[_i]] - _p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 rb = Rb * ( 1e3 ) * XMTR ;
   rbMg = RbMg * ( 1e3 ) * XMTR ;
   rmb = Rmb * mg * ( 1e3 ) * exp ( ( v - 40.0 + vshift ) * valence * memb_fraction / 25.0 ) ;
   rmu = Rmu * exp ( ( - 1.0 ) * ( v - 40.0 + vshift ) * valence * ( 1.0 - memb_fraction ) / 25.0 ) ;
   rmc1b = Rmc1b * mg * ( 1e3 ) * exp ( ( v - 40.0 + vshift ) * valence * memb_fraction / 25.0 ) ;
   rmc1u = Rmc1u * exp ( ( - 1.0 ) * ( v - 40.0 + vshift ) * valence * ( 1.0 - memb_fraction ) / 25.0 ) ;
   rmc2b = Rmc2b * mg * ( 1e3 ) * exp ( ( v - 40.0 + vshift ) * valence * memb_fraction / 25.0 ) ;
   rmc2u = Rmc2u * exp ( ( - 1.0 ) * ( v - 40.0 + vshift ) * valence * ( 1.0 - memb_fraction ) / 25.0 ) ;
   rmd1b = Rmd1b * mg * ( 1e3 ) * exp ( ( v - 40.0 + vshift ) * valence * memb_fraction / 25.0 ) ;
   rmd1u = Rmd1u * exp ( ( - 1.0 ) * ( v - 40.0 + vshift ) * valence * ( 1.0 - memb_fraction ) / 25.0 ) ;
   rmd2b = Rmd2b * mg * ( 1e3 ) * exp ( ( v - 40.0 + vshift ) * valence * memb_fraction / 25.0 ) ;
   rmd2u = Rmd2u * exp ( ( - 1.0 ) * ( v - 40.0 + vshift ) * valence * ( 1.0 - memb_fraction ) / 25.0 ) ;
   /* ~ U <-> Cl ( rb * qfac , Ru * qfac )*/
 f_flux =  rb * qfac * U ;
 b_flux =  Ru * qfac * Cl ;
 _RHS1( 9) -= (f_flux - b_flux);
 _RHS1( 2) += (f_flux - b_flux);
 
 _term =  rb * qfac ;
 _MATELM1( 9 ,9)  += _term;
 _MATELM1( 2 ,9)  -= _term;
 _term =  Ru * qfac ;
 _MATELM1( 9 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  /* ~ Cl <-> Open ( Ro * qfac , Rc * qfac )*/
 f_flux =  Ro * qfac * Cl ;
 b_flux =  Rc * qfac * Open ;
 _RHS1( 2) -= (f_flux - b_flux);
 _RHS1( 7) += (f_flux - b_flux);
 
 _term =  Ro * qfac ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 7 ,2)  -= _term;
 _term =  Rc * qfac ;
 _MATELM1( 2 ,7)  -= _term;
 _MATELM1( 7 ,7)  += _term;
 /*REACTION*/
  /* ~ Cl <-> D1 ( Rd1 * qfac , Rr1 * qfac )*/
 f_flux =  Rd1 * qfac * Cl ;
 b_flux =  Rr1 * qfac * D1 ;
 _RHS1( 2) -= (f_flux - b_flux);
 _RHS1( 6) += (f_flux - b_flux);
 
 _term =  Rd1 * qfac ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 6 ,2)  -= _term;
 _term =  Rr1 * qfac ;
 _MATELM1( 2 ,6)  -= _term;
 _MATELM1( 6 ,6)  += _term;
 /*REACTION*/
  /* ~ D1 <-> D2 ( Rd2 * qfac , Rr2 * qfac )*/
 f_flux =  Rd2 * qfac * D1 ;
 b_flux =  Rr2 * qfac * D2 ;
 _RHS1( 6) -= (f_flux - b_flux);
 _RHS1( 5) += (f_flux - b_flux);
 
 _term =  Rd2 * qfac ;
 _MATELM1( 6 ,6)  += _term;
 _MATELM1( 5 ,6)  -= _term;
 _term =  Rr2 * qfac ;
 _MATELM1( 6 ,5)  -= _term;
 _MATELM1( 5 ,5)  += _term;
 /*REACTION*/
  /* ~ Open <-> OMg ( rmb * qfac , rmu * qfac )*/
 f_flux =  rmb * qfac * Open ;
 b_flux =  rmu * qfac * OMg ;
 _RHS1( 7) -= (f_flux - b_flux);
 
 _term =  rmb * qfac ;
 _MATELM1( 7 ,7)  += _term;
 _term =  rmu * qfac ;
 _MATELM1( 7 ,0)  -= _term;
 /*REACTION*/
  /* ~ UMg <-> ClMg ( rbMg * qfac , RuMg * qfac )*/
 f_flux =  rbMg * qfac * UMg ;
 b_flux =  RuMg * qfac * ClMg ;
 _RHS1( 8) -= (f_flux - b_flux);
 _RHS1( 1) += (f_flux - b_flux);
 
 _term =  rbMg * qfac ;
 _MATELM1( 8 ,8)  += _term;
 _MATELM1( 1 ,8)  -= _term;
 _term =  RuMg * qfac ;
 _MATELM1( 8 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ ClMg <-> OMg ( RoMg * qfac , RcMg * qfac )*/
 f_flux =  RoMg * qfac * ClMg ;
 b_flux =  RcMg * qfac * OMg ;
 _RHS1( 1) -= (f_flux - b_flux);
 
 _term =  RoMg * qfac ;
 _MATELM1( 1 ,1)  += _term;
 _term =  RcMg * qfac ;
 _MATELM1( 1 ,0)  -= _term;
 /*REACTION*/
  /* ~ ClMg <-> D1Mg ( Rd1Mg * qfac , Rr1Mg * qfac )*/
 f_flux =  Rd1Mg * qfac * ClMg ;
 b_flux =  Rr1Mg * qfac * D1Mg ;
 _RHS1( 1) -= (f_flux - b_flux);
 _RHS1( 4) += (f_flux - b_flux);
 
 _term =  Rd1Mg * qfac ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 4 ,1)  -= _term;
 _term =  Rr1Mg * qfac ;
 _MATELM1( 1 ,4)  -= _term;
 _MATELM1( 4 ,4)  += _term;
 /*REACTION*/
  /* ~ D1Mg <-> D2Mg ( Rd2Mg * qfac , Rr2Mg * qfac )*/
 f_flux =  Rd2Mg * qfac * D1Mg ;
 b_flux =  Rr2Mg * qfac * D2Mg ;
 _RHS1( 4) -= (f_flux - b_flux);
 _RHS1( 3) += (f_flux - b_flux);
 
 _term =  Rd2Mg * qfac ;
 _MATELM1( 4 ,4)  += _term;
 _MATELM1( 3 ,4)  -= _term;
 _term =  Rr2Mg * qfac ;
 _MATELM1( 4 ,3)  -= _term;
 _MATELM1( 3 ,3)  += _term;
 /*REACTION*/
  /* ~ U <-> UMg ( rmc1b * qfac , rmc1u * qfac )*/
 f_flux =  rmc1b * qfac * U ;
 b_flux =  rmc1u * qfac * UMg ;
 _RHS1( 9) -= (f_flux - b_flux);
 _RHS1( 8) += (f_flux - b_flux);
 
 _term =  rmc1b * qfac ;
 _MATELM1( 9 ,9)  += _term;
 _MATELM1( 8 ,9)  -= _term;
 _term =  rmc1u * qfac ;
 _MATELM1( 9 ,8)  -= _term;
 _MATELM1( 8 ,8)  += _term;
 /*REACTION*/
  /* ~ Cl <-> ClMg ( rmc2b * qfac , rmc2u * qfac )*/
 f_flux =  rmc2b * qfac * Cl ;
 b_flux =  rmc2u * qfac * ClMg ;
 _RHS1( 2) -= (f_flux - b_flux);
 _RHS1( 1) += (f_flux - b_flux);
 
 _term =  rmc2b * qfac ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 1 ,2)  -= _term;
 _term =  rmc2u * qfac ;
 _MATELM1( 2 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ D1 <-> D1Mg ( rmd1b * qfac , rmd1u * qfac )*/
 f_flux =  rmd1b * qfac * D1 ;
 b_flux =  rmd1u * qfac * D1Mg ;
 _RHS1( 6) -= (f_flux - b_flux);
 _RHS1( 4) += (f_flux - b_flux);
 
 _term =  rmd1b * qfac ;
 _MATELM1( 6 ,6)  += _term;
 _MATELM1( 4 ,6)  -= _term;
 _term =  rmd1u * qfac ;
 _MATELM1( 6 ,4)  -= _term;
 _MATELM1( 4 ,4)  += _term;
 /*REACTION*/
  /* ~ D2 <-> D2Mg ( rmd2b * qfac , rmd2u * qfac )*/
 f_flux =  rmd2b * qfac * D2 ;
 b_flux =  rmd2u * qfac * D2Mg ;
 _RHS1( 5) -= (f_flux - b_flux);
 _RHS1( 3) += (f_flux - b_flux);
 
 _term =  rmd2b * qfac ;
 _MATELM1( 5 ,5)  += _term;
 _MATELM1( 3 ,5)  -= _term;
 _term =  rmd2u * qfac ;
 _MATELM1( 5 ,3)  -= _term;
 _MATELM1( 3 ,3)  += _term;
 /*REACTION*/
   /* U + Cl + D1 + D2 + Open + UMg + ClMg + D1Mg + D2Mg + OMg = 1.0 */
 _RHS1(0) =  1.0;
 _MATELM1(0, 0) = 1;
 _RHS1(0) -= OMg ;
 _MATELM1(0, 3) = 1;
 _RHS1(0) -= D2Mg ;
 _MATELM1(0, 4) = 1;
 _RHS1(0) -= D1Mg ;
 _MATELM1(0, 1) = 1;
 _RHS1(0) -= ClMg ;
 _MATELM1(0, 8) = 1;
 _RHS1(0) -= UMg ;
 _MATELM1(0, 7) = 1;
 _RHS1(0) -= Open ;
 _MATELM1(0, 5) = 1;
 _RHS1(0) -= D2 ;
 _MATELM1(0, 6) = 1;
 _RHS1(0) -= D1 ;
 _MATELM1(0, 2) = 1;
 _RHS1(0) -= Cl ;
 _MATELM1(0, 9) = 1;
 _RHS1(0) -= U ;
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE ode begin*/
 static int _ode_spec1() {_reset=0;{
 double b_flux, f_flux, _term; int _i;
 {int _i; for(_i=0;_i<10;_i++) _p[_dlist1[_i]] = 0.0;}
 rb = Rb * ( 1e3 ) * XMTR ;
 rbMg = RbMg * ( 1e3 ) * XMTR ;
 rmb = Rmb * mg * ( 1e3 ) * exp ( ( v - 40.0 + vshift ) * valence * memb_fraction / 25.0 ) ;
 rmu = Rmu * exp ( ( - 1.0 ) * ( v - 40.0 + vshift ) * valence * ( 1.0 - memb_fraction ) / 25.0 ) ;
 rmc1b = Rmc1b * mg * ( 1e3 ) * exp ( ( v - 40.0 + vshift ) * valence * memb_fraction / 25.0 ) ;
 rmc1u = Rmc1u * exp ( ( - 1.0 ) * ( v - 40.0 + vshift ) * valence * ( 1.0 - memb_fraction ) / 25.0 ) ;
 rmc2b = Rmc2b * mg * ( 1e3 ) * exp ( ( v - 40.0 + vshift ) * valence * memb_fraction / 25.0 ) ;
 rmc2u = Rmc2u * exp ( ( - 1.0 ) * ( v - 40.0 + vshift ) * valence * ( 1.0 - memb_fraction ) / 25.0 ) ;
 rmd1b = Rmd1b * mg * ( 1e3 ) * exp ( ( v - 40.0 + vshift ) * valence * memb_fraction / 25.0 ) ;
 rmd1u = Rmd1u * exp ( ( - 1.0 ) * ( v - 40.0 + vshift ) * valence * ( 1.0 - memb_fraction ) / 25.0 ) ;
 rmd2b = Rmd2b * mg * ( 1e3 ) * exp ( ( v - 40.0 + vshift ) * valence * memb_fraction / 25.0 ) ;
 rmd2u = Rmd2u * exp ( ( - 1.0 ) * ( v - 40.0 + vshift ) * valence * ( 1.0 - memb_fraction ) / 25.0 ) ;
 /* ~ U <-> Cl ( rb * qfac , Ru * qfac )*/
 f_flux =  rb * qfac * U ;
 b_flux =  Ru * qfac * Cl ;
 DU -= (f_flux - b_flux);
 DCl += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ Cl <-> Open ( Ro * qfac , Rc * qfac )*/
 f_flux =  Ro * qfac * Cl ;
 b_flux =  Rc * qfac * Open ;
 DCl -= (f_flux - b_flux);
 DOpen += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ Cl <-> D1 ( Rd1 * qfac , Rr1 * qfac )*/
 f_flux =  Rd1 * qfac * Cl ;
 b_flux =  Rr1 * qfac * D1 ;
 DCl -= (f_flux - b_flux);
 DD1 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ D1 <-> D2 ( Rd2 * qfac , Rr2 * qfac )*/
 f_flux =  Rd2 * qfac * D1 ;
 b_flux =  Rr2 * qfac * D2 ;
 DD1 -= (f_flux - b_flux);
 DD2 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ Open <-> OMg ( rmb * qfac , rmu * qfac )*/
 f_flux =  rmb * qfac * Open ;
 b_flux =  rmu * qfac * OMg ;
 DOpen -= (f_flux - b_flux);
 DOMg += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ UMg <-> ClMg ( rbMg * qfac , RuMg * qfac )*/
 f_flux =  rbMg * qfac * UMg ;
 b_flux =  RuMg * qfac * ClMg ;
 DUMg -= (f_flux - b_flux);
 DClMg += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ ClMg <-> OMg ( RoMg * qfac , RcMg * qfac )*/
 f_flux =  RoMg * qfac * ClMg ;
 b_flux =  RcMg * qfac * OMg ;
 DClMg -= (f_flux - b_flux);
 DOMg += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ ClMg <-> D1Mg ( Rd1Mg * qfac , Rr1Mg * qfac )*/
 f_flux =  Rd1Mg * qfac * ClMg ;
 b_flux =  Rr1Mg * qfac * D1Mg ;
 DClMg -= (f_flux - b_flux);
 DD1Mg += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ D1Mg <-> D2Mg ( Rd2Mg * qfac , Rr2Mg * qfac )*/
 f_flux =  Rd2Mg * qfac * D1Mg ;
 b_flux =  Rr2Mg * qfac * D2Mg ;
 DD1Mg -= (f_flux - b_flux);
 DD2Mg += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ U <-> UMg ( rmc1b * qfac , rmc1u * qfac )*/
 f_flux =  rmc1b * qfac * U ;
 b_flux =  rmc1u * qfac * UMg ;
 DU -= (f_flux - b_flux);
 DUMg += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ Cl <-> ClMg ( rmc2b * qfac , rmc2u * qfac )*/
 f_flux =  rmc2b * qfac * Cl ;
 b_flux =  rmc2u * qfac * ClMg ;
 DCl -= (f_flux - b_flux);
 DClMg += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ D1 <-> D1Mg ( rmd1b * qfac , rmd1u * qfac )*/
 f_flux =  rmd1b * qfac * D1 ;
 b_flux =  rmd1u * qfac * D1Mg ;
 DD1 -= (f_flux - b_flux);
 DD1Mg += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ D2 <-> D2Mg ( rmd2b * qfac , rmd2u * qfac )*/
 f_flux =  rmd2b * qfac * D2 ;
 b_flux =  rmd2u * qfac * D2Mg ;
 DD2 -= (f_flux - b_flux);
 DD2Mg += (f_flux - b_flux);
 
 /*REACTION*/
   /* U + Cl + D1 + D2 + Open + UMg + ClMg + D1Mg + D2Mg + OMg = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE matsol*/
 static int _ode_matsol1() {_reset=0;{
 double b_flux, f_flux, _term; int _i;
   b_flux = f_flux = 0.;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<10;_i++){
  	_RHS1(_i) = _dt1*(_p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 rb = Rb * ( 1e3 ) * XMTR ;
 rbMg = RbMg * ( 1e3 ) * XMTR ;
 rmb = Rmb * mg * ( 1e3 ) * exp ( ( v - 40.0 + vshift ) * valence * memb_fraction / 25.0 ) ;
 rmu = Rmu * exp ( ( - 1.0 ) * ( v - 40.0 + vshift ) * valence * ( 1.0 - memb_fraction ) / 25.0 ) ;
 rmc1b = Rmc1b * mg * ( 1e3 ) * exp ( ( v - 40.0 + vshift ) * valence * memb_fraction / 25.0 ) ;
 rmc1u = Rmc1u * exp ( ( - 1.0 ) * ( v - 40.0 + vshift ) * valence * ( 1.0 - memb_fraction ) / 25.0 ) ;
 rmc2b = Rmc2b * mg * ( 1e3 ) * exp ( ( v - 40.0 + vshift ) * valence * memb_fraction / 25.0 ) ;
 rmc2u = Rmc2u * exp ( ( - 1.0 ) * ( v - 40.0 + vshift ) * valence * ( 1.0 - memb_fraction ) / 25.0 ) ;
 rmd1b = Rmd1b * mg * ( 1e3 ) * exp ( ( v - 40.0 + vshift ) * valence * memb_fraction / 25.0 ) ;
 rmd1u = Rmd1u * exp ( ( - 1.0 ) * ( v - 40.0 + vshift ) * valence * ( 1.0 - memb_fraction ) / 25.0 ) ;
 rmd2b = Rmd2b * mg * ( 1e3 ) * exp ( ( v - 40.0 + vshift ) * valence * memb_fraction / 25.0 ) ;
 rmd2u = Rmd2u * exp ( ( - 1.0 ) * ( v - 40.0 + vshift ) * valence * ( 1.0 - memb_fraction ) / 25.0 ) ;
 /* ~ U <-> Cl ( rb * qfac , Ru * qfac )*/
 _term =  rb * qfac ;
 _MATELM1( 9 ,9)  += _term;
 _MATELM1( 2 ,9)  -= _term;
 _term =  Ru * qfac ;
 _MATELM1( 9 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  /* ~ Cl <-> Open ( Ro * qfac , Rc * qfac )*/
 _term =  Ro * qfac ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 7 ,2)  -= _term;
 _term =  Rc * qfac ;
 _MATELM1( 2 ,7)  -= _term;
 _MATELM1( 7 ,7)  += _term;
 /*REACTION*/
  /* ~ Cl <-> D1 ( Rd1 * qfac , Rr1 * qfac )*/
 _term =  Rd1 * qfac ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 6 ,2)  -= _term;
 _term =  Rr1 * qfac ;
 _MATELM1( 2 ,6)  -= _term;
 _MATELM1( 6 ,6)  += _term;
 /*REACTION*/
  /* ~ D1 <-> D2 ( Rd2 * qfac , Rr2 * qfac )*/
 _term =  Rd2 * qfac ;
 _MATELM1( 6 ,6)  += _term;
 _MATELM1( 5 ,6)  -= _term;
 _term =  Rr2 * qfac ;
 _MATELM1( 6 ,5)  -= _term;
 _MATELM1( 5 ,5)  += _term;
 /*REACTION*/
  /* ~ Open <-> OMg ( rmb * qfac , rmu * qfac )*/
 _term =  rmb * qfac ;
 _MATELM1( 7 ,7)  += _term;
 _MATELM1( 0 ,7)  -= _term;
 _term =  rmu * qfac ;
 _MATELM1( 7 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
  /* ~ UMg <-> ClMg ( rbMg * qfac , RuMg * qfac )*/
 _term =  rbMg * qfac ;
 _MATELM1( 8 ,8)  += _term;
 _MATELM1( 1 ,8)  -= _term;
 _term =  RuMg * qfac ;
 _MATELM1( 8 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ ClMg <-> OMg ( RoMg * qfac , RcMg * qfac )*/
 _term =  RoMg * qfac ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 0 ,1)  -= _term;
 _term =  RcMg * qfac ;
 _MATELM1( 1 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
  /* ~ ClMg <-> D1Mg ( Rd1Mg * qfac , Rr1Mg * qfac )*/
 _term =  Rd1Mg * qfac ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 4 ,1)  -= _term;
 _term =  Rr1Mg * qfac ;
 _MATELM1( 1 ,4)  -= _term;
 _MATELM1( 4 ,4)  += _term;
 /*REACTION*/
  /* ~ D1Mg <-> D2Mg ( Rd2Mg * qfac , Rr2Mg * qfac )*/
 _term =  Rd2Mg * qfac ;
 _MATELM1( 4 ,4)  += _term;
 _MATELM1( 3 ,4)  -= _term;
 _term =  Rr2Mg * qfac ;
 _MATELM1( 4 ,3)  -= _term;
 _MATELM1( 3 ,3)  += _term;
 /*REACTION*/
  /* ~ U <-> UMg ( rmc1b * qfac , rmc1u * qfac )*/
 _term =  rmc1b * qfac ;
 _MATELM1( 9 ,9)  += _term;
 _MATELM1( 8 ,9)  -= _term;
 _term =  rmc1u * qfac ;
 _MATELM1( 9 ,8)  -= _term;
 _MATELM1( 8 ,8)  += _term;
 /*REACTION*/
  /* ~ Cl <-> ClMg ( rmc2b * qfac , rmc2u * qfac )*/
 _term =  rmc2b * qfac ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 1 ,2)  -= _term;
 _term =  rmc2u * qfac ;
 _MATELM1( 2 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ D1 <-> D1Mg ( rmd1b * qfac , rmd1u * qfac )*/
 _term =  rmd1b * qfac ;
 _MATELM1( 6 ,6)  += _term;
 _MATELM1( 4 ,6)  -= _term;
 _term =  rmd1u * qfac ;
 _MATELM1( 6 ,4)  -= _term;
 _MATELM1( 4 ,4)  += _term;
 /*REACTION*/
  /* ~ D2 <-> D2Mg ( rmd2b * qfac , rmd2u * qfac )*/
 _term =  rmd2b * qfac ;
 _MATELM1( 5 ,5)  += _term;
 _MATELM1( 3 ,5)  -= _term;
 _term =  rmd2u * qfac ;
 _MATELM1( 5 ,3)  -= _term;
 _MATELM1( 3 ,3)  += _term;
 /*REACTION*/
   /* U + Cl + D1 + D2 + Open + UMg + ClMg + D1Mg + D2Mg + OMg = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE end*/
 
static int _ode_count(_type) int _type;{ return 10;}
 
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
	for (_i=0; _i < 10; ++_i) {
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
 _cvode_sparse(&_cvsparseobj1, 10, _dlist1, _p, _ode_matsol1, &_coef1);
 }}

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  ClMg = ClMg0;
  Cl = Cl0;
  D2Mg = D2Mg0;
  D1Mg = D1Mg0;
  D2 = D20;
  D1 = D10;
  OMg = OMg0;
  Open = Open0;
  UMg = UMg0;
  U = U0;
 {
   U = 1.0 ;
   qfac = pow( Q10 , ( ( celsius - 23.0 ) / 10.0 ) ) ;
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
   g = gNAR * gmax * Open ;
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
 error = sparse(&_sparseobj1, 10, _slist1, _dlist1, _p, &t, dt, kstates,&_coef1, _linmat1);
 if(error){fprintf(stderr,"at line 170 in file NMDA_Kampa.mod:\n\n"); nrn_complain(_p); abort_run(error);}
 
}}
 t = _save;
 }}}

}

static terminal(){}

static _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(OMg) - _p;  _dlist1[0] = &(DOMg) - _p;
 _slist1[1] = &(ClMg) - _p;  _dlist1[1] = &(DClMg) - _p;
 _slist1[2] = &(Cl) - _p;  _dlist1[2] = &(DCl) - _p;
 _slist1[3] = &(D2Mg) - _p;  _dlist1[3] = &(DD2Mg) - _p;
 _slist1[4] = &(D1Mg) - _p;  _dlist1[4] = &(DD1Mg) - _p;
 _slist1[5] = &(D2) - _p;  _dlist1[5] = &(DD2) - _p;
 _slist1[6] = &(D1) - _p;  _dlist1[6] = &(DD1) - _p;
 _slist1[7] = &(Open) - _p;  _dlist1[7] = &(DOpen) - _p;
 _slist1[8] = &(UMg) - _p;  _dlist1[8] = &(DUMg) - _p;
 _slist1[9] = &(U) - _p;  _dlist1[9] = &(DU) - _p;
_first = 0;
}
