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
#define gkltbar _p[0]
#define ik _p[1]
#define gklt _p[2]
#define w _p[3]
#define z _p[4]
#define ek _p[5]
#define Dw _p[6]
#define Dz _p[7]
#define _g _p[8]
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
 "setdata_klt", _hoc_setdata,
 "rates_klt", _hoc_rates,
 "states_klt", _hoc_states,
 "trates_klt", _hoc_trates,
 "vtrap_klt", _hoc_vtrap,
 0, 0
};
#define vtrap vtrap_klt
 extern double vtrap();
 /* declare global and static user variables */
#define usetable usetable_klt
 double usetable = 1;
#define wtau wtau_klt
 double wtau = 0;
#define winf winf_klt
 double winf = 0;
#define zss zss_klt
 double zss = 0.5;
#define ztau ztau_klt
 double ztau = 0;
#define zinf zinf_klt
 double zinf = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "gkltbar_klt", 0, 1e+09,
 "usetable_klt", 0, 1,
 "zss_klt", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "wtau_klt", "ms",
 "ztau_klt", "ms",
 "gkltbar_klt", "mho/cm2",
 "ik_klt", "mA/cm2",
 "gklt_klt", "mho/cm2",
 0,0
};
 static double delta_t = 1;
 static double v = 0;
 static double w0 = 0;
 static double z0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "zss_klt", &zss_klt,
 "winf_klt", &winf_klt,
 "zinf_klt", &zinf_klt,
 "wtau_klt", &wtau_klt,
 "ztau_klt", &ztau_klt,
 "usetable_klt", &usetable_klt,
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
"klt",
 "gkltbar_klt",
 0,
 "ik_klt",
 "gklt_klt",
 0,
 "w_klt",
 "z_klt",
 0,
 0};
 static Symbol* _k_sym;
 
static void nrn_alloc(_prop)
	Prop *_prop;
{
	Prop *prop_ion, *need_memb();
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 9, _prop);
 	/*initialize range parameters*/
 	gkltbar = 0.01592;
 	_prop->param = _p;
 	_prop->param_size = 9;
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
 _klt_reg() {
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
 	ivoc_help("help ?1 klt /Users/pbmanis/Desktop/Python/nrnlibrary/nrnlibrary/i386/klt.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double _zwexp , _zzexp ;
 static double _zq10 ;
 static double *_t_winf;
 static double *_t__zwexp;
 static double *_t_zinf;
 static double *_t__zzexp;
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
   w = w + _zwexp * ( winf - w ) ;
   z = z + _zzexp * ( zinf - z ) ;
   
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
   winf = pow( ( 1.0 / ( 1.0 + exp ( - ( _lv + 48.0 ) / 6.0 ) ) ) , 0.25 ) ;
   zinf = zss + ( ( 1.0 - zss ) / ( 1.0 + exp ( ( _lv + 71.0 ) / 10.0 ) ) ) ;
   wtau = ( 100.0 / ( 6.0 * exp ( ( _lv + 60.0 ) / 6.0 ) + 16.0 * exp ( - ( _lv + 60.0 ) / 45.0 ) ) ) + 1.5 ;
   ztau = ( 1000.0 / ( exp ( ( _lv + 60.0 ) / 20.0 ) + exp ( - ( _lv + 60.0 ) / 8.0 ) ) ) + 50.0 ;
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
    _t_winf[_i] = winf;
    _t__zwexp[_i] = _zwexp;
    _t_zinf[_i] = zinf;
    _t__zzexp[_i] = _zzexp;
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
 winf = _t_winf[0];
 _zwexp = _t__zwexp[0];
 zinf = _t_zinf[0];
 _zzexp = _t__zzexp[0];
 return; }
 if (_i >= 300) {
 winf = _t_winf[300];
 _zwexp = _t__zwexp[300];
 zinf = _t_zinf[300];
 _zzexp = _t__zzexp[300];
 return; }
 _theta = _xi - (double)_i;
 winf = _t_winf[_i] + _theta*(_t_winf[_i+1] - _t_winf[_i]);
 _zwexp = _t__zwexp[_i] + _theta*(_t__zwexp[_i+1] - _t__zwexp[_i]);
 zinf = _t_zinf[_i] + _theta*(_t_zinf[_i+1] - _t_zinf[_i]);
 _zzexp = _t__zzexp[_i] + _theta*(_t__zzexp[_i+1] - _t__zzexp[_i]);
 }

 
static int  _f_trates (  _lv )  
	double _lv ;
 {
   double _ltinc ;
 rates ( _threadargscomma_ _lv ) ;
   _ltinc = - dt * _zq10 ;
   _zwexp = 1.0 - exp ( _ltinc / wtau ) ;
   _zzexp = 1.0 - exp ( _ltinc / ztau ) ;
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
 
static int _ode_count(_type)int _type; { hoc_execerror("klt", "cannot be used with CVODE");}
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
  w = w0;
  z = z0;
 {
   trates ( _threadargscomma_ v ) ;
   w = winf ;
   z = zinf ;
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
   gklt = gkltbar * ( pow( w , 4.0 ) ) * z ;
   ik = gklt * ( v - ek ) ;
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
 if(error){fprintf(stderr,"at line 73 in file klt.mod:\n\n"); nrn_complain(_p); abort_run(error);}
 
}}
 t = _save;
 } }}

}

static terminal(){}

static _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
   _t_winf = makevector(301*sizeof(double));
   _t__zwexp = makevector(301*sizeof(double));
   _t_zinf = makevector(301*sizeof(double));
   _t__zzexp = makevector(301*sizeof(double));
_first = 0;
}
