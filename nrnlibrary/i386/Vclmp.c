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
#define dur (_p + 0)
#define amp (_p + 5)
#define gain _p[10]
#define rstim _p[11]
#define tau1 _p[12]
#define tau2 _p[13]
#define e0 _p[14]
#define vo0 _p[15]
#define vi0 _p[16]
#define fac _p[17]
#define i _p[18]
#define e _p[19]
#define vo _p[20]
#define vi _p[21]
#define stim _p[22]
#define tc _p[23]
#define De _p[24]
#define Dvo _p[25]
#define Dvi _p[26]
#define _g _p[27]
#define _nd_area  *_ppvar[0]._pval
 
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
 static double _hoc_icur();
 static double _hoc_update();
 static double _hoc_vstim();
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
 "icur", _hoc_icur,
 "update", _hoc_update,
 "vstim", _hoc_vstim,
 0, 0
};
#define icur icur_VClamp2
 extern double icur();
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "dur", 0, 1e+09,
 "fac", 1, 10,
 "gain", 0, 1e+09,
 "rstim", 1e-09, 1e+09,
 "tau2", 0, 1e+09,
 "tau1", 0, 1e+09,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "dur", "ms",
 "amp", "mV",
 "rstim", "megohm",
 "tau1", "ms",
 "tau2", "ms",
 "e0", "mV",
 "vo0", "mV",
 "vi0", "mV",
 "fac", "1.0",
 "e", "mV",
 "vo", "mV",
 "vi", "mV",
 "i", "nA",
 0,0
};
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
 
static int _ode_count();
 /* connect range variables in _p that hoc is supposed to know about */
 static char *_mechanism[] = {
 "6.2.0",
"VClamp2",
 "dur[5]",
 "amp[5]",
 "gain",
 "rstim",
 "tau1",
 "tau2",
 "e0",
 "vo0",
 "vi0",
 "fac",
 0,
 "i",
 0,
 "e",
 "vo",
 "vi",
 0,
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
 	_p = nrn_prop_data_alloc(_mechtype, 28, _prop);
 	/*initialize range parameters*/
 	gain = 100000;
 	rstim = 1;
 	tau1 = 0.001;
 	tau2 = 0;
 	e0 = 0;
 	vo0 = 0;
 	vi0 = 0;
 	fac = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 28;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 2, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static _initlists();
 _Vclmp_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func,
	 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
  hoc_register_dparam_size(_mechtype, 2);
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 VClamp2 /Users/pbmanis/Desktop/Python/nrnlibrary/nrnlibrary/i386/Vclmp.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "vclmp2.mod";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
static update();
static vstim();
 
#define _RHS1(arg) _coef1[arg][3]
 static int _slist1[3];static double **_coef1;
 static int clamp();
 
static int  vstim (  )  {
   tc = 0.0 ;
   {int  _lj ;for ( _lj = 0 ; _lj <= 5 - 1 ; _lj ++ ) {
     stim = amp [ _lj ] ;
     tc = tc + dur [ _lj ] ;
     if ( t < tc ) {
       tc = tc + 100.0 ;
       
/*VERBATIM*/
			break;
 }
     } }
    return 0; }
 
static double _hoc_vstim(_vptr) void* _vptr; {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 vstim (  ) ;
 return(_r);
}
 
double icur (  )  {
   double _licur;
 double _lvout ;
 if ( t > tc ) {
     e0 = 0.0 ;
     vo0 = 0.0 ;
     _licur = 0.0 ;
     }
   else {
     error =  0; clamp();
 error = simeq(3, _coef1, _p, _slist1);
 if(error){fprintf(stderr,"at line 134 in file Vclmp.mod:\n		SOLVE clamp\n"); nrn_complain(_p); abort_run(error);}
 _licur = ( vo - v ) / rstim ;
     }
   
return _licur;
 }
 
static double _hoc_icur(_vptr) void* _vptr; {
 double _r;
    _hoc_setdata(_vptr);
 _r =  icur (  ) ;
 return(_r);
}
 
static int clamp ()
 {
   zero_matrix(_coef1, 3, 4);
{
  int _counte = -1;
 double _lt1 , _lt2 ;
 _lt1 = tau1 / dt ;
   _lt2 = tau2 / dt ;
    ++_counte;
 _coef1[_counte][0] -=  1.0 ;
 _RHS1(_counte) -=  v ;
 _coef1[_counte][1] +=  fac * 1.0 ;
 _RHS1(_counte) +=  fac * v ;
   ;
  ++_counte;
 _coef1[_counte][1] -=  _lt2 * 1.0 ;
 _RHS1(_counte) -=  _lt2 * vo0 ;
 _coef1[_counte][1] -=  1.0 ;
 _coef1[_counte][2] -=  gain * 1.0 ;
   ;
  ++_counte;
 _RHS1(_counte) -=  stim ;
 _coef1[_counte][2] +=  1.0 ;
 _coef1[_counte][0] -=  1.0 ;
 _coef1[_counte][2] +=  1.0 ;
 _coef1[_counte][0] -=  _lt1 * 1.0 ;
 _coef1[_counte][2] +=  _lt1 * 1.0 ;
 _RHS1(_counte) -=  _lt1 * ( vi0 - e0 ) ;
 _RHS1(_counte) -=  0.0 ;
   ;
 
}
 }
 
static int  update (  )  {
   i = icur ( _threadargs_ ) ;
   e0 = e ;
   vo0 = vo ;
   vi0 = vi ;
   
/*VERBATIM*/
	return 0;
  return 0; }
 
static double _hoc_update(_vptr) void* _vptr; {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 update (  ) ;
 return(_r);
}
 
static int _ode_count(_type)int _type; { hoc_execerror("VClamp2", "cannot be used with CVODE");}

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  e = e0;
  vi = vi0;
  vo = vo0;
 {
   e0 = 0.0 ;
   vo = v ;
   vo0 = v ;
   vi = v ;
   vi0 = v ;
   e = 0.0 ;
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
#if EXTRACELLULAR
 _nd = _ml->_nodelist[_iml];
 if (_nd->_extnode) {
    _v = NODEV(_nd) +_nd->_extnode->_v[0];
 }else
#endif
 {
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 }
 v = _v;
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   vstim ( _threadargs_ ) ;
   i = icur ( _threadargs_ ) ;
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
#if EXTRACELLULAR
 _nd = _ml->_nodelist[_iml];
 if (_nd->_extnode) {
    _v = NODEV(_nd) +_nd->_extnode->_v[0];
 }else
#endif
 {
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 }
 _g = _nrn_current(_v + .001);
 	{ _rhs = _nrn_current(_v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) += _rhs;
  }else
#endif
  {
	NODERHS(_nd) += _rhs;
  }
#if EXTRACELLULAR
 if (_nd->_extnode) {
   *_nd->_extnode->_rhs[0] += _rhs;
 }
#endif
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) -= _g;
  }else
#endif
  {
	NODED(_nd) -= _g;
  }
#if EXTRACELLULAR
 if (_nd->_extnode) {
   *_nd->_extnode->_d[0] += _g;
 }
#endif
 
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
#if EXTRACELLULAR
 _nd = _ml->_nodelist[_iml];
 if (_nd->_extnode) {
    _v = NODEV(_nd) +_nd->_extnode->_v[0];
 }else
#endif
 {
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 }
 _break = t + .5*dt; _save = t;
 v=_v;
{
 { {
 for (; t < _break; t += dt) {
 error =  update();
 if(error){fprintf(stderr,"at line 103 in file Vclmp.mod:\n	SOLVE update\n"); nrn_complain(_p); abort_run(error);}
 
}}
 t = _save;
 }}}

}

static terminal(){}

static _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(vi) - _p;
 _slist1[1] = &(vo) - _p;
 _slist1[2] = &(e) - _p;
 if (_first) _coef1 = makematrix(3, 4);
_first = 0;
}
