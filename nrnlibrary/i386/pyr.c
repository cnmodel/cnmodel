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
#define gnabar _p[0]
#define gnapbar _p[1]
#define gkbar _p[2]
#define gkifbar _p[3]
#define gkisbar _p[4]
#define ghbar _p[5]
#define ghvshift _p[6]
#define gl _p[7]
#define el _p[8]
#define ntau _p[9]
#define kd_avh _p[10]
#define nap_A _p[11]
#define nap_B _p[12]
#define nap_C _p[13]
#define nap_D _p[14]
#define nap_E _p[15]
#define nap_F _p[16]
#define nap_G _p[17]
#define nap_H _p[18]
#define kif_ivh _p[19]
#define kif_avh _p[20]
#define kis_ivh _p[21]
#define kis_avh _p[22]
#define kif_a_start _p[23]
#define kif_i_start _p[24]
#define kis_a_start _p[25]
#define kis_i_start _p[26]
#define gna _p[27]
#define gk _p[28]
#define gkif _p[29]
#define gkis _p[30]
#define gh _p[31]
#define il _p[32]
#define minf _p[33]
#define hinf _p[34]
#define mtau _p[35]
#define htau _p[36]
#define ninf _p[37]
#define nap_inf _p[38]
#define nap_tau _p[39]
#define napi_inf _p[40]
#define napi_tau _p[41]
#define kif_a_inf _p[42]
#define kif_i_inf _p[43]
#define kis_a_inf _p[44]
#define kis_i_inf _p[45]
#define kif_a_tau _p[46]
#define kif_i_tau _p[47]
#define kis_a_tau _p[48]
#define kis_i_tau _p[49]
#define kh_m_inf _p[50]
#define kh_n_inf _p[51]
#define kh_m_tau _p[52]
#define kh_n_tau _p[53]
#define akif _p[54]
#define akis _p[55]
#define aih _p[56]
#define m _p[57]
#define h _p[58]
#define n _p[59]
#define kifa _p[60]
#define kifi _p[61]
#define kisa _p[62]
#define kisi _p[63]
#define khm _p[64]
#define khn _p[65]
#define nap _p[66]
#define napi _p[67]
#define ek _p[68]
#define ena _p[69]
#define eh _p[70]
#define Dm _p[71]
#define Dh _p[72]
#define Dn _p[73]
#define Dkifa _p[74]
#define Dkifi _p[75]
#define Dkisa _p[76]
#define Dkisi _p[77]
#define Dkhm _p[78]
#define Dkhn _p[79]
#define Dnap _p[80]
#define Dnapi _p[81]
#define gnap _p[82]
#define ina _p[83]
#define ik _p[84]
#define ih _p[85]
#define v _p[86]
#define _g _p[87]
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
#define _ion_ek	*_ppvar[3]._pval
#define _ion_ik	*_ppvar[4]._pval
#define _ion_dikdv	*_ppvar[5]._pval
#define _ion_eh	*_ppvar[6]._pval
#define _ion_ih	*_ppvar[7]._pval
#define _ion_dihdv	*_ppvar[8]._pval
 
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
 static int _hoc_alphbet();
 static int _hoc_kh_nt();
 static int _hoc_kh_mt();
 static int _hoc_kh_n();
 static int _hoc_kh_m();
 static int _hoc_kis_ht();
 static int _hoc_kis_mt();
 static int _hoc_kis_h();
 static int _hoc_kis_m();
 static int _hoc_kif_ht();
 static int _hoc_kif_mt();
 static int _hoc_kif_h();
 static int _hoc_kif_m();
 static int _hoc_kd_m();
 static int _hoc_na_pit();
 static int _hoc_na_pi();
 static int _hoc_na_pt();
 static int _hoc_na_p();
 static int _hoc_na_ht();
 static int _hoc_na_h();
 static int _hoc_na_mt();
 static int _hoc_na_m();
 static int _hoc_rates();
 static int _hoc_vtrap();
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
 "setdata_pyr", _hoc_setdata,
 "alphbet_pyr", _hoc_alphbet,
 "kh_nt_pyr", _hoc_kh_nt,
 "kh_mt_pyr", _hoc_kh_mt,
 "kh_n_pyr", _hoc_kh_n,
 "kh_m_pyr", _hoc_kh_m,
 "kis_ht_pyr", _hoc_kis_ht,
 "kis_mt_pyr", _hoc_kis_mt,
 "kis_h_pyr", _hoc_kis_h,
 "kis_m_pyr", _hoc_kis_m,
 "kif_ht_pyr", _hoc_kif_ht,
 "kif_mt_pyr", _hoc_kif_mt,
 "kif_h_pyr", _hoc_kif_h,
 "kif_m_pyr", _hoc_kif_m,
 "kd_m_pyr", _hoc_kd_m,
 "na_pit_pyr", _hoc_na_pit,
 "na_pi_pyr", _hoc_na_pi,
 "na_pt_pyr", _hoc_na_pt,
 "na_p_pyr", _hoc_na_p,
 "na_ht_pyr", _hoc_na_ht,
 "na_h_pyr", _hoc_na_h,
 "na_mt_pyr", _hoc_na_mt,
 "na_m_pyr", _hoc_na_m,
 "rates_pyr", _hoc_rates,
 "vtrap_pyr", _hoc_vtrap,
 0, 0
};
#define alphbet alphbet_pyr
#define kh_nt kh_nt_pyr
#define kh_mt kh_mt_pyr
#define kh_n kh_n_pyr
#define kh_m kh_m_pyr
#define kis_ht kis_ht_pyr
#define kis_mt kis_mt_pyr
#define kis_h kis_h_pyr
#define kis_m kis_m_pyr
#define kif_ht kif_ht_pyr
#define kif_mt kif_mt_pyr
#define kif_h kif_h_pyr
#define kif_m kif_m_pyr
#define kd_m kd_m_pyr
#define na_pit na_pit_pyr
#define na_pi na_pi_pyr
#define na_pt na_pt_pyr
#define na_p na_p_pyr
#define na_ht na_ht_pyr
#define na_h na_h_pyr
#define na_mt na_mt_pyr
#define na_m na_m_pyr
#define vtrap vtrap_pyr
 extern double alphbet();
 extern double kh_nt();
 extern double kh_mt();
 extern double kh_n();
 extern double kh_m();
 extern double kis_ht();
 extern double kis_mt();
 extern double kis_h();
 extern double kis_m();
 extern double kif_ht();
 extern double kif_mt();
 extern double kif_h();
 extern double kif_m();
 extern double kd_m();
 extern double na_pit();
 extern double na_pi();
 extern double na_pt();
 extern double na_p();
 extern double na_ht();
 extern double na_h();
 extern double na_mt();
 extern double na_m();
 extern double vtrap();
 
static void _check_rates(double*, Datum*, Datum*, _NrnThread*); 
static void _check_table_thread(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, int _type) {
   _check_rates(_p, _ppvar, _thread, _nt);
 }
 #define _zmexp _thread[0]._pval[0]
 #define _zhexp _thread[0]._pval[1]
 #define _znexp _thread[0]._pval[2]
 #define _znap_exp _thread[0]._pval[3]
 #define _znapi_exp _thread[0]._pval[4]
 #define _zkif_a_exp _thread[0]._pval[5]
 #define _zkif_i_exp _thread[0]._pval[6]
 #define _zkis_a_exp _thread[0]._pval[7]
 #define _zkis_i_exp _thread[0]._pval[8]
 #define _zkh_m_exp _thread[0]._pval[9]
 #define _zkh_n_exp _thread[0]._pval[10]
 #define _zq10 _thread[0]._pval[11]
 /* declare global and static user variables */
#define htau0 htau0_pyr
 double htau0 = 0.5;
#define kif_hivh kif_hivh_pyr
 double kif_hivh = -87;
#define mtau0 mtau0_pyr
 double mtau0 = 0.05;
#define nap_shift nap_shift_pyr
 double nap_shift = 0;
#define usetable usetable_pyr
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "gl_pyr", 0, 1e+09,
 "ghbar_pyr", 0, 1e+09,
 "gkifbar_pyr", 0, 1e+09,
 "gnapbar_pyr", 0, 1e+09,
 "gkbar_pyr", 0, 1e+09,
 "gnabar_pyr", 0, 1e+09,
 "htau0_pyr", 0.1, 100,
 "mtau0_pyr", 0.01, 100,
 "ntau_pyr", 0.1, 100,
 "usetable_pyr", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "mtau0_pyr", "ms",
 "htau0_pyr", "ms",
 "gnabar_pyr", "mho/cm2",
 "gnapbar_pyr", "mho/cm2",
 "gkbar_pyr", "mho/cm2",
 "gkifbar_pyr", "mho/cm2",
 "gkisbar_pyr", "mho/cm2",
 "ghbar_pyr", "mho/cm2",
 "ghvshift_pyr", "mV",
 "gl_pyr", "mho/cm2",
 "el_pyr", "mV",
 "ntau_pyr", "ms",
 "gna_pyr", "mho/cm2",
 "gk_pyr", "mho/cm2",
 "gkif_pyr", "mho/cm2",
 "gkis_pyr", "mho/cm2",
 "gh_pyr", "mho/cm2",
 "il_pyr", "mA/cm2",
 0,0
};
 static double delta_t = 1;
 static double h0 = 0;
 static double khn0 = 0;
 static double khm0 = 0;
 static double kisi0 = 0;
 static double kisa0 = 0;
 static double kifi0 = 0;
 static double kifa0 = 0;
 static double m0 = 0;
 static double napi0 = 0;
 static double nap0 = 0;
 static double n0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "mtau0_pyr", &mtau0_pyr,
 "htau0_pyr", &htau0_pyr,
 "nap_shift_pyr", &nap_shift_pyr,
 "kif_hivh_pyr", &kif_hivh_pyr,
 "usetable_pyr", &usetable_pyr,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(), nrn_init(), nrn_state();
 static void nrn_cur(), nrn_jacob();
 
static int _ode_count(), _ode_map(), _ode_spec(), _ode_matsol();
 
#define _cvode_ieq _ppvar[9]._i
 /* connect range variables in _p that hoc is supposed to know about */
 static char *_mechanism[] = {
 "6.2.0",
"pyr",
 "gnabar_pyr",
 "gnapbar_pyr",
 "gkbar_pyr",
 "gkifbar_pyr",
 "gkisbar_pyr",
 "ghbar_pyr",
 "ghvshift_pyr",
 "gl_pyr",
 "el_pyr",
 "ntau_pyr",
 "kd_avh_pyr",
 "nap_A_pyr",
 "nap_B_pyr",
 "nap_C_pyr",
 "nap_D_pyr",
 "nap_E_pyr",
 "nap_F_pyr",
 "nap_G_pyr",
 "nap_H_pyr",
 "kif_ivh_pyr",
 "kif_avh_pyr",
 "kis_ivh_pyr",
 "kis_avh_pyr",
 "kif_a_start_pyr",
 "kif_i_start_pyr",
 "kis_a_start_pyr",
 "kis_i_start_pyr",
 0,
 "gna_pyr",
 "gk_pyr",
 "gkif_pyr",
 "gkis_pyr",
 "gh_pyr",
 "il_pyr",
 "minf_pyr",
 "hinf_pyr",
 "mtau_pyr",
 "htau_pyr",
 "ninf_pyr",
 "nap_inf_pyr",
 "nap_tau_pyr",
 "napi_inf_pyr",
 "napi_tau_pyr",
 "kif_a_inf_pyr",
 "kif_i_inf_pyr",
 "kis_a_inf_pyr",
 "kis_i_inf_pyr",
 "kif_a_tau_pyr",
 "kif_i_tau_pyr",
 "kis_a_tau_pyr",
 "kis_i_tau_pyr",
 "kh_m_inf_pyr",
 "kh_n_inf_pyr",
 "kh_m_tau_pyr",
 "kh_n_tau_pyr",
 "akif_pyr",
 "akis_pyr",
 "aih_pyr",
 0,
 "m_pyr",
 "h_pyr",
 "n_pyr",
 "kifa_pyr",
 "kifi_pyr",
 "kisa_pyr",
 "kisi_pyr",
 "khm_pyr",
 "khn_pyr",
 "nap_pyr",
 "napi_pyr",
 0,
 0};
 static Symbol* _na_sym;
 static Symbol* _k_sym;
 static Symbol* _h_sym;
 
static void nrn_alloc(_prop)
	Prop *_prop;
{
	Prop *prop_ion, *need_memb();
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 88, _prop);
 	/*initialize range parameters*/
 	gnabar = 0.02857;
 	gnapbar = 1e-05;
 	gkbar = 0.006667;
 	gkifbar = 0.0125;
 	gkisbar = 0.0033333;
 	ghbar = 0.00025;
 	ghvshift = 0;
 	gl = 0.00025;
 	el = -57.7;
 	ntau = 0.5;
 	kd_avh = -40;
 	nap_A = 200;
 	nap_B = 1;
 	nap_C = -18;
 	nap_D = -16;
 	nap_E = 25;
 	nap_F = 1;
 	nap_G = 58;
 	nap_H = 8;
 	kif_ivh = -89.6;
 	kif_avh = -57;
 	kis_ivh = -40.9;
 	kis_avh = -38.4;
 	kif_a_start = -1;
 	kif_i_start = -1;
 	kis_a_start = -1;
 	kis_i_start = -1;
 	_prop->param = _p;
 	_prop->param_size = 88;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 10, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[3]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[4]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[5]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 prop_ion = need_memb(_h_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[6]._pval = &prop_ion->param[0]; /* eh */
 	_ppvar[7]._pval = &prop_ion->param[3]; /* ih */
 	_ppvar[8]._pval = &prop_ion->param[4]; /* _ion_dihdv */
 
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
 _pyr_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("na", -10000.);
 	ion_reg("k", -10000.);
 	ion_reg("h", 1.0);
 	_na_sym = hoc_lookup("na_ion");
 	_k_sym = hoc_lookup("k_ion");
 	_h_sym = hoc_lookup("h_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 2);
  _extcall_thread = (Datum*)ecalloc(1, sizeof(Datum));
  _thread_mem_init(_extcall_thread);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_thread_reg(_mechtype, 1, _thread_mem_init);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
     _nrn_thread_table_reg(_mechtype, _check_table_thread);
  hoc_register_dparam_size(_mechtype, 10);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 pyr /Users/pbmanis/Desktop/Python/nrnlibrary/nrnlibrary/i386/pyr.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 /*Top LOCAL _zmexp , _zhexp , _znexp , _znap_exp , _znapi_exp , _zkif_a_exp , _zkif_i_exp , _zkis_a_exp , _zkis_i_exp , _zkh_m_exp , _zkh_n_exp */
 /*Top LOCAL _zq10 */
 static double *_t_minf;
 static double *_t_mtau;
 static double *_t_hinf;
 static double *_t_htau;
 static double *_t_nap_inf;
 static double *_t_nap_tau;
 static double *_t_napi_inf;
 static double *_t_napi_tau;
 static double *_t_ninf;
 static double *_t_ntau;
 static double *_t_kif_a_inf;
 static double *_t_kif_a_tau;
 static double *_t_kif_i_inf;
 static double *_t_kif_i_tau;
 static double *_t_kis_a_inf;
 static double *_t_kis_a_tau;
 static double *_t_kis_i_inf;
 static double *_t_kis_i_tau;
 static double *_t_kh_m_inf;
 static double *_t_kh_n_inf;
 static double *_t_kh_m_tau;
 static double *_t_kh_n_tau;
static int _reset;
static char *modelname = "pyr.mod   DCN pyramidal cell model  ";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
static _f_rates();
static rates();
 
static int _ode_spec1(), _ode_matsol1();
 static _n_rates();
 static int _slist1[11], _dlist1[11];
 static int states();
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / mtau ;
   Dh = ( hinf - h ) / htau ;
   Dnap = ( nap_inf - nap ) / nap_tau ;
   Dnapi = ( napi_inf - napi ) / napi_tau ;
   Dn = ( ninf - n ) / ntau ;
   Dkifa = ( kif_a_inf - kifa ) / kif_a_tau ;
   Dkifi = ( kif_i_inf - kifi ) / kif_i_tau ;
   Dkisa = ( kis_a_inf - kisa ) / kis_a_tau ;
   Dkisi = ( kis_i_inf - kisi ) / kis_i_tau ;
   Dkhm = ( kh_m_inf - khm ) / kh_m_tau ;
   Dkhn = ( kh_n_inf - khn ) / kh_n_tau ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 rates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau )) ;
 Dnap = Dnap  / (1. - dt*( ( ( ( - 1.0 ) ) ) / nap_tau )) ;
 Dnapi = Dnapi  / (1. - dt*( ( ( ( - 1.0 ) ) ) / napi_tau )) ;
 Dn = Dn  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ntau )) ;
 Dkifa = Dkifa  / (1. - dt*( ( ( ( - 1.0 ) ) ) / kif_a_tau )) ;
 Dkifi = Dkifi  / (1. - dt*( ( ( ( - 1.0 ) ) ) / kif_i_tau )) ;
 Dkisa = Dkisa  / (1. - dt*( ( ( ( - 1.0 ) ) ) / kis_a_tau )) ;
 Dkisi = Dkisi  / (1. - dt*( ( ( ( - 1.0 ) ) ) / kis_i_tau )) ;
 Dkhm = Dkhm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / kh_m_tau )) ;
 Dkhn = Dkhn  / (1. - dt*( ( ( ( - 1.0 ) ) ) / kh_n_tau )) ;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   rates ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / mtau)))*(- ( ( ( minf ) ) / mtau ) / ( ( ( ( - 1.0) ) ) / mtau ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / htau)))*(- ( ( ( hinf ) ) / htau ) / ( ( ( ( - 1.0) ) ) / htau ) - h) ;
    nap = nap + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / nap_tau)))*(- ( ( ( nap_inf ) ) / nap_tau ) / ( ( ( ( - 1.0) ) ) / nap_tau ) - nap) ;
    napi = napi + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / napi_tau)))*(- ( ( ( napi_inf ) ) / napi_tau ) / ( ( ( ( - 1.0) ) ) / napi_tau ) - napi) ;
    n = n + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ntau)))*(- ( ( ( ninf ) ) / ntau ) / ( ( ( ( - 1.0) ) ) / ntau ) - n) ;
    kifa = kifa + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / kif_a_tau)))*(- ( ( ( kif_a_inf ) ) / kif_a_tau ) / ( ( ( ( - 1.0) ) ) / kif_a_tau ) - kifa) ;
    kifi = kifi + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / kif_i_tau)))*(- ( ( ( kif_i_inf ) ) / kif_i_tau ) / ( ( ( ( - 1.0) ) ) / kif_i_tau ) - kifi) ;
    kisa = kisa + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / kis_a_tau)))*(- ( ( ( kis_a_inf ) ) / kis_a_tau ) / ( ( ( ( - 1.0) ) ) / kis_a_tau ) - kisa) ;
    kisi = kisi + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / kis_i_tau)))*(- ( ( ( kis_i_inf ) ) / kis_i_tau ) / ( ( ( ( - 1.0) ) ) / kis_i_tau ) - kisi) ;
    khm = khm + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / kh_m_tau)))*(- ( ( ( kh_m_inf ) ) / kh_m_tau ) / ( ( ( ( - 1.0) ) ) / kh_m_tau ) - khm) ;
    khn = khn + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / kh_n_tau)))*(- ( ( ( kh_n_inf ) ) / kh_n_tau ) / ( ( ( ( - 1.0) ) ) / kh_n_tau ) - khn) ;
   }
  return 0;
}
 static double _mfac_rates, _tmin_rates;
  static void _check_rates(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_celsius;
  static double _sav_kif_avh;
  static double _sav_kif_ivh;
  static double _sav_kis_avh;
  static double _sav_kis_ivh;
  static double _sav_kd_avh;
  static double _sav_nap_A;
  static double _sav_nap_C;
  static double _sav_nap_D;
  static double _sav_nap_E;
  static double _sav_nap_G;
  static double _sav_nap_H;
  if (!usetable) {return;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_sav_kif_avh != kif_avh) { _maktable = 1;}
  if (_sav_kif_ivh != kif_ivh) { _maktable = 1;}
  if (_sav_kis_avh != kis_avh) { _maktable = 1;}
  if (_sav_kis_ivh != kis_ivh) { _maktable = 1;}
  if (_sav_kd_avh != kd_avh) { _maktable = 1;}
  if (_sav_nap_A != nap_A) { _maktable = 1;}
  if (_sav_nap_C != nap_C) { _maktable = 1;}
  if (_sav_nap_D != nap_D) { _maktable = 1;}
  if (_sav_nap_E != nap_E) { _maktable = 1;}
  if (_sav_nap_G != nap_G) { _maktable = 1;}
  if (_sav_nap_H != nap_H) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_rates =  - 200.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_rates)/400.; _mfac_rates = 1./_dx;
   for (_i=0, _x=_tmin_rates; _i < 401; _x += _dx, _i++) {
    _f_rates(_p, _ppvar, _thread, _nt, _x);
    _t_minf[_i] = minf;
    _t_mtau[_i] = mtau;
    _t_hinf[_i] = hinf;
    _t_htau[_i] = htau;
    _t_nap_inf[_i] = nap_inf;
    _t_nap_tau[_i] = nap_tau;
    _t_napi_inf[_i] = napi_inf;
    _t_napi_tau[_i] = napi_tau;
    _t_ninf[_i] = ninf;
    _t_ntau[_i] = ntau;
    _t_kif_a_inf[_i] = kif_a_inf;
    _t_kif_a_tau[_i] = kif_a_tau;
    _t_kif_i_inf[_i] = kif_i_inf;
    _t_kif_i_tau[_i] = kif_i_tau;
    _t_kis_a_inf[_i] = kis_a_inf;
    _t_kis_a_tau[_i] = kis_a_tau;
    _t_kis_i_inf[_i] = kis_i_inf;
    _t_kis_i_tau[_i] = kis_i_tau;
    _t_kh_m_inf[_i] = kh_m_inf;
    _t_kh_n_inf[_i] = kh_n_inf;
    _t_kh_m_tau[_i] = kh_m_tau;
    _t_kh_n_tau[_i] = kh_n_tau;
   }
   _sav_celsius = celsius;
   _sav_kif_avh = kif_avh;
   _sav_kif_ivh = kif_ivh;
   _sav_kis_avh = kis_avh;
   _sav_kis_ivh = kis_ivh;
   _sav_kd_avh = kd_avh;
   _sav_nap_A = nap_A;
   _sav_nap_C = nap_C;
   _sav_nap_D = nap_D;
   _sav_nap_E = nap_E;
   _sav_nap_G = nap_G;
   _sav_nap_H = nap_H;
  }
 }

 static rates(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv) { 
#if 0
_check_rates(_p, _ppvar, _thread, _nt);
#endif
 _n_rates(_p, _ppvar, _thread, _nt, _lv);
 return;
 }

 static _n_rates(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_rates(_p, _ppvar, _thread, _nt, _lv); return; 
}
 _xi = _mfac_rates * (_lv - _tmin_rates);
 _i = (int) _xi;
 if (_xi <= 0.) {
 minf = _t_minf[0];
 mtau = _t_mtau[0];
 hinf = _t_hinf[0];
 htau = _t_htau[0];
 nap_inf = _t_nap_inf[0];
 nap_tau = _t_nap_tau[0];
 napi_inf = _t_napi_inf[0];
 napi_tau = _t_napi_tau[0];
 ninf = _t_ninf[0];
 ntau = _t_ntau[0];
 kif_a_inf = _t_kif_a_inf[0];
 kif_a_tau = _t_kif_a_tau[0];
 kif_i_inf = _t_kif_i_inf[0];
 kif_i_tau = _t_kif_i_tau[0];
 kis_a_inf = _t_kis_a_inf[0];
 kis_a_tau = _t_kis_a_tau[0];
 kis_i_inf = _t_kis_i_inf[0];
 kis_i_tau = _t_kis_i_tau[0];
 kh_m_inf = _t_kh_m_inf[0];
 kh_n_inf = _t_kh_n_inf[0];
 kh_m_tau = _t_kh_m_tau[0];
 kh_n_tau = _t_kh_n_tau[0];
 return; }
 if (_i >= 400) {
 minf = _t_minf[400];
 mtau = _t_mtau[400];
 hinf = _t_hinf[400];
 htau = _t_htau[400];
 nap_inf = _t_nap_inf[400];
 nap_tau = _t_nap_tau[400];
 napi_inf = _t_napi_inf[400];
 napi_tau = _t_napi_tau[400];
 ninf = _t_ninf[400];
 ntau = _t_ntau[400];
 kif_a_inf = _t_kif_a_inf[400];
 kif_a_tau = _t_kif_a_tau[400];
 kif_i_inf = _t_kif_i_inf[400];
 kif_i_tau = _t_kif_i_tau[400];
 kis_a_inf = _t_kis_a_inf[400];
 kis_a_tau = _t_kis_a_tau[400];
 kis_i_inf = _t_kis_i_inf[400];
 kis_i_tau = _t_kis_i_tau[400];
 kh_m_inf = _t_kh_m_inf[400];
 kh_n_inf = _t_kh_n_inf[400];
 kh_m_tau = _t_kh_m_tau[400];
 kh_n_tau = _t_kh_n_tau[400];
 return; }
 _theta = _xi - (double)_i;
 minf = _t_minf[_i] + _theta*(_t_minf[_i+1] - _t_minf[_i]);
 mtau = _t_mtau[_i] + _theta*(_t_mtau[_i+1] - _t_mtau[_i]);
 hinf = _t_hinf[_i] + _theta*(_t_hinf[_i+1] - _t_hinf[_i]);
 htau = _t_htau[_i] + _theta*(_t_htau[_i+1] - _t_htau[_i]);
 nap_inf = _t_nap_inf[_i] + _theta*(_t_nap_inf[_i+1] - _t_nap_inf[_i]);
 nap_tau = _t_nap_tau[_i] + _theta*(_t_nap_tau[_i+1] - _t_nap_tau[_i]);
 napi_inf = _t_napi_inf[_i] + _theta*(_t_napi_inf[_i+1] - _t_napi_inf[_i]);
 napi_tau = _t_napi_tau[_i] + _theta*(_t_napi_tau[_i+1] - _t_napi_tau[_i]);
 ninf = _t_ninf[_i] + _theta*(_t_ninf[_i+1] - _t_ninf[_i]);
 ntau = _t_ntau[_i] + _theta*(_t_ntau[_i+1] - _t_ntau[_i]);
 kif_a_inf = _t_kif_a_inf[_i] + _theta*(_t_kif_a_inf[_i+1] - _t_kif_a_inf[_i]);
 kif_a_tau = _t_kif_a_tau[_i] + _theta*(_t_kif_a_tau[_i+1] - _t_kif_a_tau[_i]);
 kif_i_inf = _t_kif_i_inf[_i] + _theta*(_t_kif_i_inf[_i+1] - _t_kif_i_inf[_i]);
 kif_i_tau = _t_kif_i_tau[_i] + _theta*(_t_kif_i_tau[_i+1] - _t_kif_i_tau[_i]);
 kis_a_inf = _t_kis_a_inf[_i] + _theta*(_t_kis_a_inf[_i+1] - _t_kis_a_inf[_i]);
 kis_a_tau = _t_kis_a_tau[_i] + _theta*(_t_kis_a_tau[_i+1] - _t_kis_a_tau[_i]);
 kis_i_inf = _t_kis_i_inf[_i] + _theta*(_t_kis_i_inf[_i+1] - _t_kis_i_inf[_i]);
 kis_i_tau = _t_kis_i_tau[_i] + _theta*(_t_kis_i_tau[_i+1] - _t_kis_i_tau[_i]);
 kh_m_inf = _t_kh_m_inf[_i] + _theta*(_t_kh_m_inf[_i+1] - _t_kh_m_inf[_i]);
 kh_n_inf = _t_kh_n_inf[_i] + _theta*(_t_kh_n_inf[_i+1] - _t_kh_n_inf[_i]);
 kh_m_tau = _t_kh_m_tau[_i] + _theta*(_t_kh_m_tau[_i+1] - _t_kh_m_tau[_i]);
 kh_n_tau = _t_kh_n_tau[_i] + _theta*(_t_kh_n_tau[_i+1] - _t_kh_n_tau[_i]);
 }

 
static int  _f_rates ( _p, _ppvar, _thread, _nt, _lv ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lv ;
 {
   double _lalpha , _lbeta , _lsum ;
  _zq10 = pow( 3.0 , ( ( celsius - 22.0 ) / 10.0 ) ) ;
   minf = na_m ( _threadargscomma_ _lv ) ;
   mtau = na_mt ( _threadargscomma_ _lv ) ;
   hinf = na_h ( _threadargscomma_ _lv ) ;
   htau = na_ht ( _threadargscomma_ _lv ) ;
   nap_inf = na_p ( _threadargscomma_ _lv + nap_shift ) ;
   nap_tau = na_pt ( _threadargscomma_ _lv + nap_shift ) ;
   napi_inf = na_pi ( _threadargscomma_ _lv + nap_shift ) ;
   napi_tau = na_pit ( _threadargscomma_ _lv + nap_shift ) ;
   ninf = kd_m ( _threadargscomma_ _lv ) ;
   kif_a_inf = kif_m ( _threadargscomma_ _lv ) ;
   kif_i_inf = kif_h ( _threadargscomma_ _lv ) ;
   kif_a_tau = kif_mt ( _threadargscomma_ _lv ) ;
   kif_i_tau = kif_ht ( _threadargscomma_ _lv ) ;
   kis_a_inf = kis_m ( _threadargscomma_ _lv ) ;
   kis_i_inf = kis_h ( _threadargscomma_ _lv ) ;
   kis_a_tau = kis_mt ( _threadargscomma_ _lv ) ;
   kis_i_tau = kis_ht ( _threadargscomma_ _lv ) ;
   kh_m_inf = kh_m ( _threadargscomma_ _lv ) ;
   kh_n_inf = kh_n ( _threadargscomma_ _lv ) ;
   kh_m_tau = kh_mt ( _threadargscomma_ _lv ) ;
   kh_n_tau = kh_nt ( _threadargscomma_ _lv ) ;
    return 0; }
 
static int _hoc_rates() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 
#if 1
 _check_rates(_p, _ppvar, _thread, _nt);
#endif
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double na_m ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lna_m;
 _lna_m = 1.0 / ( 1.0 + exp ( - ( _lx + 38.0 ) / 3.0 ) ) ;
   
return _lna_m;
 }
 
static int _hoc_na_m() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  na_m ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double na_mt ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lna_mt;
 _lna_mt = mtau0 ;
   
return _lna_mt;
 }
 
static int _hoc_na_mt() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  na_mt ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double na_h ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lna_h;
 _lna_h = 1.0 / ( 1.0 + exp ( ( _lx + 43.0 ) / 3.0 ) ) ;
   
return _lna_h;
 }
 
static int _hoc_na_h() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  na_h ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double na_ht ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lna_ht;
 _lna_ht = htau0 ;
   
return _lna_ht;
 }
 
static int _hoc_na_ht() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  na_ht ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double na_p ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lna_p;
 _lna_p = alphbet ( _threadargscomma_ _lx , nap_A , nap_B , nap_C , nap_D ) ;
   _lna_p = _lna_p / ( _lna_p + alphbet ( _threadargscomma_ _lx , nap_E , nap_F , nap_G , nap_H ) ) ;
   
return _lna_p;
 }
 
static int _hoc_na_p() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  na_p ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double na_pt ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lna_pt;
 _lna_pt = alphbet ( _threadargscomma_ _lx , nap_A , nap_B , nap_C , nap_D ) ;
   _lna_pt = 1.0 / ( _lna_pt + alphbet ( _threadargscomma_ _lx , nap_E , nap_F , nap_G , nap_H ) ) ;
   
return _lna_pt;
 }
 
static int _hoc_na_pt() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  na_pt ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double na_pi ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lna_pi;
 _lna_pi = 0.06435 / ( 1.0 + exp ( ( _lx + 73.26415 ) / 3.71928 ) ) ;
   _lna_pi = _lna_pi / ( _lna_pi + ( 0.13496 / ( 1.0 + exp ( ( v + 10.27853 ) / ( - 9.09334 ) ) ) ) ) ;
   
return _lna_pi;
 }
 
static int _hoc_na_pi() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  na_pi ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double na_pit ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lna_pit;
 _lna_pit = 0.06435 / ( 1.0 + exp ( ( _lx + 73.26415 ) / 3.71928 ) ) ;
   _lna_pit = 1.0 / ( _lna_pit + ( 0.13496 / ( 1.0 + exp ( ( v + 10.27853 ) / ( - 9.09334 ) ) ) ) ) ;
   
return _lna_pit;
 }
 
static int _hoc_na_pit() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  na_pit ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double kd_m ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lkd_m;
 _lkd_m = 1.0 / ( 1.0 + exp ( - ( _lx - kd_avh ) / 3.0 ) ) ;
   
return _lkd_m;
 }
 
static int _hoc_kd_m() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  kd_m ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double kif_m ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lkif_m;
 _lkif_m = 1.0 / ( 1.0 + exp ( - ( _lx - kif_avh ) / 25.8 ) ) ;
   
return _lkif_m;
 }
 
static int _hoc_kif_m() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  kif_m ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double kif_h ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lkif_h;
 _lkif_h = 1.0 / ( 1.0 + exp ( ( _lx - kif_ivh ) / 6.7 ) ) ;
   
return _lkif_h;
 }
 
static int _hoc_kif_h() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  kif_h ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double kif_mt ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lkif_mt;
 _lkif_mt = 0.15 * exp ( ( _lx - kif_avh ) / 10.0 ) + 0.3 * exp ( - ( _lx - kif_avh ) / 10.0 ) ;
   _lkif_mt = 0.5 + ( 1.0 / _lkif_mt ) ;
   
return _lkif_mt;
 }
 
static int _hoc_kif_mt() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  kif_mt ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double kif_ht ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lkif_ht;
 _lkif_ht = 0.015 * exp ( ( _lx - kif_hivh ) / 20.0 ) + 0.03 * exp ( - ( _lx - kif_hivh ) / 20.0 ) ;
   _lkif_ht = 10.0 + ( 1.0 / _lkif_ht ) ;
   
return _lkif_ht;
 }
 
static int _hoc_kif_ht() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  kif_ht ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double kis_m ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lkis_m;
 _lkis_m = 1.0 / ( 1.0 + exp ( - ( _lx - kis_avh ) / 23.7 ) ) ;
   
return _lkis_m;
 }
 
static int _hoc_kis_m() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  kis_m ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double kis_h ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lkis_h;
 _lkis_h = 1.0 / ( 1.0 + exp ( ( _lx - kis_ivh ) / 9.0 ) ) ;
   
return _lkis_h;
 }
 
static int _hoc_kis_h() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  kis_h ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double kis_mt ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lkis_mt;
 _lkis_mt = 0.15 * exp ( ( _lx - kis_avh ) / 10.0 ) + 0.3 * exp ( - ( _lx - kis_avh ) / 10.0 ) ;
   _lkis_mt = 0.5 + 1.0 / _lkis_mt ;
   
return _lkis_mt;
 }
 
static int _hoc_kis_mt() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  kis_mt ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double kis_ht ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lkis_ht;
 _lkis_ht = 200.0 ;
   
return _lkis_ht;
 }
 
static int _hoc_kis_ht() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  kis_ht ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double kh_m ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lkh_m;
 _lkh_m = 1.0 / ( 1.0 + exp ( ( _lx + 68.9 + ghvshift ) / 6.5 ) ) ;
   
return _lkh_m;
 }
 
static int _hoc_kh_m() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  kh_m ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double kh_n ( _p, _ppvar, _thread, _nt, _lx ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx ;
 {
   double _lkh_n;
 _lkh_n = 1.0 / ( 1.0 + exp ( ( _lx + 68.9 + ghvshift ) / 6.5 ) ) ;
   
return _lkh_n;
 }
 
static int _hoc_kh_n() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  kh_n ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double kh_mt ( _p, _ppvar, _thread, _nt, _lv ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lv ;
 {
   double _lkh_mt;
 _lkh_mt = exp ( ( _lv + 183.6 + ghvshift ) / 15.24 ) ;
   
return _lkh_mt;
 }
 
static int _hoc_kh_mt() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  kh_mt ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double kh_nt ( _p, _ppvar, _thread, _nt, _lv ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lv ;
 {
   double _lkh_nt;
 _lkh_nt = exp ( ( _lv + 158.6 + ghvshift ) / 11.2 ) / ( 1.0 + exp ( ( _lv + 75.0 + ghvshift ) / 5.5 ) ) ;
   
return _lkh_nt;
 }
 
static int _hoc_kh_nt() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  kh_nt ( _p, _ppvar, _thread, _nt, *getarg(1) ) ;
 ret(_r);
}
 
double vtrap ( _p, _ppvar, _thread, _nt, _lx , _ly ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
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
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  vtrap ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) ) ;
 ret(_r);
}
 
double alphbet ( _p, _ppvar, _thread, _nt, _lx , _lA , _lB , _lC , _lD ) double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt; 
	double _lx , _lA , _lB , _lC , _lD ;
 {
   double _lalphbet;
 _lalphbet = _lA / ( _lB + exp ( ( _lx + _lC ) / _lD ) ) ;
   
return _lalphbet;
 }
 
static int _hoc_alphbet() {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  alphbet ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) , *getarg(5) ) ;
 ret(_r);
}
 
static int _ode_count(_type) int _type;{ return 11;}
 
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
  ek = _ion_ek;
  eh = _ion_eh;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
    }}
 
static int _ode_map(_ieq, _pv, _pvdot, _pp, _ppd, _atol, _type) int _ieq, _type; double** _pv, **_pvdot, *_pp, *_atol; Datum* _ppd; { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 11; ++_i) {
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
  ek = _ion_ek;
  eh = _ion_eh;
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _thread_mem_init(Datum* _thread) {
   _thread[0]._pval = (double*)ecalloc(12, sizeof(double));
 }
 
static void _thread_cleanup(Datum* _thread) {
   free((void*)(_thread[0]._pval));
 }
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_k_sym, _ppvar, 3, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 4, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 5, 4);
   nrn_update_ion_pointer(_h_sym, _ppvar, 6, 0);
   nrn_update_ion_pointer(_h_sym, _ppvar, 7, 3);
   nrn_update_ion_pointer(_h_sym, _ppvar, 8, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  h = h0;
  khn = khn0;
  khm = khm0;
  kisi = kisi0;
  kisa = kisa0;
  kifi = kifi0;
  kifa = kifa0;
  m = m0;
  napi = napi0;
  nap = nap0;
  n = n0;
 {
   rates ( _threadargscomma_ v ) ;
   m = minf ;
   h = hinf ;
   nap = nap_inf ;
   napi = napi_inf ;
   n = ninf ;
   khm = kh_m_inf ;
   khn = kh_n_inf ;
   if ( kif_a_start < 0.0 ) {
     kifa = kif_a_inf ;
     }
   else {
     kifa = kif_a_start ;
     }
   if ( kif_i_start < 0.0 ) {
     kifi = kif_i_inf ;
     }
   else {
     kifi = kif_i_start ;
     }
   if ( kis_a_start < 0.0 ) {
     kisa = kis_a_inf ;
     }
   else {
     kisa = kis_a_start ;
     }
   if ( kis_i_start < 0.0 ) {
     kisi = kis_i_inf ;
     }
   else {
     kisi = kis_i_start ;
     }
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
 _check_rates(_p, _ppvar, _thread, _nt);
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
  ek = _ion_ek;
  eh = _ion_eh;
 initmodel(_p, _ppvar, _thread, _nt);
   }}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   gna = gnabar * m * m * h ;
   gnap = gnapbar * nap * nap * nap * napi ;
   ina = gna * ( v - ena ) + gnap * ( v - ena ) + 0.25 * gh * ( v - ena ) ;
   gk = gkbar * n * n ;
   akif = kifa * kifa * kifa * kifa * kifi ;
   gkif = gkifbar * akif ;
   akis = kisa * kisa * kisa * kisa * kisi ;
   gkis = gkisbar * akis ;
   aih = khm * khn ;
   gh = ghbar * aih ;
   ik = gk * ( v - ek ) + gkif * ( v - ek ) + gkis * ( v - ek ) + 0.75 * gh * ( v - ek ) ;
   ih = gh * ( v - eh ) ;
   il = gl * ( v - el ) ;
   }
 _current += ina;
 _current += ik;
 _current += ih;
 _current += il;

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
  ek = _ion_ek;
  eh = _ion_eh;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dih;
 double _dik;
 double _dina;
  _dina = ina;
  _dik = ik;
  _dih = ih;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dinadv += (_dina - ina)/.001 ;
  _ion_dikdv += (_dik - ik)/.001 ;
  _ion_dihdv += (_dih - ih)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
  _ion_ik += ik ;
  _ion_ih += ih ;
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
  ek = _ion_ek;
  eh = _ion_eh;
 { {
 for (; t < _break; t += dt) {
   states(_p, _ppvar, _thread, _nt);
  
}}
 t = _save;
 }   }}

}

static terminal(){}

static _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
 _slist1[2] = &(nap) - _p;  _dlist1[2] = &(Dnap) - _p;
 _slist1[3] = &(napi) - _p;  _dlist1[3] = &(Dnapi) - _p;
 _slist1[4] = &(n) - _p;  _dlist1[4] = &(Dn) - _p;
 _slist1[5] = &(kifa) - _p;  _dlist1[5] = &(Dkifa) - _p;
 _slist1[6] = &(kifi) - _p;  _dlist1[6] = &(Dkifi) - _p;
 _slist1[7] = &(kisa) - _p;  _dlist1[7] = &(Dkisa) - _p;
 _slist1[8] = &(kisi) - _p;  _dlist1[8] = &(Dkisi) - _p;
 _slist1[9] = &(khm) - _p;  _dlist1[9] = &(Dkhm) - _p;
 _slist1[10] = &(khn) - _p;  _dlist1[10] = &(Dkhn) - _p;
   _t_minf = makevector(401*sizeof(double));
   _t_mtau = makevector(401*sizeof(double));
   _t_hinf = makevector(401*sizeof(double));
   _t_htau = makevector(401*sizeof(double));
   _t_nap_inf = makevector(401*sizeof(double));
   _t_nap_tau = makevector(401*sizeof(double));
   _t_napi_inf = makevector(401*sizeof(double));
   _t_napi_tau = makevector(401*sizeof(double));
   _t_ninf = makevector(401*sizeof(double));
   _t_ntau = makevector(401*sizeof(double));
   _t_kif_a_inf = makevector(401*sizeof(double));
   _t_kif_a_tau = makevector(401*sizeof(double));
   _t_kif_i_inf = makevector(401*sizeof(double));
   _t_kif_i_tau = makevector(401*sizeof(double));
   _t_kis_a_inf = makevector(401*sizeof(double));
   _t_kis_a_tau = makevector(401*sizeof(double));
   _t_kis_i_inf = makevector(401*sizeof(double));
   _t_kis_i_tau = makevector(401*sizeof(double));
   _t_kh_m_inf = makevector(401*sizeof(double));
   _t_kh_n_inf = makevector(401*sizeof(double));
   _t_kh_m_tau = makevector(401*sizeof(double));
   _t_kh_n_tau = makevector(401*sizeof(double));
_first = 0;
}
