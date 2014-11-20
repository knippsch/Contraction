/*******************************************************************************
*
* File measure.h
*
* Copyright (C) 2008-11 Andreas Juettner
*                       Bastian Knippschild     
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/
#define MAX_INFILE_SIZE 500000
#define NMAX_TWIST 100
#define NMAX_PARM 1000
#define NMAX_DIRAC 16

#ifdef SMEARLINKS
  #define EXTERN extern
#else
  #define EXTERN 
#endif

EXTERN su3_dble *pud_sm1[VOLUME][4];
EXTERN su3_dble *pud_sm2[VOLUME][4];

typedef struct {
 int p[3],musrc,musnk __attribute__((__packed__));
 complex_dble corr[NPROC0*L0] __attribute__((__packed__));
} twopt_mom;

typedef struct {
  su3_dble c11,c12,c13,c14,c21,c22,c23,c24,c31,c32,c33,c34,
	c41,c42,c43,c44;
} full_spinor_dble;

typedef struct {
  int smearid;
  char type[7];
  double alpha[3];
} smearparm;
typedef struct {
  char type[7];
  int  pos[4];
  double  kappa;
  int  n;
  smearparm smear;
} srcdef;
typedef struct {
  double kappa;
  double csw;
  srcdef src;
  double beta; 
  double theta[4]; 
  int    nmx;
  double res;
  int    t_snk;
  int    *mom_ins;
  int    gam_ins;
} propinfo;

#ifndef RNDGAUGE_C
  extern void transform_sd(spinor_dble *pk, spinor_dble *pl);
  extern void transform_ud(su3_dble *input[VOLUME][4]);
#endif
#ifndef PROPAGATORS_C
  extern void propagator(full_spinor_dble *sv, spinor_dble *rnd, propinfo prop);
  extern void seq_propagator(full_spinor_dble *sv, full_spinor_dble *src, 
 		propinfo prop);
  extern void seq_propagator_b(full_spinor_dble *sv_out, full_spinor_dble *svk,
			       full_spinor_dble *svl,    full_spinor_dble *svm, 
			       propinfo prop, int snk_sm, int pol);
  extern void propagator_col(spinor_dble *sv, spinor_dble *rnd, propinfo prop,
			int idirac, int icolor);
#endif
#ifndef CORRELATORS_C
  extern void C2(full_spinor_dble *k, full_spinor_dble *l, int num_c, int *mu, 
		               int nmom_max, int *momenta, complex_dble *sst);
  extern void C2a2a(int nhits,spinor_dble **k, spinor_dble **kpsi, 
 	                  spinor_dble **l, spinor_dble **lpsi,int num_c, int *mu, 
 	                  int *momenta, complex_dble *sst);
  extern void C2_b(full_spinor_dble *svk, full_spinor_dble *svl, full_spinor_dble *svm, 
	                  int numb,int *mub, int nmom_max,int *momenta, complex_dble **s_pointer, int pol);
  extern void C3_b(full_spinor_dble *sv0, full_spinor_dble *sv1, 
                  	int numb, int *mub, int nmom_max, int *momenta, 
                  	complex_dble *s_pointer, propinfo prop);
  extern void seq_src_test(full_spinor_dble *k, int gamma_insertion, int *pos);
  extern int gen_momenta(int comp_max,int *momenta);
#endif

#ifndef ALGEBRA_C
  extern void copy_sv_fs(int vol,spinor_dble *s,full_spinor_dble *sv,int id,int ic);
  extern void copy_fs_sv(int vol,full_spinor_dble *sv,spinor_dble *s,int id,int ic);
  extern void set_fsd2zero(int vol,full_spinor_dble *pk);
  extern void mul_cmplx_fsv_dble(complex_dble *phase,full_spinor_dble *in,full_spinor_dble *out);
  extern void full_mat_prod(full_spinor_dble *in1,full_spinor_dble *in2,full_spinor_dble *out);
  extern su3_vector_dble mul_cplx(complex_dble z,su3_vector_dble s);
  extern spinor_dble mul_cplx_spinor_dble(complex_dble z,spinor_dble s);
  extern void mul_cplx_su3_dble(complex_dble *z,su3_dble *in, su3_dble *out);
  extern void mul_gamma_l(int mu,full_spinor_dble *in, full_spinor_dble *out);
  extern void mul_gamma_r(int mu,full_spinor_dble *in, full_spinor_dble *out);
  extern void meson_trace(full_spinor_dble *in1,full_spinor_dble *in2,complex_dble *out);
  extern void adj_full_spinor(full_spinor_dble *in, full_spinor_dble *out);
  extern complex_dble spinor_prod_dble_local(spinor_dble *k,spinor_dble *l);
  extern void mul_s_gamma_l(int mu,spinor_dble *in,spinor_dble *out);
  extern void mul_s_gamma_r(int mu,spinor_dble *in,spinor_dble *out);
  extern void accumulate_fs(full_spinor_dble *sk, int nhits,spinor_dble *k,spinor_dble *l);
  extern void fsv_trace(full_spinor_dble *in,complex_dble *out);
  extern full_spinor_dble su3_times_fsd(su3_dble u, full_spinor_dble s);
  extern full_spinor_dble fsd_times_su3(su3_dble u, full_spinor_dble s);
  extern su3_dble mul_cplx_su3_dble_out(complex_dble z,su3_dble s);
  extern full_spinor_dble mul_one_pm_gamma_l(double sign, int mu,full_spinor_dble s);
#endif

#ifndef SOURCE_C
  typedef struct {int x; int y; int z; int t;} site;

  site coords[VOLUME];

  void init_coords();
  extern void random_Z4(int vol,spinor_dble *pk);
  extern void SourceRadius(srcdef src,spinor_dble *sd);
  extern void srcfld(spinor_dble *rnd, spinor_dble *sd, srcdef src, int idirac,
 		                  int icolor);
  extern void src_field_random(int vol, spinor_dble *pk, int t);
  extern void src_field_dt_random(int vol, spinor_dble *pk, int d, int t);
  extern void extended_src(int igamma, int *mom, int x0,
 		                        full_spinor_dble *k, full_spinor_dble *l);
  extern void jacobi_source(spinor_dble *sd, srcdef src, int idirac, int icolor);
  extern void jacobi_sink(spinor_dble *sd, srcdef src, int idirac, int icolor);
#endif

#ifndef LINK_SMEARING_C
  extern void approx_project_to_su3_dble(su3_dble *u, int iter);
  extern void copy_bnd_ud_reverse(su3_dble *u1[VOLUME][4]);
  extern void free_ucom_bufs_reverse(void);
  extern void alloc_pud_sm1(void);
  extern void alloc_pud_sm2(void);
  extern void free_pud_sm1(void);
  extern void free_pud_sm2(void);
  extern void su3_dble_to_su3(su3_dble *um1, su3 *um2);
  extern void su3_to_su3_dble(su3_dble *um1, su3 *um2);
  extern void pud_output(su3_dble *doublelinks[VOLUME][4]);
  extern void pu_output(su3 *singlelinks[VOLUME][4]);
  extern void cmp_2pud(su3_dble *doublelinks1[VOLUME][4], 
		                     su3_dble *doublelinks2[VOLUME][4]);
  extern void pud_copy(su3_dble *u1[VOLUME][4], su3_dble *u2[VOLUME][4]);
  extern void APE_smearing(su3_dble *u1[VOLUME][4], su3_dble *u2[VOLUME][4], 
		                         smearparm parm);
  extern void HYP_smearing(su3_dble *u1[VOLUME][4], su3_dble *u2[VOLUME][4], 
		                         smearparm parm);
  extern void check_SU3(su3_dble *u1[VOLUME][4]);
#endif


#ifndef TWISTBC_C
  extern void twistbc(double *theta,double sign);
#endif

#ifndef ARCH_ILDG_C
  extern void import_cnfg_ildg(char *in);
#endif

#ifndef IOARCHIVE_C
  extern long file_size(FILE *fp);
  extern void export_init(char *out_file, int argc,char *argv[]);
  extern void meson_IOloop(char *tagname,char *out_file, complex_dble *corr, 
 			                       int nmom_max, int *mom, int *mu, int num_c);
  extern void baryon_IOloop(char *tagname,char *out_file, complex_dble *corr,  			
                            int nmom_max, int *mom, int *mu, int num_c);
  extern void meson_ASCII_IO(char *filename,int rec_seek, int msg_seek,int nmom_max);
  extern void read_fsv(char *in, full_spinor_dble *sv);
  extern void write_fsv(char *in, full_spinor_dble *sv);
  extern void read_fsv2(char *stem, int prop_dump,int *pim, int nprop,int *isv,
	                       full_spinor_dble **prop_pt,full_spinor_dble **sv,full_spinor_dble **wsd);
  
#endif

#ifndef BARYONS_C
  extern void C3_check(full_spinor_dble *sv0, propinfo prop); 
  extern void cc_gamma_row(int mu, spinor_dble  s[3][4], spinor_dble r[3][4]);
  extern void cc_gamma_column(int mu, spinor_dble s[3][4], spinor_dble r[3][4]);
  extern void twopt_baryons_dble(int gam, int gamp, int num_c, int *mu, int nmom_max, int *momenta, 
				                             full_spinor_dble *svk, full_spinor_dble *svl, double pm, 
				                             complex_dble *uds, complex_dble *uud, complex_dble *uuu); 
  extern void twopt_omega_dble(int gam, int gamp, int num_c, int *mu, int num_mom,
                               int *momenta, full_spinor_dble *svk, double pm,
                               complex_dble *sst);
  extern void d_source_dble(int t, int gam, int gamp, full_spinor_dble *svk, double pm, 
				                        spinor_dble *grand_source);
  extern void u_source_dble(int t, int gam, int gamp, full_spinor_dble *svk, full_spinor_dble *svl, 
				                        double pm, spinor_dble *grand_source);
  extern void check_closure_dble(spinor_dble *grand_source, full_spinor_dble *propagator, 
                                  full_spinor_dble *ext_prop, propinfo prop);
#endif



#define _trace_full_spinor_re(s,r)\
  s=\
  r.c11.c11.re+r.c11.c22.re+r.c11.c33.re+\
  r.c12.c11.re+r.c12.c22.re+r.c12.c33.re+\
  r.c13.c11.re+r.c13.c22.re+r.c13.c33.re+\
  r.c14.c11.re+r.c14.c22.re+r.c14.c33.re+\
  r.c21.c11.re+r.c21.c22.re+r.c21.c33.re+\
  r.c22.c11.re+r.c22.c22.re+r.c22.c33.re+\
  r.c23.c11.re+r.c23.c22.re+r.c23.c33.re+\
  r.c24.c11.re+r.c24.c22.re+r.c24.c33.re+\
  r.c31.c11.re+r.c31.c22.re+r.c31.c33.re+\
  r.c32.c11.re+r.c32.c22.re+r.c32.c33.re+\
  r.c33.c11.re+r.c33.c22.re+r.c33.c33.re+\
  r.c34.c11.re+r.c34.c22.re+r.c34.c33.re+\
  r.c41.c11.re+r.c41.c22.re+r.c41.c33.re+\
  r.c42.c11.re+r.c42.c22.re+r.c42.c33.re+\
  r.c43.c11.re+r.c43.c22.re+r.c43.c33.re+\
  r.c44.c11.re+r.c44.c22.re+r.c44.c33.re;

#define _trace_full_spinor_im(s,r)\
  s=\
  r.c11.c11.im+r.c11.c22.im+r.c11.c33.im+\
  r.c12.c11.im+r.c12.c22.im+r.c12.c33.im+\
  r.c13.c11.im+r.c13.c22.im+r.c13.c33.im+\
  r.c14.c11.im+r.c14.c22.im+r.c14.c33.im+\
  r.c21.c11.im+r.c21.c22.im+r.c21.c33.im+\
  r.c22.c11.im+r.c22.c22.im+r.c22.c33.im+\
  r.c23.c11.im+r.c23.c22.im+r.c23.c33.im+\
  r.c24.c11.im+r.c24.c22.im+r.c24.c33.im+\
  r.c31.c11.im+r.c31.c22.im+r.c31.c33.im+\
  r.c32.c11.im+r.c32.c22.im+r.c32.c33.im+\
  r.c33.c11.im+r.c33.c22.im+r.c33.c33.im+\
  r.c34.c11.im+r.c34.c22.im+r.c34.c33.im+\
  r.c41.c11.im+r.c41.c22.im+r.c41.c33.im+\
  r.c42.c11.im+r.c42.c22.im+r.c42.c33.im+\
  r.c43.c11.im+r.c43.c22.im+r.c43.c33.im+\
  r.c44.c11.im+r.c44.c22.im+r.c44.c33.im;

/*
* adds two su3 matrices
*/
#define _su3_add(r, u, s) \
(r).c11.re = (u).c11.re + (s).c11.re; \
(r).c11.im = (u).c11.im + (s).c11.im; \
(r).c12.re = (u).c12.re + (s).c12.re; \
(r).c12.im = (u).c12.im + (s).c12.im; \
(r).c13.re = (u).c13.re + (s).c13.re; \
(r).c13.im = (u).c13.im + (s).c13.im; \
(r).c21.re = (u).c21.re + (s).c21.re; \
(r).c21.im = (u).c21.im + (s).c21.im; \
(r).c22.re = (u).c22.re + (s).c22.re; \
(r).c22.im = (u).c22.im + (s).c22.im; \
(r).c23.re = (u).c23.re + (s).c23.re; \
(r).c23.im = (u).c23.im + (s).c23.im; \
(r).c31.re = (u).c31.re + (s).c31.re; \
(r).c31.im = (u).c31.im + (s).c31.im; \
(r).c32.re = (u).c32.re + (s).c32.re; \
(r).c32.im = (u).c32.im + (s).c32.im; \
(r).c33.re = (u).c33.re + (s).c33.re; \
(r).c33.im = (u).c33.im + (s).c33.im



