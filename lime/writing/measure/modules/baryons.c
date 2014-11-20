/******************************************************************************
 *
 *  File baryons.c
 *
 *  Computation of 2-pt correlation functions of the baryons 
 *  (octect and decuplet).
 *
 *  Computes also the 3-point correlators of the nucleon.
 *
 *  More details and explicit formulae can be found in the notes:
 *  "Correlation functions of the baryons"
 *  by  S. Capitani
 *
 *
 *  The externally accessible functions are:
 *
 * 
 *    void cc_gamma_row(int mu, spinor_dble s[3][4], spinor_dble r[3][4])
 * 
 *      Multiplies the 12 spinors s with C * gamma_mu, where mu=0,3 for 
 *      the standard Dirac matrices, mu=4 is the identity, and mu=5 is gamma_5
 *      These gammas operate on the ROW indices of the propagators
 *
 *      We use  C=gamma0*gamma2
 *
 *      Notice that the representation of the Dirac matrices must be the chiral
 *      one as in the notes... If this is changed, then this and the following
 *      routine must also be changed accordingly
 * 
 *
 *    void cc_gamma_column(int mu, spinor_dble s[3][4], spinor_dble r[3][4])
 *
 *      Like before, but: 
 *        - the gammas operate on the COLUMN indices of the propagators
 *        - the atual matrices are  gamma_0 times \Gamma^\* \gamma_0,
 *          where Gamma = C * gamma_mu
 *
 *
 *    void twopt_baryons_dble(int gam, int gamp, int num_c, int *mu, 
 *              int *momenta, spinor_dble const *grand_fprop, double pm, 
 *              complex_dble *uds, complex_dble *uud, complex_dble *uuu)
 *
 *      Returns the two-point correlator of baryons composed of:
 *        - three quarks of different flavor  [ *uds ]
 *        - two quarks of one flavor, a third quark of another flavor  [ *uud ]
 *             (in particular the nucleon, channel 5 5)
 *        - all three quarks of the same flavor  [ *uuu]
 *             (in particular the omega, channels 2 2 etc.)    
 *          NOTE: the case uuu does not work if a gamma_5 is present, in which
 *          case the function "twopt_omega_dble" should be used (see below)
 *      
 *      All three quarks are always degenerate in mass
 *
 *      gam and gamp are the spinor indices of the baryonic currents J 
 *        and bar{J}
 *
 *      The calculation is done as follows:
 *      - the Dirac matrices are worked out as soon as possible 
 *      - then the row Dirac indices of two propagators are summed over,
 *        and at the same time the cross product of their row color indices
 *        IS TAKEN 
 *      - then the rest of the Dirac indices are treated as needed in the
 *        particular contraction, while at the same time the cross product 
 *        of the column color indices is taken, and all color indices
 *        are summed over 
 *
 *      The general rule is: first the non-primed indices (row) are done, then 
 *        the primed indices (column); but the Dirac matrices are worked out
 *        as soon as possible
 *
 *
 *    void twopt_omega_dble(int gam, int gamp, int num_c, int *mu, int *momenta,
 *              spinor_dble const *grand_fprop, double pm, complex_dble *sst)
 *       
 *      This functions contains the 6 contractions for the case "uuu"
 *
 *      If the channels containing a gamma_5 are not needed, the function
 *      "twopt_baryons_dble" should be used, which is much faster as it
 *      contains only the 2 contractions to which in this case can be reduced
 *
 * 
 *    void d_source_dble(int t, int gam, int gamp, 
 *              spinor_dble const *grand_fprop, double pm, 
 *              spinor_dble *grand_source)
 *
 *      Generates the source to be used in the inversion which computes the 
 *      extended propagator, when the operator is put on a "d" quark line)
 *
 *      Returns "gamma_5 * the adjoint" of the source, which is the actual
 *      object that has to be used when the Dirac operator acts on its right
 *
 * 
 *    void u_source_dble(int t, int gam, int gamp, 
 *              spinor_dble const *grand_fprop, double pm, 
 *              spinor_dble *grand_source)
 *
 *      Same as above, but when the operator is put on a "u" quark line
 * 
 * 
 *    void check_closure_dble(spinor_dble *grand_source)
 * 
 *       If the calculation of the extended propagator is correct, 
 *       the numbers which are printed by this subroutine are equal
 *       to 3 times the calues of the 2-point correlator at the timeslice 
 *       where the sink is located
 *
 *
 *  Author:  Stefano Capitani  <capitan@kph.uni-mainz.de>
 *           Bastian Knippschild <knippsch@kph.uni-mainz.de>
 *
 *  Date:  November 2008-11
 *
 *****************************************************************************/

#define BARYONS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "global.h"
#include "su3.h"
#include "start.h"
#include "measure.h"
#include "mpi.h"
#include "flags.h"
#include "random.h"
#include "misc.h"
#include "sw_term.h"
#include "dfl.h"
#include "sap_gcr.h"
#include "version.h"
#include "version_qcd-measure.h"
#include "update.h"
#include "linalg.h"


#if (defined SSE3)
#include "sse2.h"

static const complex_dble cplxd0={0};

static sse_int _sse_sgn_dble_2 __attribute__ ((unused)) ={0x0,0x0,0x0,0x80000000};

#define _sse_load_3_dble(s) \
__asm__ __volatile__ ("movapd %0, %%xmm6 \n\t" \
                      "movapd %1, %%xmm7 \n\t" \
                      "movapd %2, %%xmm8" \
                      : \
                      : \
                      "m" ((s).c1), \
                      "m" ((s).c2), \
                      "m" ((s).c3) \
                      : \
                      "xmm6", "xmm7", "xmm8")

#define _sse_load_4_dble(s) \
__asm__ __volatile__ ("movapd %0, %%xmm9 \n\t" \
                      "movapd %1, %%xmm10 \n\t" \
                      "movapd %2, %%xmm11" \
                      : \
                      : \
                      "m" ((s).c1), \
                      "m" ((s).c2), \
                      "m" ((s).c3) \
                      : \
                      "xmm9", "xmm10", "xmm11")
/* calculates the cross product of the two su3-vectors */
/* vectors r and s. The result is stored in xmm0-2     */
#define _sse_cross_color_add(r,s) \
 __asm__ __volatile__ ("movddup %1, %%xmm3 \n\t" \
                       "mulpd %5, %%xmm3 \n\t" \
                       "movddup %3, %%xmm4 \n\t" \
                       "mulpd %4, %%xmm4 \n\n" \
                       "subpd %%xmm4, %%xmm3 \n\t" \
                       "shufpd $0x1, %%xmm3, %%xmm3 \n\t" \
                       "movddup %2, %%xmm4 \n\t" \
                       "mulpd %4, %%xmm4 \n\t" \
                       "movddup %0, %%xmm5 \n\t" \
                       "mulpd %5, %%xmm5 \n\t" \
                       "subpd %%xmm4, %%xmm5 \n\t" \
                       "addsubpd %%xmm3, %%xmm5 \n\t" \
                       "addpd %%xmm5, %%xmm2" \
                       : \
                       : \
                       "m" ((r).c1.re), \
                       "m" ((r).c1.im), \
                       "m" ((r).c2.re), \
                       "m" ((r).c2.im), \
                       "m" ((s).c1), \
                       "m" ((s).c2)); \
 __asm__ __volatile__ ("movddup %1, %%xmm3 \n\t" \
                       "mulpd %5, %%xmm3 \n\t" \
                       "movddup %3, %%xmm4 \n\t" \
                       "mulpd %4, %%xmm4 \n\n" \
                       "subpd %%xmm4, %%xmm3 \n\t" \
                       "shufpd $0x1, %%xmm3, %%xmm3 \n\t" \
                       "movddup %2, %%xmm4 \n\t" \
                       "mulpd %4, %%xmm4 \n\t" \
                       "movddup %0, %%xmm5 \n\t" \
                       "mulpd %5, %%xmm5 \n\t" \
                       "subpd %%xmm4, %%xmm5 \n\t" \
                       "addsubpd %%xmm3, %%xmm5 \n\t" \
                       "addpd %%xmm5, %%xmm1" \
                       : \
                       : \
                       "m" ((r).c3.re), \
                       "m" ((r).c3.im), \
                       "m" ((r).c1.re), \
                       "m" ((r).c1.im), \
                       "m" ((s).c3), \
                       "m" ((s).c1)); \
 __asm__ __volatile__ ("movddup %1, %%xmm3 \n\t" \
                       "mulpd %5, %%xmm3 \n\t" \
                       "movddup %3, %%xmm4 \n\t" \
                       "mulpd %4, %%xmm4 \n\n" \
                       "subpd %%xmm4, %%xmm3 \n\t" \
                       "shufpd $0x1, %%xmm3, %%xmm3 \n\t" \
                       "movddup %2, %%xmm4 \n\t" \
                       "mulpd %4, %%xmm4 \n\t" \
                       "movddup %0, %%xmm5 \n\t" \
                       "mulpd %5, %%xmm5 \n\t" \
                       "subpd %%xmm4, %%xmm5 \n\t" \
                       "addsubpd %%xmm3, %%xmm5 \n\t" \
                       "addpd %%xmm5, %%xmm0" \
                       : \
                       : \
                       "m" ((r).c2.re), \
                       "m" ((r).c2.im), \
                       "m" ((r).c3.re), \
                       "m" ((r).c3.im), \
                       "m" ((s).c2), \
                       "m" ((s).c3) \
                       : \
                       "xmm0", "xmm1", "xmm2", \
                       "xmm3", "xmm4", "xmm5")
#endif

static void flush_corr(int num_c,int num_mom, complex_dble *s)
{
  int i;
  for (i=0;i<num_c*num_c*num_mom;i++) {
    s[i].re=0.0;
    s[i].im=0.0;
    }
}
static void flush_corr2(int num_c,int num_mom, complex_dble *s)
{
  int i;
  for (i=0;i<num_c*num_mom;i++) {
    s[i].re=0.0;
    s[i].im=0.0;
    }
}

static void add_corr(int num_c, int num_mom, complex_dble *si, complex_dble *so){
  int i;
  for (i=0;i<num_c*num_c*num_mom;i++){
   so[i].re += si[i].re;
   so[i].im += si[i].im;
   }
}
static void add_corr2(int num_c, int num_mom, complex_dble *si, complex_dble *so){
  int i;
  for (i=0;i<num_c*num_mom;i++){
   so[i].re += si[i].re;
   so[i].im += si[i].im;
   }
}

static const spinor_dble sd0={{{0.0}}};
static const su3_vector_dble svd0={{0.0}};


#define _resprod_re(r,s) \
   (r).c1.re*(s).c1.re-(r).c1.im*(s).c1.im+ \
   (r).c2.re*(s).c2.re-(r).c2.im*(s).c2.im+ \
   (r).c3.re*(s).c3.re-(r).c3.im*(s).c3.im

#define _resprod_im(r,s) \
   (r).c1.re*(s).c1.im+(r).c1.im*(s).c1.re+ \
   (r).c2.re*(s).c2.im+(r).c2.im*(s).c2.re+ \
   (r).c3.re*(s).c3.im+(r).c3.im*(s).c3.re


#define _cross_color(v,w,z) \
   (v).c1.re= (w).c2.re*(z).c3.re-(w).c2.im*(z).c3.im  \
             -(w).c3.re*(z).c2.re+(w).c3.im*(z).c2.im; \
   (v).c1.im=-(w).c3.re*(z).c2.im-(w).c3.im*(z).c2.re  \
             +(w).c2.re*(z).c3.im+(w).c2.im*(z).c3.re; \
   (v).c2.re= (w).c3.re*(z).c1.re-(w).c3.im*(z).c1.im  \
             -(w).c1.re*(z).c3.re+(w).c1.im*(z).c3.im; \
   (v).c2.im=-(w).c1.re*(z).c3.im-(w).c1.im*(z).c3.re  \
             +(w).c3.re*(z).c1.im+(w).c3.im*(z).c1.re; \
   (v).c3.re= (w).c1.re*(z).c2.re-(w).c1.im*(z).c2.im  \
             -(w).c2.re*(z).c1.re+(w).c2.im*(z).c1.im; \
   (v).c3.im=-(w).c2.re*(z).c1.im-(w).c2.im*(z).c1.re  \
             +(w).c1.re*(z).c2.im+(w).c1.im*(z).c2.re


void cc_gamma_row(int mu, spinor_dble s[3][4], spinor_dble r[3][4])
{
	complex_dble i,m_i,m_1,p_1;
	int id,ic;
	
	i.re=0.0;
	i.im=1.0;
	
	m_i.re=0.0;
	m_i.im=-1.0;
	
	m_1.re=-1.0;
	m_1.im=0.0;
	
	p_1.re=1.0;
	p_1.im=0.0;
	
	if (mu==0)
	{
		for (id=0;id<4;id++)
		{
			for (ic=0;ic<3;ic++)
			{
				_vector_mulc(r[ic][id].c1,p_1,s[ic][id].c4);
				_vector_mulc(r[ic][id].c2,m_1,s[ic][id].c3);
				_vector_mulc(r[ic][id].c3,m_1,s[ic][id].c2);
				_vector_mulc(r[ic][id].c4,p_1,s[ic][id].c1);
			}
		}
	}
	else if (mu==1)
	{
		for (id=0;id<4;id++)
		{
			for (ic=0;ic<3;ic++)
			{
				_vector_mulc(r[ic][id].c1,i,s[ic][id].c3);
				_vector_mulc(r[ic][id].c2,m_i,s[ic][id].c4);
				_vector_mulc(r[ic][id].c3,i,s[ic][id].c1);
				_vector_mulc(r[ic][id].c4,m_i,s[ic][id].c2);
			}
		}
	}
	else if (mu==2)
	{
		for (id=0;id<4;id++)
		{
			for (ic=0;ic<3;ic++)
			{
				_vector_mulc(r[ic][id].c1,m_1,s[ic][id].c3);
				_vector_mulc(r[ic][id].c2,m_1,s[ic][id].c4);
				_vector_mulc(r[ic][id].c3,m_1,s[ic][id].c1);
				_vector_mulc(r[ic][id].c4,m_1,s[ic][id].c2);
			}
		}
	}
	else if (mu==3)
	{
		for (id=0;id<4;id++)
		{
			for (ic=0;ic<3;ic++)
			{
				_vector_mulc(r[ic][id].c1,m_i,s[ic][id].c4);
				_vector_mulc(r[ic][id].c2,m_i,s[ic][id].c3);
				_vector_mulc(r[ic][id].c3,m_i,s[ic][id].c2);
				_vector_mulc(r[ic][id].c4,m_i,s[ic][id].c1);
			}
		}
	}
	else if (mu==4)
	{
		for (id=0;id<4;id++)
		{
			for (ic=0;ic<3;ic++)
			{
				_vector_mulc(r[ic][id].c1,m_i,s[ic][id].c2);
				_vector_mulc(r[ic][id].c2,m_i,s[ic][id].c1);
				_vector_mulc(r[ic][id].c3,i,s[ic][id].c4);
				_vector_mulc(r[ic][id].c4,i,s[ic][id].c3);
			}
		}
	}
	else if (mu==5)
	{
		for (id=0;id<4;id++)
		{
			for (ic=0;ic<3;ic++)
			{
				_vector_mulc(r[ic][id].c1,m_1,s[ic][id].c2);
				_vector_mulc(r[ic][id].c2,p_1,s[ic][id].c1);
				_vector_mulc(r[ic][id].c3,m_1,s[ic][id].c4);
				_vector_mulc(r[ic][id].c4,p_1,s[ic][id].c3);
			}
		}
	}
	else
	{
      error_root(1,1,"cc_gamma_row [baryons.c]",
                 "Gamma matrix %d not defined ",mu);
	}
}


void cc_gamma_column(int mu, spinor_dble s[3][4], spinor_dble r[3][4])
{
	complex_dble i,m_i,m_1,p_1;
	int ic;
	
	i.re=0.0;
	i.im=1.0;
	
	m_i.re=0.0;
	m_i.im=-1.0;
	
	m_1.re=-1.0;
	m_1.im=0.0;
	
	p_1.re=1.0;
	p_1.im=0.0;
	
	for (ic=0;ic<3;ic++)
	{
		r[ic][0]=sd0;
		r[ic][1]=sd0;
		r[ic][2]=sd0;
		r[ic][3]=sd0;
	}
	
	if (mu==0)
	{
		for (ic=0;ic<3;ic++)
		{
			_spinor_mulc_add_assign(r[ic][0],m_1,s[ic][3]);
			_spinor_mulc_add_assign(r[ic][1],p_1,s[ic][2]);
			_spinor_mulc_add_assign(r[ic][2],p_1,s[ic][1]);
			_spinor_mulc_add_assign(r[ic][3],m_1,s[ic][0]);
		}
	}
	else if (mu==1)
	{
		for (ic=0;ic<3;ic++)
		{
			_spinor_mulc_add_assign(r[ic][0],m_i,s[ic][2]);
			_spinor_mulc_add_assign(r[ic][1],i,s[ic][3]);
			_spinor_mulc_add_assign(r[ic][2],m_i,s[ic][0]);
			_spinor_mulc_add_assign(r[ic][3],i,s[ic][1]);
		}
	}
	else if (mu==2)
	{
		for (ic=0;ic<3;ic++)
		{
			_spinor_mulc_add_assign(r[ic][0],m_1,s[ic][2]);
			_spinor_mulc_add_assign(r[ic][1],m_1,s[ic][3]);
			_spinor_mulc_add_assign(r[ic][2],m_1,s[ic][0]);
			_spinor_mulc_add_assign(r[ic][3],m_1,s[ic][1]);
		}
	}
	else if (mu==3)
	{
		for (ic=0;ic<3;ic++)
		{
			_spinor_mulc_add_assign(r[ic][0],m_i,s[ic][3]);
			_spinor_mulc_add_assign(r[ic][1],m_i,s[ic][2]);
			_spinor_mulc_add_assign(r[ic][2],m_i,s[ic][1]);
			_spinor_mulc_add_assign(r[ic][3],m_i,s[ic][0]);
		}
	}
	else if (mu==4)
	{
		for (ic=0;ic<3;ic++)
		{
			_spinor_mulc_add_assign(r[ic][0],m_i,s[ic][1]);
			_spinor_mulc_add_assign(r[ic][1],m_i,s[ic][0]);
			_spinor_mulc_add_assign(r[ic][2],i,s[ic][3]);
			_spinor_mulc_add_assign(r[ic][3],i,s[ic][2]);
		}
	}
	else if (mu==5)
	{
		for (ic=0;ic<3;ic++)
		{
			_spinor_mulc_add_assign(r[ic][0],m_1,s[ic][1]);
			_spinor_mulc_add_assign(r[ic][1],p_1,s[ic][0]);
			_spinor_mulc_add_assign(r[ic][2],m_1,s[ic][3]);
			_spinor_mulc_add_assign(r[ic][3],p_1,s[ic][2]);
		}
	}
	else
	{
      error_root(1,1,"cc_gamma_column [baryons.c]",
                 "Gamma matrix %d not defined ",mu);
	}
	
}


extern void twopt_baryons_dble(int gam, int gamp, int num_c, int *mu, int nmom_max, int *momenta,
                        full_spinor_dble *svk, full_spinor_dble *svl, double pm,
                        complex_dble *uds, complex_dble *uud, complex_dble *uuu)
{
  int ic,ic1,id,id1;
  int help1, help2, help3, help4, help5, help6;
  double L[4];
  complex_dble corr0,corr1,corr2,corr3;
  spinor_dble fprop1[3][4],fprop2[3][4],fprop3[3][4];
  su3_vector_dble fcadd;
  su3_vector_dble fcontr1[3][4][3][4];
  su3_vector_dble fcontr2[3];
  spinor_dble *f_help1,*f_help2;

  su3_vector_dble *sqq0=NULL,*sqq1=NULL,*sqq2=NULL;
  int x0,x1,x2,x3,ix;
  int i,ip,mu1,mu2,imu1,imu2;
  double phase;
  complex_dble *st1=NULL,*sx1=NULL,*sy1=NULL,*sz1=NULL,*xx1=NULL;
  complex_dble *st2=NULL,*sx2=NULL,*sy2=NULL,*sz2=NULL,*xx2=NULL;
  complex_dble *st3=NULL,*sx3=NULL,*sy3=NULL,*sz3=NULL,*xx3=NULL;

  int num_mom=nmom_max;

  double PI;
  PI=0.5*6.283185307179586476925286;

  L[0] = (double) NPROC0 * L0;
  L[1] = (double) NPROC1 * L1;
  L[2] = (double) NPROC2 * L2;
  L[3] = (double) NPROC3 * L3;

  help1 = num_c*num_mom;
  help2 = help1*NPROC0*L0;
  xx1=amalloc(help1*sizeof(complex_dble),3);
  xx2=amalloc(help1*sizeof(complex_dble),3);
  xx3=amalloc(help1*sizeof(complex_dble),3);
  sx1=amalloc(help1*sizeof(complex_dble),3);
  sx2=amalloc(help1*sizeof(complex_dble),3);
  sx3=amalloc(help1*sizeof(complex_dble),3);
  sy1=amalloc(help1*sizeof(complex_dble),3);
  sy2=amalloc(help1*sizeof(complex_dble),3);
  sy3=amalloc(help1*sizeof(complex_dble),3);
  sz1=amalloc(help1*sizeof(complex_dble),3);
  sz2=amalloc(help1*sizeof(complex_dble),3);
  sz3=amalloc(help1*sizeof(complex_dble),3);

  st1=amalloc(help2*sizeof(complex_dble),3);
  st2=amalloc(help2*sizeof(complex_dble),3);
  st3=amalloc(help2*sizeof(complex_dble),3);


  for (i=0;i<help2;i++)
  {
   (*(st1+i)).re=(*(uds+i)).re;
   (*(st1+i)).im=(*(uds+i)).im;
   (*(st2+i)).re=(*(uud+i)).re;
   (*(st2+i)).im=(*(uud+i)).im;
   (*(st3+i)).re=(*(uuu+i)).re;
   (*(st3+i)).im=(*(uuu+i)).im;
  }

  for (x0=0;x0<L0;x0++)
  {
   flush_corr2(num_c,num_mom,sx1);
   flush_corr2(num_c,num_mom,sx2);
   flush_corr2(num_c,num_mom,sx3);
   for (x1=0;x1<L1;x1++)
   {
    flush_corr2(num_c,num_mom,sy1);
    flush_corr2(num_c,num_mom,sy2);
    flush_corr2(num_c,num_mom,sy3);
    for (x2=0;x2<L2;x2++)
    {
     flush_corr2(num_c,num_mom,sz1);
     flush_corr2(num_c,num_mom,sz2);
     flush_corr2(num_c,num_mom,sz3);
     for (x3=0;x3<L3;x3++)
     {
      ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];

      i=0;
      for(imu1=0;imu1<num_c;imu1++)
      {
   /*    for(imu2=0;imu2<num_c;imu2++)
       {*/

        imu2=imu1;
        
        
        mu1 = mu[imu1];
        mu2 = mu[imu2];

        f_help1 = (spinor_dble*)&fprop1[0][0];
        f_help2 = (spinor_dble*)&fprop2[0][0];

        for (id=0;id<4;id++){
	         copy_fs_sv(1,(svl+ix),f_help1+id  ,id+1,1);
          copy_fs_sv(1,(svl+ix),f_help1+id+4,id+1,2);
          copy_fs_sv(1,(svl+ix),f_help1+id+8,id+1,3);

          copy_fs_sv(1,(svk+ix),f_help2+id  ,id+1,1);
          copy_fs_sv(1,(svk+ix),f_help2+id+4,id+1,2);
          copy_fs_sv(1,(svk+ix),f_help2+id+8,id+1,3);
        }


        corr1.re=0.0;
        corr1.im=0.0;
        corr0.re=0.0;
        corr0.im=0.0;

       	cc_gamma_row(imu1,fprop1,fprop3);
	       cc_gamma_column(imu2,fprop3,fprop1); 

        for (ic=0;ic<3;ic++){
         for (id=0;id<4;id++){
          sqq1=(su3_vector_dble*)(&fprop1[ic][id]);

/*          #if (defined SSE2)
           _prefetch_spinor(sqq1);
          #endif
 */
          fcontr1[ic][id][0][0] = svd0;
          fcontr1[ic][id][0][1] = svd0;
          fcontr1[ic][id][0][2] = svd0;
          fcontr1[ic][id][0][3] = svd0;

          fcontr1[ic][id][1][0] = svd0;
          fcontr1[ic][id][1][1] = svd0;
          fcontr1[ic][id][1][2] = svd0;
          fcontr1[ic][id][1][3] = svd0;

          fcontr1[ic][id][2][0] = svd0;
          fcontr1[ic][id][2][1] = svd0;
          fcontr1[ic][id][2][2] = svd0;
          fcontr1[ic][id][2][3] = svd0;

         #if (defined SSE3)
           _prefetch_spinor(sqq1);
          #endif

          for (ic1=0;ic1<3;ic1++){

           #if (defined SSE3)
            _prefetch_spinor(&fcontr1[ic][id][ic1][0]);
            _prefetch_spinor((&fprop2[ic1][0]));
            _prefetch_spinor((&fprop2[ic1][1]));
            sqq2=(su3_vector_dble*)(&fprop2[ic1][0]);
            _sse_load_dble(fcontr1[ic][id][ic1][0]);

            _sse_cross_color_add(*(sqq1+0),*(sqq2+0));
            _sse_cross_color_add(*(sqq1+1),*(sqq2+1));
            _sse_cross_color_add(*(sqq1+2),*(sqq2+2));
            _sse_cross_color_add(*(sqq1+3),*(sqq2+3));

            _sse_store_dble(fcontr1[ic][id][ic1][0]);

            _prefetch_spinor((&fprop2[ic1][2]));
            sqq2=(su3_vector_dble*)(&fprop2[ic1][1]);
            _sse_load_dble(fcontr1[ic][id][ic1][1]);

            _sse_cross_color_add(*(sqq1+0),*(sqq2+0));
            _sse_cross_color_add(*(sqq1+1),*(sqq2+1));
            _sse_cross_color_add(*(sqq1+2),*(sqq2+2));
            _sse_cross_color_add(*(sqq1+3),*(sqq2+3));

            _sse_store_dble(fcontr1[ic][id][ic1][1]);

            _prefetch_spinor((&fprop2[ic1][3]));
            sqq2=(su3_vector_dble*)(&fprop2[ic1][2]);
            _sse_load_dble(fcontr1[ic][id][ic1][2]);

            _sse_cross_color_add(*(sqq1+0),*(sqq2+0));
            _sse_cross_color_add(*(sqq1+1),*(sqq2+1));
            _sse_cross_color_add(*(sqq1+2),*(sqq2+2));
            _sse_cross_color_add(*(sqq1+3),*(sqq2+3));

            _sse_store_dble(fcontr1[ic][id][ic1][2]);

            sqq2=(su3_vector_dble*)(&fprop2[ic1][3]);
            _sse_load_dble(fcontr1[ic][id][ic1][3]);

            _sse_cross_color_add(*(sqq1+0),*(sqq2+0));
            _sse_cross_color_add(*(sqq1+1),*(sqq2+1));
            _sse_cross_color_add(*(sqq1+2),*(sqq2+2));
            _sse_cross_color_add(*(sqq1+3),*(sqq2+3));

            _sse_store_dble(fcontr1[ic][id][ic1][3]);
           #else
            for (id1=0;id1<4;id1++){
             sqq2=(su3_vector_dble*)(&fprop2[ic1][id1]);

             _cross_color(fcadd,*(sqq1+0),*(sqq2+0));
             _vector_add_assign(fcontr1[ic][id][ic1][id1],fcadd);
             _cross_color(fcadd,*(sqq1+1),*(sqq2+1));
             _vector_add_assign(fcontr1[ic][id][ic1][id1],fcadd);
             _cross_color(fcadd,*(sqq1+2),*(sqq2+2));
             _vector_add_assign(fcontr1[ic][id][ic1][id1],fcadd);
             _cross_color(fcadd,*(sqq1+3),*(sqq2+3));
	            _vector_add_assign(fcontr1[ic][id][ic1][id1],fcadd);
            } 
           #endif        
          }
         }
        }


        /** This is the first term **/
        fcontr2[0] = svd0;
        fcontr2[1] = svd0;
        fcontr2[2] = svd0;

        for (id=0;id<4;id++)
        {
	         _vector_add_assign(fcontr2[0],fcontr1[1][id][2][id]);
	         _vector_sub_assign(fcontr2[0],fcontr1[2][id][1][id]);
        	 _vector_add_assign(fcontr2[1],fcontr1[2][id][0][id]);
	         _vector_sub_assign(fcontr2[1],fcontr1[0][id][2][id]);
	         _vector_add_assign(fcontr2[2],fcontr1[0][id][1][id]);
	         _vector_sub_assign(fcontr2[2],fcontr1[1][id][0][id]);
        }

        sqq1=(su3_vector_dble*)(&fprop2[0][gamp]);
        corr1.re += _resprod_re(fcontr2[0],*(sqq1+gam));
        corr1.im += _resprod_im(fcontr2[0],*(sqq1+gam));

        sqq1=(su3_vector_dble*)(&fprop2[1][gamp]);
        corr1.re += _resprod_re(fcontr2[1],*(sqq1+gam));
        corr1.im += _resprod_im(fcontr2[1],*(sqq1+gam));

        sqq1=(su3_vector_dble*)(&fprop2[2][gamp]);
        corr1.re += _resprod_re(fcontr2[2],*(sqq1+gam));
       	corr1.im += _resprod_im(fcontr2[2],*(sqq1+gam));


        /**  And now the second term...                **/
        for (id=0;id<4;id++)
        {
         sqq0=(su3_vector_dble*)(&fprop2[0][id]);
         sqq1=(su3_vector_dble*)(&fprop2[1][id]);
         sqq2=(su3_vector_dble*)(&fprop2[2][id]);

         corr0.re += _resprod_re(fcontr1[0][id][1][gamp],*(sqq2+gam));
         corr0.re -= _resprod_re(fcontr1[0][id][2][gamp],*(sqq1+gam));
         corr0.re += _resprod_re(fcontr1[1][id][2][gamp],*(sqq0+gam));
         corr0.re -= _resprod_re(fcontr1[1][id][0][gamp],*(sqq2+gam));
         corr0.re += _resprod_re(fcontr1[2][id][0][gamp],*(sqq1+gam));
         corr0.re -= _resprod_re(fcontr1[2][id][1][gamp],*(sqq0+gam));

         corr0.im += _resprod_im(fcontr1[0][id][1][gamp],*(sqq2+gam));
         corr0.im -= _resprod_im(fcontr1[0][id][2][gamp],*(sqq1+gam));
         corr0.im += _resprod_im(fcontr1[1][id][2][gamp],*(sqq0+gam));
         corr0.im -= _resprod_im(fcontr1[1][id][0][gamp],*(sqq2+gam));
         corr0.im += _resprod_im(fcontr1[2][id][0][gamp],*(sqq1+gam));
         corr0.im -= _resprod_im(fcontr1[2][id][1][gamp],*(sqq0+gam));
        }

        corr0.re = corr0.re*pm;
        corr0.im = corr0.im*pm;

        /**  the 3 channels **/

        corr1.re = corr1.re*pm;
        corr1.im = corr1.im*pm;

        corr2.re = corr1.re + corr0.re;  
        corr2.im = corr1.im + corr0.im;

        corr3.re = 2.0*corr1.re + 4.0*corr0.re;  
        corr3.im = 2.0*corr1.im + 4.0*corr0.im;


        /** final part **/
        help1 = cpr[1]*L1+x1;
        help2 = cpr[2]*L2+x2;
        help3 = cpr[3]*L3+x3;
        for (ip=0;ip<num_mom;ip++)
        {

         phase = 2.0*PI*(
	          (double)((*(momenta+ip*3  ))*help1)/L[1] +
           (double)((*(momenta+ip*3+1))*help2)/L[2] +
           (double)((*(momenta+ip*3+2))*help3)/L[2]);

         (*(xx1+i)).re = corr1.re*cos(phase) - corr1.im*sin(phase);
         (*(xx1+i)).im = corr1.im*cos(phase) + corr1.re*sin(phase);
         (*(xx2+i)).re = corr2.re*cos(phase) - corr2.im*sin(phase);
         (*(xx2+i)).im = corr2.im*cos(phase) + corr2.re*sin(phase);
         (*(xx3+i)).re = corr3.re*cos(phase) - corr3.im*sin(phase);
         (*(xx3+i)).im = corr3.im*cos(phase) + corr3.re*sin(phase);
         i=i+1;

        }
      /* }*/
      }

      add_corr2(num_c,num_mom,xx1,sz1);
      add_corr2(num_c,num_mom,xx2,sz2);
      add_corr2(num_c,num_mom,xx3,sz3);
     }

     add_corr2(num_c,num_mom,sz1,sy1);
     add_corr2(num_c,num_mom,sz2,sy2);
     add_corr2(num_c,num_mom,sz3,sy3);
    }

    add_corr2(num_c,num_mom,sy1,sx1);
    add_corr2(num_c,num_mom,sy2,sx2);
    add_corr2(num_c,num_mom,sy3,sx3);
   }

   help4 = cpr[0]*L0 + x0;
   help5 = L0*NPROC0;
   for (i=0;i<num_c*num_mom;i++)
   {
    (*(st1+help4 + i*help5)).re+=(*(sx1+i)).re;
    (*(st1+help4 + i*help5)).im+=(*(sx1+i)).im;
    (*(st2+help4 + i*help5)).re+=(*(sx2+i)).re;
    (*(st2+help4 + i*help5)).im+=(*(sx2+i)).im;
    (*(st3+help4 + i*help5)).re+=(*(sx3+i)).re;
    (*(st3+help4 + i*help5)).im+=(*(sx3+i)).im;
   }

  }

  for (i=0;i<NPROC0*L0*num_c*num_mom;i++)
  {
   (*(uds+i)).re=(*(st1+i)).re;
   (*(uds+i)).im=(*(st1+i)).im;
   (*(uud+i)).re=(*(st2+i)).re;
   (*(uud+i)).im=(*(st2+i)).im;
   (*(uuu+i)).re=(*(st3+i)).re;
   (*(uuu+i)).im=(*(st3+i)).im;
  }

  afree(xx1);
  afree(xx2);
  afree(xx3);
  afree(sx1);
  afree(sx2);
  afree(sx3);
  afree(sy1);
  afree(sy2);
  afree(sy3);
  afree(sz1);
  afree(sz2);
  afree(sz3);
  afree(st1);
  afree(st2);
  afree(st3);

}


extern void twopt_omega_dble(int gam, int gamp, int num_c, int *mu, int num_mom,
                             int *momenta, full_spinor_dble *svk, double pm,
                             complex_dble *sst)
{
  int ic,ic1,id,id1,id2;
  complex_dble corr;
  spinor_dble fprop1[3][4],fprop2[3][4],fprop3[3][4],fprop4[3][4],fprop5[3][4];
  su3_vector_dble fcadd;
  su3_vector_dble fcontr1[3][4][3][4];
  su3_vector_dble fcontr2[3];
  spinor_dble *f_help;

  su3_vector_dble *sqq0=NULL,*sqq1=NULL,*sqq2=NULL;
  int x0,x1,x2,x3,ix;
  int i,ip,mu1,mu2,imu1,imu2;
  double phase;
  complex_dble *st=NULL,*sx=NULL,*sy=NULL,*sz=NULL,*x=NULL;

  double PI;
  PI=0.5*6.283185307179586476925286;


  x =amalloc(num_c*num_c*num_mom*sizeof(complex_dble),3);
  sx=amalloc(num_c*num_c*num_mom*sizeof(complex_dble),3);
  sy=amalloc(num_c*num_c*num_mom*sizeof(complex_dble),3);
  sz=amalloc(num_c*num_c*num_mom*sizeof(complex_dble),3);
  st=amalloc(num_c*num_c*num_mom*NPROC0*L0*sizeof(complex_dble),3);


  for (i=0;i<NPROC0*L0*num_c*num_c*num_mom;i++)
  {
   (*(st+i)).re=(*(sst+i)).re;
   (*(st+i)).im=(*(sst+i)).im;
  }

  for (x0=0;x0<L0;x0++)
  {
   flush_corr(num_c,num_mom,sx);
   for (x1=0;x1<L1;x1++)
   {
    flush_corr(num_c,num_mom,sy);
    for (x2=0;x2<L2;x2++)
    {
     flush_corr(num_c,num_mom,sz);
     for (x3=0;x3<L3;x3++)
     {
      ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];

      i=0;

      for(imu1=0;imu1<num_c;imu1++)
      {
       for(imu2=0;imu2<num_c;imu2++)
       {

        mu1 = mu[imu1];
        mu2 = mu[imu2];

        f_help = (spinor_dble*)&fprop1[0][0];

        for (id=0;id<4;id++)
        {
	 for (ic=0;ic<3;ic++)
	 {
	  copy_fs_sv(1,(svk+ix),f_help+id+4*ic,id+1,ic+1);
          fprop2[ic][id] = fprop1[ic][id];
          fprop3[ic][id] = fprop1[ic][id];
         }
        }

        corr.re=0.0;
        corr.im=0.0;


        cc_gamma_row(imu1,fprop1,fprop4);
	cc_gamma_column(imu2,fprop2,fprop5); 
        for (id=0;id<4;id++)
        {
         for (ic=0;ic<3;ic++)
         {
          fprop1[ic][id] = fprop4[ic][id];
          fprop2[ic][id] = fprop5[ic][id];
         }
        }

/** ... preparation for the third and sixth terms ... **/

        for (id=0;id<4;id++)
        {
         for (ic=0;ic<3;ic++)
         {
          for (id1=0;id1<4;id1++)
          {
	   for (ic1=0;ic1<3;ic1++)
           {
	      
            fcontr1[ic][id][ic1][id1] = svd0;

            sqq1=(su3_vector_dble*)(&fprop1[ic][id]);
            sqq2=(su3_vector_dble*)(&fprop2[ic1][id1]);

            for (id2=0;id2<4;id2++)
            {                            
             _cross_color(fcadd,*(sqq1+id2),*(sqq2+id2));
	     _vector_add_assign(fcontr1[ic][id][ic1][id1],fcadd);
            }  

  	   }
          }
         }
        }

/** This is the third term **/

        fcontr2[0] = svd0;
        fcontr2[1] = svd0;
        fcontr2[2] = svd0;

        for (id=0;id<4;id++)
        {
	 _vector_add_assign(fcontr2[0],fcontr1[1][id][2][id]);
	 _vector_sub_assign(fcontr2[0],fcontr1[2][id][1][id]);
	 _vector_add_assign(fcontr2[1],fcontr1[2][id][0][id]);
	 _vector_sub_assign(fcontr2[1],fcontr1[0][id][2][id]);
	 _vector_add_assign(fcontr2[2],fcontr1[0][id][1][id]);
	 _vector_sub_assign(fcontr2[2],fcontr1[1][id][0][id]);
        }

        for (ic=0;ic<3;ic++)
        {
         sqq1=(su3_vector_dble*)(&fprop3[ic][gamp]);

         corr.re += _resprod_re(fcontr2[ic],*(sqq1+gam));
	 corr.im += _resprod_im(fcontr2[ic],*(sqq1+gam));
        }

/** This is the sixth term **/

        for (id=0;id<4;id++)
        {
         sqq0=(su3_vector_dble*)(&fprop3[0][id]);
         sqq1=(su3_vector_dble*)(&fprop3[1][id]);
         sqq2=(su3_vector_dble*)(&fprop3[2][id]);

         corr.re += _resprod_re(fcontr1[0][gamp][1][id],*(sqq2+gam));
         corr.re -= _resprod_re(fcontr1[0][gamp][2][id],*(sqq1+gam));
         corr.re += _resprod_re(fcontr1[1][gamp][2][id],*(sqq0+gam));
         corr.re -= _resprod_re(fcontr1[1][gamp][0][id],*(sqq2+gam));
         corr.re += _resprod_re(fcontr1[2][gamp][0][id],*(sqq1+gam));
         corr.re -= _resprod_re(fcontr1[2][gamp][1][id],*(sqq0+gam));

         corr.im += _resprod_im(fcontr1[0][gamp][1][id],*(sqq2+gam));
         corr.im -= _resprod_im(fcontr1[0][gamp][2][id],*(sqq1+gam));
         corr.im += _resprod_im(fcontr1[1][gamp][2][id],*(sqq0+gam));
         corr.im -= _resprod_im(fcontr1[1][gamp][0][id],*(sqq2+gam));
         corr.im += _resprod_im(fcontr1[2][gamp][0][id],*(sqq1+gam));
         corr.im -= _resprod_im(fcontr1[2][gamp][1][id],*(sqq0+gam));
        }

/** ... preparation for the fourth and fifth terms ... **/

        for (id=0;id<4;id++)
        {
         for (ic=0;ic<3;ic++)
         {
          for (id1=0;id1<4;id1++)
          {
	   for (ic1=0;ic1<3;ic1++)
           {
	      
	    fcontr1[ic][id][ic1][id1] = svd0;

            sqq1=(su3_vector_dble*)(&fprop1[ic][id]);
            sqq2=(su3_vector_dble*)(&fprop3[ic1][id1]);

            for (id2=0;id2<4;id2++)
            {                            
             _cross_color(fcadd,*(sqq1+id2),*(sqq2+id2));
	     _vector_add_assign(fcontr1[ic][id][ic1][id1],fcadd);
            }  

  	   }
          }
         }
        }

/** This is the fourth term **/

        for (id=0;id<4;id++)
        {
         sqq0=(su3_vector_dble*)(&fprop2[0][id]);
         sqq1=(su3_vector_dble*)(&fprop2[1][id]);
         sqq2=(su3_vector_dble*)(&fprop2[2][id]);

         corr.re += _resprod_re(fcontr1[0][id][1][gamp],*(sqq2+gam));
         corr.re -= _resprod_re(fcontr1[0][id][2][gamp],*(sqq1+gam));
         corr.re += _resprod_re(fcontr1[1][id][2][gamp],*(sqq0+gam));
         corr.re -= _resprod_re(fcontr1[1][id][0][gamp],*(sqq2+gam));
         corr.re += _resprod_re(fcontr1[2][id][0][gamp],*(sqq1+gam));
         corr.re -= _resprod_re(fcontr1[2][id][1][gamp],*(sqq0+gam));

         corr.im += _resprod_im(fcontr1[0][id][1][gamp],*(sqq2+gam));
         corr.im -= _resprod_im(fcontr1[0][id][2][gamp],*(sqq1+gam));
         corr.im += _resprod_im(fcontr1[1][id][2][gamp],*(sqq0+gam));
         corr.im -= _resprod_im(fcontr1[1][id][0][gamp],*(sqq2+gam));
         corr.im += _resprod_im(fcontr1[2][id][0][gamp],*(sqq1+gam));
         corr.im -= _resprod_im(fcontr1[2][id][1][gamp],*(sqq0+gam));
        }

/** This is the fifth term **/

        for (id=0;id<4;id++)
        {
         sqq0=(su3_vector_dble*)(&fprop2[0][id]);
         sqq1=(su3_vector_dble*)(&fprop2[1][id]);
         sqq2=(su3_vector_dble*)(&fprop2[2][id]);

         corr.re += _resprod_re(fcontr1[0][gamp][1][id],*(sqq2+gam));
         corr.re -= _resprod_re(fcontr1[0][gamp][2][id],*(sqq1+gam));
         corr.re += _resprod_re(fcontr1[1][gamp][2][id],*(sqq0+gam));
         corr.re -= _resprod_re(fcontr1[1][gamp][0][id],*(sqq2+gam));
         corr.re += _resprod_re(fcontr1[2][gamp][0][id],*(sqq1+gam));
         corr.re -= _resprod_re(fcontr1[2][gamp][1][id],*(sqq0+gam));

         corr.im += _resprod_im(fcontr1[0][gamp][1][id],*(sqq2+gam));
         corr.im -= _resprod_im(fcontr1[0][gamp][2][id],*(sqq1+gam));
         corr.im += _resprod_im(fcontr1[1][gamp][2][id],*(sqq0+gam));
         corr.im -= _resprod_im(fcontr1[1][gamp][0][id],*(sqq2+gam));
         corr.im += _resprod_im(fcontr1[2][gamp][0][id],*(sqq1+gam));
         corr.im -= _resprod_im(fcontr1[2][gamp][1][id],*(sqq0+gam));
        }

/** ... preparation for the first and second terms ... **/

	cc_gamma_column(imu2,fprop1,fprop4); 
        for (id=0;id<4;id++)
        {
         for (ic=0;ic<3;ic++)
         {
          fprop1[ic][id] = fprop4[ic][id];
         }
        }


        for (id=0;id<4;id++)
        {
         for (ic=0;ic<3;ic++)
         {
          for (id1=0;id1<4;id1++)
          {
	   for (ic1=0;ic1<3;ic1++)
           {
	      
	    fcontr1[ic][id][ic1][id1] = svd0;

            sqq1=(su3_vector_dble*)(&fprop1[ic][id]);
            sqq2=(su3_vector_dble*)(&fprop3[ic1][id1]);

            for (id2=0;id2<4;id2++)
            {                            
             _cross_color(fcadd,*(sqq1+id2),*(sqq2+id2));
	     _vector_add_assign(fcontr1[ic][id][ic1][id1],fcadd);
            }  

  	   }
          }
         }
        }
        

/** This is the first term **/

        fcontr2[0] = svd0;
        fcontr2[1] = svd0;
        fcontr2[2] = svd0;

        for (id=0;id<4;id++)
        {
	 _vector_add_assign(fcontr2[0],fcontr1[1][id][2][id]);
	 _vector_sub_assign(fcontr2[0],fcontr1[2][id][1][id]);
	 _vector_add_assign(fcontr2[1],fcontr1[2][id][0][id]);
	 _vector_sub_assign(fcontr2[1],fcontr1[0][id][2][id]);
	 _vector_add_assign(fcontr2[2],fcontr1[0][id][1][id]);
	 _vector_sub_assign(fcontr2[2],fcontr1[1][id][0][id]);
        }

        for (ic=0;ic<3;ic++)
        {
         sqq1=(su3_vector_dble*)(&fprop3[ic][gamp]);

         corr.re += _resprod_re(fcontr2[ic],*(sqq1+gam));
	 corr.im += _resprod_im(fcontr2[ic],*(sqq1+gam));
        }


/** This is the second term **/

        for (id=0;id<4;id++)
        {
         sqq0=(su3_vector_dble*)(&fprop3[0][id]);
         sqq1=(su3_vector_dble*)(&fprop3[1][id]);
         sqq2=(su3_vector_dble*)(&fprop3[2][id]);

         corr.re += _resprod_re(fcontr1[0][id][1][gamp],*(sqq2+gam));
         corr.re -= _resprod_re(fcontr1[0][id][2][gamp],*(sqq1+gam));
         corr.re += _resprod_re(fcontr1[1][id][2][gamp],*(sqq0+gam));
         corr.re -= _resprod_re(fcontr1[1][id][0][gamp],*(sqq2+gam));
         corr.re += _resprod_re(fcontr1[2][id][0][gamp],*(sqq1+gam));
         corr.re -= _resprod_re(fcontr1[2][id][1][gamp],*(sqq0+gam));

         corr.im += _resprod_im(fcontr1[0][id][1][gamp],*(sqq2+gam));
         corr.im -= _resprod_im(fcontr1[0][id][2][gamp],*(sqq1+gam));
         corr.im += _resprod_im(fcontr1[1][id][2][gamp],*(sqq0+gam));
         corr.im -= _resprod_im(fcontr1[1][id][0][gamp],*(sqq2+gam));
         corr.im += _resprod_im(fcontr1[2][id][0][gamp],*(sqq1+gam));
         corr.im -= _resprod_im(fcontr1[2][id][1][gamp],*(sqq0+gam));
        }


        corr.re = corr.re*pm;
        corr.im = corr.im*pm;

        for (ip=0;ip<num_mom;ip++)
        {

         phase = 2.0*PI*(
	   (double)((*(momenta+ip*3  ))*(cpr[1]*L1+x1))/(double)(NPROC1*L1) +
           (double)((*(momenta+ip*3+1))*(cpr[2]*L2+x2))/(double)(NPROC2*L2) +
           (double)((*(momenta+ip*3+2))*(cpr[3]*L3+x3))/(double)(NPROC3*L3));

         (*(x+i)).re = corr.re*cos(phase) - corr.im*sin(phase);
         (*(x+i)).im = corr.im*cos(phase) + corr.re*sin(phase);
         i=i+1;

        }
       }
      }

      add_corr(num_c,num_mom,x,sz);
     }

     add_corr(num_c,num_mom,sz,sy);
    }

    add_corr(num_c,num_mom,sy,sx);
   }

   for (i=0;i<num_c*num_c*num_mom;i++)
   {
    (*(st+cpr[0]*L0 + x0 + i*L0*NPROC0)).re+=(*(sx+i)).re;
    (*(st+cpr[0]*L0 + x0 + i*L0*NPROC0)).im+=(*(sx+i)).im;
   }

  }

  for (i=0;i<NPROC0*L0*num_c*num_c*num_mom;i++)
  {
   (*(sst+i)).re=(*(st+i)).re;
   (*(sst+i)).im=(*(st+i)).im;
  }

  afree(x);
  afree(sx);
  afree(sy);
  afree(sz);
  afree(st);

}


extern void d_source_dble(int t, int gam, int gamp, full_spinor_dble *svk, double pm, 
								  spinor_dble *grand_source)
{
  int ic,ic1,id,id1;
  spinor_dble fprop1[3][4],fprop2[3][4],fprop3[3][4],fprop4[3][4],fprop5[3][4];
  spinor_dble fcontr1[3][3][4];
  spinor_dble fcontr2[3][4];
  spinor_dble *f_help;

  su3_vector_dble *sqq0=NULL,*sqq1=NULL,*sqq2=NULL;
  complex_dble *cqq0=NULL,*cqq1=NULL;

  int ix;

  for (ix=0;ix<VOLUME;ix++)
  {
   if (coords[ix].t==t)
   {

    f_help = (spinor_dble*)&fprop1[0][0];

    for (id=0;id<4;id++)
    {
     for (ic=0;ic<3;ic++)
     {
      copy_fs_sv(1,(svk+ix),f_help+id+4*ic,id+1,ic+1);
      fprop2[ic][id] = fprop1[ic][id];
      fprop3[ic][id] = fprop1[ic][id];
      fcontr2[ic][id] = sd0;
     }
    }


/** This is the first term **/

    cc_gamma_row(5,fprop1,fprop4);
    cc_gamma_column(5,fprop4,fprop1); 

    for (id=0;id<4;id++)
    {
     for (ic=0;ic<3;ic++)
     {
      for (ic1=0;ic1<3;ic1++)
      {
	      
       sqq0=(su3_vector_dble*)(&fcontr1[ic][ic1][id]);
       sqq1=(su3_vector_dble*)(&fprop1[ic][id]);
       sqq2=(su3_vector_dble*)(&fprop2[ic1][gamp]);

       for (id1=0;id1<4;id1++)
       {
        _cross_color(*(sqq0+id1),*(sqq1+id1),*(sqq2+gam));
       }

      }
     }
    }

    for (id=0;id<4;id++)
    {
     _spinor_sub_assign(fcontr2[0][id],fcontr1[1][2][id]);
     _spinor_add_assign(fcontr2[0][id],fcontr1[2][1][id]);
     _spinor_sub_assign(fcontr2[1][id],fcontr1[2][0][id]);
     _spinor_add_assign(fcontr2[1][id],fcontr1[0][2][id]);
     _spinor_sub_assign(fcontr2[2][id],fcontr1[0][1][id]);
     _spinor_add_assign(fcontr2[2][id],fcontr1[1][0][id]);
    }      
      

/** This is the second term **/

    cc_gamma_row(5,fprop3,fprop4);
    cc_gamma_column(5,fprop2,fprop5);  
    for (id=0;id<4;id++)
    {
     for (ic=0;ic<3;ic++)
     {
      fprop3[ic][id] = fprop4[ic][id];
      fprop2[ic][id] = fprop5[ic][id];
     }
    }

    for (id=0;id<4;id++)
    {
     for (ic=0;ic<3;ic++)
     {
      for (ic1=0;ic1<3;ic1++)
      {
	      
       sqq0=(su3_vector_dble*)(&fcontr1[ic][ic1][id]);
       sqq1=(su3_vector_dble*)(&fprop3[ic][gamp]);
       sqq2=(su3_vector_dble*)(&fprop2[ic1][id]);

       for (id1=0;id1<4;id1++)
       {
        _cross_color(*(sqq0+id1),*(sqq1+id1),*(sqq2+gam));
       }

      }
     }
    }

    for (id=0;id<4;id++)
    {
     _spinor_sub_assign(fcontr2[0][id],fcontr1[1][2][id]);
     _spinor_add_assign(fcontr2[0][id],fcontr1[2][1][id]);
     _spinor_sub_assign(fcontr2[1][id],fcontr1[2][0][id]);
     _spinor_add_assign(fcontr2[1][id],fcontr1[0][2][id]);
     _spinor_sub_assign(fcontr2[2][id],fcontr1[0][1][id]);
     _spinor_add_assign(fcontr2[2][id],fcontr1[1][0][id]);
    }      


/**   this is gamma_5 * adjoint  (in fact, only imaginary conjugation)   **/

    for (id=0;id<4;id++)
    {
     for (ic=0;ic<3;ic++)
     {

      cqq0=(complex_dble*)(grand_source+ix+VOLUME*(ic+3*id));
      cqq1=(complex_dble*)(&fcontr2[ic][id]); 

      for (id1=0;id1<4;id1++)
      {
       for (ic1=0;ic1<3;ic1++)
       {

        if (id1>1)
	{         
         (*(cqq0+ic1+3*id1)).re -= (*(cqq1+ic1+3*id1)).re*pm;
         (*(cqq0+ic1+3*id1)).im += (*(cqq1+ic1+3*id1)).im*pm;
        }
        else
        {
         (*(cqq0+ic1+3*id1)).re += (*(cqq1+ic1+3*id1)).re*pm;
         (*(cqq0+ic1+3*id1)).im -= (*(cqq1+ic1+3*id1)).im*pm;
        }

       }
      }
     }
    }


   }
  }

}


extern void u_source_dble(int t, int gam, int gamp, full_spinor_dble *svk, full_spinor_dble *svl, 
					      	  double pm, spinor_dble *grand_source)
{
  int ic,ic1,id,id1;
  spinor_dble fprop1[3][4],fprop2[3][4],fprop3[3][4];
  su3_vector_dble fcadd;
  spinor_dble fcontr1[3][3][4];
  spinor_dble fcontr2[3][4];
  spinor_dble *f_help1,*f_help2;

  su3_vector_dble *sqq0=NULL,*sqq1=NULL,*sqq2=NULL;
  complex_dble *cqq0=NULL,*cqq1=NULL;

  int ix; 

  for (ix=0;ix<VOLUME;ix++)
  {
   if (coords[ix].t==t)
   {

    f_help1 = (spinor_dble*)&fprop1[0][0];
    f_help2 = (spinor_dble*)&fprop2[0][0];

    for (id=0;id<4;id++)
    {
     for (ic=0;ic<3;ic++)
     {
      copy_fs_sv(1,(svl+ix),f_help1+id+4*ic,id+1,ic+1);
      copy_fs_sv(1,(svk+ix),f_help2+id+4*ic,id+1,ic+1);
      fcontr2[ic][id] = sd0;
     }
    }

    cc_gamma_row(5,fprop1,fprop3);
    cc_gamma_column(5,fprop3,fprop1);


/** This is the first term **/

    for (id=0;id<4;id++)
    {
     for (ic=0;ic<3;ic++)
     {
      for (ic1=0;ic1<3;ic1++)
      {
	      
       sqq0=(su3_vector_dble*)(&fcontr1[ic][ic1][id]);
       sqq1=(su3_vector_dble*)(&fprop1[ic][id]);
       sqq2=(su3_vector_dble*)(&fprop2[ic1][gamp]);

       for (id1=0;id1<4;id1++)
       {
        _cross_color(*(sqq0+id1),*(sqq1+id1),*(sqq2+gam));
       }

      }
     }
    }

    for (id=0;id<4;id++)
    {
     _spinor_sub_assign(fcontr2[0][id],fcontr1[1][2][id]);
     _spinor_add_assign(fcontr2[0][id],fcontr1[2][1][id]);
     _spinor_sub_assign(fcontr2[1][id],fcontr1[2][0][id]);
     _spinor_add_assign(fcontr2[1][id],fcontr1[0][2][id]);
     _spinor_sub_assign(fcontr2[2][id],fcontr1[0][1][id]);
     _spinor_add_assign(fcontr2[2][id],fcontr1[1][0][id]);
    }      


/** This is the second term **/

    for (id=0;id<4;id++)
    {
     for (ic=0;ic<3;ic++)
     {
      for (ic1=0;ic1<3;ic1++)
      {

       fcontr1[ic][ic1][id] = sd0; 

       sqq0=(su3_vector_dble*)(&fcontr1[ic][ic1][id]);
       sqq1=(su3_vector_dble*)(&fprop1[ic][id]);
       sqq2=(su3_vector_dble*)(&fprop2[ic1][gamp]);

       for (id1=0;id1<4;id1++)
       {                            
        _cross_color(fcadd,*(sqq1+id1),*(sqq2+id1));
        _vector_add_assign(*(sqq0+gam),fcadd);
       }  

      }
     }
    }

    for (id=0;id<4;id++)
    {
     _spinor_sub_assign(fcontr2[0][id],fcontr1[1][2][id]);
     _spinor_add_assign(fcontr2[0][id],fcontr1[2][1][id]);
     _spinor_sub_assign(fcontr2[1][id],fcontr1[2][0][id]);
     _spinor_add_assign(fcontr2[1][id],fcontr1[0][2][id]);
     _spinor_sub_assign(fcontr2[2][id],fcontr1[0][1][id]);
     _spinor_add_assign(fcontr2[2][id],fcontr1[1][0][id]);
    }      


/** This is the third term **/

    for (ic=0;ic<3;ic++)
    {
     for (ic1=0;ic1<3;ic1++)
     {
       
      fcontr1[ic][ic1][gamp] = sd0; 

      for (id=0;id<4;id++)
      {
	     
       sqq0=(su3_vector_dble*)(&fcontr1[ic][ic1][gamp]);
       sqq1=(su3_vector_dble*)(&fprop1[ic][id]);
       sqq2=(su3_vector_dble*)(&fprop2[ic1][id]);

       for (id1=0;id1<4;id1++)
       {
        _cross_color(fcadd,*(sqq1+id1),*(sqq2+id1));
        _vector_add_assign(*(sqq0+gam),fcadd);
       }

      }
     }
    }

    _spinor_sub_assign(fcontr2[0][gamp],fcontr1[1][2][gamp]);
    _spinor_add_assign(fcontr2[0][gamp],fcontr1[2][1][gamp]);
    _spinor_sub_assign(fcontr2[1][gamp],fcontr1[2][0][gamp]);
    _spinor_add_assign(fcontr2[1][gamp],fcontr1[0][2][gamp]);
    _spinor_sub_assign(fcontr2[2][gamp],fcontr1[0][1][gamp]);
    _spinor_add_assign(fcontr2[2][gamp],fcontr1[1][0][gamp]);


/** This is the fourth term **/

    for (ic=0;ic<3;ic++)
    {
     for (ic1=0;ic1<3;ic1++)
     {

      fcontr1[ic][ic1][gamp] = sd0; 

      for (id=0;id<4;id++)
      {

       sqq0=(su3_vector_dble*)(&fcontr1[ic][ic1][gamp]);
       sqq1=(su3_vector_dble*)(&fprop1[ic][id]);
       sqq2=(su3_vector_dble*)(&fprop2[ic1][id]);

       for (id1=0;id1<4;id1++)
       {
        _cross_color(fcadd,*(sqq1+id1),*(sqq2+gam));
        _vector_add_assign(*(sqq0+id1),fcadd);
       }
      }

     }
    }
     
    _spinor_sub_assign(fcontr2[0][gamp],fcontr1[1][2][gamp]);
    _spinor_add_assign(fcontr2[0][gamp],fcontr1[2][1][gamp]);
    _spinor_sub_assign(fcontr2[1][gamp],fcontr1[2][0][gamp]);
    _spinor_add_assign(fcontr2[1][gamp],fcontr1[0][2][gamp]);
    _spinor_sub_assign(fcontr2[2][gamp],fcontr1[0][1][gamp]);
    _spinor_add_assign(fcontr2[2][gamp],fcontr1[1][0][gamp]);


/**   this is gamma_5 * adjoint  (in fact, only imaginary conjugation)   **/

    for (id=0;id<4;id++)
    {
     for (ic=0;ic<3;ic++)
     {

      cqq0=(complex_dble*)(grand_source+ix+VOLUME*(ic+3*id));
      cqq1=(complex_dble*)(&fcontr2[ic][id]); 

      for (id1=0;id1<4;id1++)
      {
       for (ic1=0;ic1<3;ic1++)
       {

        if (id1>1)
	       {         
         (*(cqq0+ic1+3*id1)).re -= (*(cqq1+ic1+3*id1)).re*pm;
         (*(cqq0+ic1+3*id1)).im += (*(cqq1+ic1+3*id1)).im*pm;
        }
        else
        {
         (*(cqq0+ic1+3*id1)).re += (*(cqq1+ic1+3*id1)).re*pm;
         (*(cqq0+ic1+3*id1)).im -= (*(cqq1+ic1+3*id1)).im*pm;
        }

       }
      }
     }
    }


   }
  }

}


extern void check_closure_dble(spinor_dble *grand_source,
                        full_spinor_dble *propagator, 
                        full_spinor_dble *ext_prop, 
                        propinfo prop){

  int id,ic,idirac,icolour;
  int my_rank, time[4], i;
  int x0,x1,x2,x3,ix;

  complex_dble resu_test1, resu_test2, resu_h, result;

  full_spinor_dble shelp1,shelp2,*proph;
  spinor_dble **wsd;

  propinfo prophelp = prop;

  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

  wsd=reserve_wsd(1);

  resu_test1.re = 0.0;
  resu_test1.im = 0.0;
  resu_test2.re = 0.0;
  resu_test2.im = 0.0;

  prophelp.src.pos[0] = 0;
  prophelp.src.pos[1] = 0;
  prophelp.src.pos[2] = 0;
  prophelp.src.pos[3] = 0;

  proph = amalloc(VOLUME*sizeof(full_spinor_dble),3);

  if(  (strcmp(prop.src.type, "J_thin") == 0) ||
       (strcmp(prop.src.type, "J_APE")  == 0) ||
       (strcmp(prop.src.type, "J_HYP")  == 0) ){

    for(idirac=1; idirac<=4; idirac++){
      for(icolour=1; icolour<=3; icolour++){
        copy_fs_sv(VOLUME, ext_prop, wsd[0], idirac, icolour);
        jacobi_source(wsd[0], prophelp.src, idirac, icolour);
        copy_sv_fs(VOLUME, wsd[0], proph, idirac, icolour);
      }
    }

  }

  for (x0=0;x0<L0;x0++){
    for (x1=0;x1<L1;x1++){
      for (x2=0;x2<L2;x2++){
        for (x3=0;x3<L3;x3++){
 
          ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];
  
          /* test 1 */
          if( (coords[ix].t==0) &&  
              (coords[ix].x==0) && 
              (coords[ix].y==0) && 
              (coords[ix].z==0) ){

                fsv_trace(proph+ix , &resu_test1);
          }
          /* test 2 */
          if( (coords[ix].t==prop.t_snk) ){
            
            resu_h.re = 0.0;
            resu_h.im = 0.0;

            for(id=0; id<4; id++){
              for(ic=0; ic<3; ic++){
                copy_sv_fs(1,grand_source+ix+VOLUME*(ic+3*id),&shelp1,id+1,ic+1);
              }
            }
            
            mul_gamma_l(5,&shelp1,&shelp2);
            adj_full_spinor(&shelp2,&shelp1);

            meson_trace(propagator+ix,&shelp1,&resu_h);
            add_corr(1,1,&resu_h,&resu_test2);

          }  

        }
      }
    } 
  }

  /* writing the stuff out */
  MPI_Reduce(&resu_test2,&result,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if (my_rank==0){
    printf("------------------------------------------------------------------------\n");
    printf("Closure of extended src. : Re: %.16e, Im: %.16e\n",result.re,-result.im);
  }

  MPI_Reduce(&resu_test1,&result,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if (my_rank==0){
    printf("Closure of extended prop.: Re: %.16e, Im: %.16e\n",result.re,result.im);
    printf("------------------------------------------------------------------------\n");
  }

  release_wsd();
  afree(proph);
}


