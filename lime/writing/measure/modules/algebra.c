/*******************************************************************************
*
* File algebtra.h
*
* Copyright (C) 2008-11 Andreas Juettner
*                       Bastian Knippschild
*
* The externally accessible functions are
*
* extern void copy_sv_fs(int vol,spinor_dble *s,full_spinor_dble *sv,
*  				int id,int ic)
* 	copies spinor_dble *s into full_spinor_dble *sv with
*	dirac column index id and color column index ic
*
* extern void copy_fs_sv(int vol,full_spinor_dble *sv,spinor_dble *s,
*				int id,int ic)
* 	copies full_spinor_dble *sv with
*	dirac column index id and color column index ic
*	into spinor_dble *s with
*
* extern void set_fsd2zero(int vol,full_spinor_dble *pk)
* 	sets entried of full_spinor_dble to zero
*
* extern full_spinor_dble adj_full_spinor(full_spinor_dble s)
*	returns adj of full spinor
*
* extern full_spinor_dble full_mat_prod(full_spinor_dble k,full_spinor_dble l)
*	returns matrix product of two full_spinor_dble
*
* extern complex_dble meson_trace(full_spinor_dble k,full_spinor_dble l)
*	returns the trace over the product of two full_spinor_dble
*
* extern complex_dble fsv_trace(full_spinor_dble r)
*	returns the the trace of a full_spinor_dble
*
* extern su3_dble mul_cplx_su3_dble(complex_dble z,su3_dble s)
*	returns z*s where z is a complex_dble and s a su3_dble
*
* extern full_spinor_dble mul_cmplx_fsv_dble(full_spinor_dble s,
*				complex_dble phase)
*	returns phase*s where s is a full_spinor_dble and phase a complex_dble
*
* extern full_spinor_dble mul_gamma_l(int mu,full_spinor_dble s)
*	left-multiplies a full_spinor_dble with gamma_mu
*
* extern complex_dble spinor_prod_dble_local(spinor_dble *k,spinor_dble *l)
*	returns the result for the local spinor product
*
* extern spinor_dble mul_s_gamma_l(int mu,spinor_dble s)
*	left-multiplies a spinor_dble with gamma_mu
* extern spinor_dble mul_s_gamma_r(int mu,spinor_dble s)
*	right-multiplies a spinor_dble with gamma_mu
*
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/
#define ALGEBRA_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "start.h"
#include "global.h"
#include "measure.h"
#include "mpi.h"

static const su3_dble ud0={{0.0}};
#define _su3_tensor_prod(r,u,s) \
  (r).c11.re=(u).c1.re*(s).c1.re-(u).c1.im*(s).c1.im; \
  (r).c11.im=(u).c1.re*(s).c1.im+(u).c1.im*(s).c1.re; \
  (r).c12.re=(u).c1.re*(s).c2.re-(u).c1.im*(s).c2.im; \
  (r).c12.im=(u).c1.re*(s).c2.im+(u).c1.im*(s).c2.re; \
  (r).c13.re=(u).c1.re*(s).c3.re-(u).c1.im*(s).c3.im; \
  (r).c13.im=(u).c1.re*(s).c3.im+(u).c1.im*(s).c3.re; \
  (r).c21.re=(u).c2.re*(s).c1.re-(u).c2.im*(s).c1.im; \
  (r).c21.im=(u).c2.re*(s).c1.im+(u).c2.im*(s).c1.re; \
  (r).c22.re=(u).c2.re*(s).c2.re-(u).c2.im*(s).c2.im; \
  (r).c22.im=(u).c2.re*(s).c2.im+(u).c2.im*(s).c2.re; \
  (r).c23.re=(u).c2.re*(s).c3.re-(u).c2.im*(s).c3.im; \
  (r).c23.im=(u).c2.re*(s).c3.im+(u).c2.im*(s).c3.re; \
  (r).c31.re=(u).c3.re*(s).c1.re-(u).c3.im*(s).c1.im; \
  (r).c31.im=(u).c3.re*(s).c1.im+(u).c3.im*(s).c1.re; \
  (r).c32.re=(u).c3.re*(s).c2.re-(u).c3.im*(s).c2.im; \
  (r).c32.im=(u).c3.re*(s).c2.im+(u).c3.im*(s).c2.re; \
  (r).c33.re=(u).c3.re*(s).c3.re-(u).c3.im*(s).c3.im; \
  (r).c33.im=(u).c3.re*(s).c3.im+(u).c3.im*(s).c3.re
             

static void mul_cplx(complex_dble *z,su3_vector_dble *in, su3_vector_dble *out)
{
   (*out).c1.re = (*z).re*(*in).c1.re-(*z).im*(*in).c1.im;
   (*out).c1.im = (*z).im*(*in).c1.re+(*z).re*(*in).c1.im;
   (*out).c2.re = (*z).re*(*in).c2.re-(*z).im*(*in).c2.im;
   (*out).c2.im = (*z).im*(*in).c2.re+(*z).re*(*in).c2.im;
   (*out).c3.re = (*z).re*(*in).c3.re-(*z).im*(*in).c3.im;
   (*out).c3.im = (*z).im*(*in).c3.re+(*z).re*(*in).c3.im;
}

static void peek_colour_column(su3_vector_dble *l, su3_dble *r, int column)
{
  switch (column)
  {
  case 1:
   (*l).c1.re=(*r).c11.re;(*l).c1.im=(*r).c11.im;
   (*l).c2.re=(*r).c21.re;(*l).c2.im=(*r).c21.im;
   (*l).c3.re=(*r).c31.re;(*l).c3.im=(*r).c31.im;
    break;
  case 2:               
   (*l).c1.re=(*r).c12.re;(*l).c1.im=(*r).c12.im;
   (*l).c2.re=(*r).c22.re;(*l).c2.im=(*r).c22.im;
   (*l).c3.re=(*r).c32.re;(*l).c3.im=(*r).c32.im;
    break;
  case 3:               
   (*l).c1.re=(*r).c13.re;(*l).c1.im=(*r).c13.im;
   (*l).c2.re=(*r).c23.re;(*l).c2.im=(*r).c23.im;
   (*l).c3.re=(*r).c33.re;(*l).c3.im=(*r).c33.im;
    break;

  }
}

static void poke_colour_column(su3_dble *r, su3_vector_dble *l, int column)
{
  switch (column)
  {
  case 1:
   (*r).c11.re=(*l).c1.re;(*r).c11.im=(*l).c1.im;
   (*r).c21.re=(*l).c2.re;(*r).c21.im=(*l).c2.im;
   (*r).c31.re=(*l).c3.re;(*r).c31.im=(*l).c3.im;
    break;
  case 2:
   (*r).c12.re=(*l).c1.re;(*r).c12.im=(*l).c1.im;
   (*r).c22.re=(*l).c2.re;(*r).c22.im=(*l).c2.im;
   (*r).c32.re=(*l).c3.re;(*r).c32.im=(*l).c3.im;
    break;
  case 3:
   (*r).c13.re=(*l).c1.re;(*r).c13.im=(*l).c1.im;
   (*r).c23.re=(*l).c2.re;(*r).c23.im=(*l).c2.im;
   (*r).c33.re=(*l).c3.re;(*r).c33.im=(*l).c3.im;
    break;

  }
}
extern void copy_sv_fs(int vol,spinor_dble *s,full_spinor_dble *sv,int id,int ic)
{
 int i=0;
 spinor_dble *rpk;
 rpk=s; 
 
 for (;rpk<(s+vol);)
 {
 
 switch (id)
  {
   case 1:
    poke_colour_column(&((*(sv+i)).c11),&((*rpk).c1),ic);
    poke_colour_column(&((*(sv+i)).c21),&((*rpk).c2),ic);
    poke_colour_column(&((*(sv+i)).c31),&((*rpk).c3),ic);
    poke_colour_column(&((*(sv+i)).c41),&((*rpk).c4),ic);
    break;
   case 2:                             
    poke_colour_column(&((*(sv+i)).c12),&((*rpk).c1),ic);
    poke_colour_column(&((*(sv+i)).c22),&((*rpk).c2),ic);
    poke_colour_column(&((*(sv+i)).c32),&((*rpk).c3),ic);
    poke_colour_column(&((*(sv+i)).c42),&((*rpk).c4),ic);
    break;
   case 3:                             
    poke_colour_column(&((*(sv+i)).c13),&((*rpk).c1),ic);
    poke_colour_column(&((*(sv+i)).c23),&((*rpk).c2),ic);
    poke_colour_column(&((*(sv+i)).c33),&((*rpk).c3),ic);
    poke_colour_column(&((*(sv+i)).c43),&((*rpk).c4),ic);
    break;
   case 4:                             
    poke_colour_column(&((*(sv+i)).c14),&((*rpk).c1),ic);
    poke_colour_column(&((*(sv+i)).c24),&((*rpk).c2),ic);
    poke_colour_column(&((*(sv+i)).c34),&((*rpk).c3),ic);
    poke_colour_column(&((*(sv+i)).c44),&((*rpk).c4),ic);
    break;
  }
 rpk+=1;
 i+=1;
 }
}
static const full_spinor_dble fsd0={{{0.0}}};

extern void set_fsd2zero(int vol,full_spinor_dble *pk)
{
   full_spinor_dble *rpk;

   for (rpk=pk;rpk<(pk+vol);rpk++)
      *rpk=fsd0;
}

extern void copy_fs_sv(int vol,full_spinor_dble *sv,spinor_dble *s,int id,int ic)
{
 int i=0;
 spinor_dble *rpk;
 rpk=s; 
 
 for (;rpk<(s+vol);)
 {
 
 switch (id)
  {
   case 1:

    peek_colour_column(&((*rpk).c1),&((*(sv+i)).c11),ic);
    peek_colour_column(&((*rpk).c2),&((*(sv+i)).c21),ic);
    peek_colour_column(&((*rpk).c3),&((*(sv+i)).c31),ic);
    peek_colour_column(&((*rpk).c4),&((*(sv+i)).c41),ic);
    break;
   case 2:                                       
    peek_colour_column(&((*rpk).c1),&((*(sv+i)).c12),ic);
    peek_colour_column(&((*rpk).c2),&((*(sv+i)).c22),ic);
    peek_colour_column(&((*rpk).c3),&((*(sv+i)).c32),ic);
    peek_colour_column(&((*rpk).c4),&((*(sv+i)).c42),ic);
    break;
   case 3:                                       
    peek_colour_column(&((*rpk).c1),&((*(sv+i)).c13),ic);
    peek_colour_column(&((*rpk).c2),&((*(sv+i)).c23),ic);
    peek_colour_column(&((*rpk).c3),&((*(sv+i)).c33),ic);
    peek_colour_column(&((*rpk).c4),&((*(sv+i)).c43),ic);
    break;
   case 4:                                       
    peek_colour_column(&((*rpk).c1),&((*(sv+i)).c14),ic);
    peek_colour_column(&((*rpk).c2),&((*(sv+i)).c24),ic);
    peek_colour_column(&((*rpk).c3),&((*(sv+i)).c34),ic);
    peek_colour_column(&((*rpk).c4),&((*(sv+i)).c44),ic);
    break;
  }
 rpk+=1;
 i+=1;
 }
}


static void mul_cplx_spinor_dble(complex_dble *z,spinor_dble *in,spinor_dble *out)
{
   mul_cplx(z,&((*in).c1),&((*out).c1));
   mul_cplx(z,&((*in).c2),&((*out).c2));
   mul_cplx(z,&((*in).c3),&((*out).c3));
   mul_cplx(z,&((*in).c4),&((*out).c4));
}


extern void adj_full_spinor(full_spinor_dble *in, full_spinor_dble *out)
{ 
  _su3_dagger((*out).c11, (*in).c11); 
  _su3_dagger((*out).c12, (*in).c21); 
  _su3_dagger((*out).c13, (*in).c31); 
  _su3_dagger((*out).c14, (*in).c41); 
  _su3_dagger((*out).c21, (*in).c12); 
  _su3_dagger((*out).c22, (*in).c22); 
  _su3_dagger((*out).c23, (*in).c32); 
  _su3_dagger((*out).c24, (*in).c42); 
  _su3_dagger((*out).c31, (*in).c13); 
  _su3_dagger((*out).c32, (*in).c23); 
  _su3_dagger((*out).c33, (*in).c33); 
  _su3_dagger((*out).c34, (*in).c43); 
  _su3_dagger((*out).c41, (*in).c14); 
  _su3_dagger((*out).c42, (*in).c24); 
  _su3_dagger((*out).c43, (*in).c34); 
  _su3_dagger((*out).c44, (*in).c44); 
}

#define _su3_times_su3_diag(u,v,w) \
   (u).c11.re= (v).c11.re*(w).c11.re-(v).c11.im*(w).c11.im  \
              +(v).c12.re*(w).c21.re-(v).c12.im*(w).c21.im  \
              +(v).c13.re*(w).c31.re-(v).c13.im*(w).c31.im; \
   (u).c11.im= (v).c11.re*(w).c11.im+(v).c11.im*(w).c11.re  \
              +(v).c12.re*(w).c21.im+(v).c12.im*(w).c21.re  \
              +(v).c13.re*(w).c31.im+(v).c13.im*(w).c31.re; \
   (u).c22.re= (v).c21.re*(w).c12.re-(v).c21.im*(w).c12.im  \
              +(v).c22.re*(w).c22.re-(v).c22.im*(w).c22.im  \
              +(v).c23.re*(w).c32.re-(v).c23.im*(w).c32.im; \
   (u).c22.im= (v).c21.re*(w).c12.im+(v).c21.im*(w).c12.re  \
              +(v).c22.re*(w).c22.im+(v).c22.im*(w).c22.re  \
              +(v).c23.re*(w).c32.im+(v).c23.im*(w).c32.re; \
   (u).c33.re= (v).c31.re*(w).c13.re-(v).c31.im*(w).c13.im  \
              +(v).c32.re*(w).c23.re-(v).c32.im*(w).c23.im  \
              +(v).c33.re*(w).c33.re-(v).c33.im*(w).c33.im; \
   (u).c33.im= (v).c31.re*(w).c13.im+(v).c31.im*(w).c13.re  \
              +(v).c32.re*(w).c23.im+(v).c32.im*(w).c23.re  \
              +(v).c33.re*(w).c33.im+(v).c33.im*(w).c33.re


extern void full_mat_prod(full_spinor_dble *in1,full_spinor_dble *in2,full_spinor_dble *out)
{
  _su3_times_su3_diag((*out).c11,(*in1).c11,(*in2).c11);
  _su3_times_su3_diag((*out).c12,(*in1).c12,(*in2).c21);
  _su3_times_su3_diag((*out).c13,(*in1).c13,(*in2).c31);
  _su3_times_su3_diag((*out).c14,(*in1).c14,(*in2).c41);
  _su3_times_su3_diag((*out).c21,(*in1).c21,(*in2).c12);
  _su3_times_su3_diag((*out).c22,(*in1).c22,(*in2).c22);
  _su3_times_su3_diag((*out).c23,(*in1).c23,(*in2).c32);
  _su3_times_su3_diag((*out).c24,(*in1).c24,(*in2).c42);
  _su3_times_su3_diag((*out).c31,(*in1).c31,(*in2).c13);
  _su3_times_su3_diag((*out).c32,(*in1).c32,(*in2).c23);
  _su3_times_su3_diag((*out).c33,(*in1).c33,(*in2).c33);
  _su3_times_su3_diag((*out).c34,(*in1).c34,(*in2).c43);
  _su3_times_su3_diag((*out).c41,(*in1).c41,(*in2).c14);
  _su3_times_su3_diag((*out).c42,(*in1).c42,(*in2).c24);
  _su3_times_su3_diag((*out).c43,(*in1).c43,(*in2).c34);
  _su3_times_su3_diag((*out).c44,(*in1).c44,(*in2).c44);
}

extern void meson_trace(full_spinor_dble *in1,full_spinor_dble *in2,complex_dble *out)
{
   full_spinor_dble r;

  /* Compute the contributions to the trace */
   full_mat_prod(in1, in2, &r);
	
  /* sum up the trace */
   (*out).re = r.c11.c11.re+r.c11.c22.re+r.c11.c33.re+
	 r.c12.c11.re+r.c12.c22.re+r.c12.c33.re+
	 r.c13.c11.re+r.c13.c22.re+r.c13.c33.re+
	 r.c14.c11.re+r.c14.c22.re+r.c14.c33.re+
	 r.c21.c11.re+r.c21.c22.re+r.c21.c33.re+
	 r.c22.c11.re+r.c22.c22.re+r.c22.c33.re+
 	 r.c23.c11.re+r.c23.c22.re+r.c23.c33.re+
 	 r.c24.c11.re+r.c24.c22.re+r.c24.c33.re+
 	 r.c31.c11.re+r.c31.c22.re+r.c31.c33.re+
 	 r.c32.c11.re+r.c32.c22.re+r.c32.c33.re+
 	 r.c33.c11.re+r.c33.c22.re+r.c33.c33.re+
 	 r.c34.c11.re+r.c34.c22.re+r.c34.c33.re+
 	 r.c41.c11.re+r.c41.c22.re+r.c41.c33.re+
 	 r.c42.c11.re+r.c42.c22.re+r.c42.c33.re+
 	 r.c43.c11.re+r.c43.c22.re+r.c43.c33.re+
 	 r.c44.c11.re+r.c44.c22.re+r.c44.c33.re;
 
   (*out).im = r.c11.c11.im+r.c11.c22.im+r.c11.c33.im+
 	 r.c12.c11.im+r.c12.c22.im+r.c12.c33.im+
 	 r.c13.c11.im+r.c13.c22.im+r.c13.c33.im+
	 r.c14.c11.im+r.c14.c22.im+r.c14.c33.im+
	 r.c21.c11.im+r.c21.c22.im+r.c21.c33.im+
	 r.c22.c11.im+r.c22.c22.im+r.c22.c33.im+
	 r.c23.c11.im+r.c23.c22.im+r.c23.c33.im+
	 r.c24.c11.im+r.c24.c22.im+r.c24.c33.im+
	 r.c31.c11.im+r.c31.c22.im+r.c31.c33.im+
	 r.c32.c11.im+r.c32.c22.im+r.c32.c33.im+
	 r.c33.c11.im+r.c33.c22.im+r.c33.c33.im+
	 r.c34.c11.im+r.c34.c22.im+r.c34.c33.im+
	 r.c41.c11.im+r.c41.c22.im+r.c41.c33.im+
	 r.c42.c11.im+r.c42.c22.im+r.c42.c33.im+
	 r.c43.c11.im+r.c43.c22.im+r.c43.c33.im+
	 r.c44.c11.im+r.c44.c22.im+r.c44.c33.im;
}

extern void fsv_trace(full_spinor_dble *in,complex_dble *out)
{
 /* sum up the trace */
 (*out).re = (*in).c11.c11.re+(*in).c11.c22.re+(*in).c11.c33.re+
	(*in).c22.c11.re+(*in).c22.c22.re+(*in).c22.c33.re+
	(*in).c33.c11.re+(*in).c33.c22.re+(*in).c33.c33.re+
	(*in).c44.c11.re+(*in).c44.c22.re+(*in).c44.c33.re;

 (*out).im = (*in).c11.c11.im+(*in).c11.c22.im+(*in).c11.c33.im+
	(*in).c22.c11.im+(*in).c22.c22.im+(*in).c22.c33.im+
	(*in).c33.c11.im+(*in).c33.c22.im+(*in).c33.c33.im+
	(*in).c44.c11.im+(*in).c44.c22.im+(*in).c44.c33.im;
}


extern void mul_cplx_su3_dble(complex_dble *z,su3_dble *in, su3_dble *out)
{
   (*out).c11.re=(*z).re*(*in).c11.re-(*z).im*(*in).c11.im;
   (*out).c12.re=(*z).re*(*in).c12.re-(*z).im*(*in).c12.im;
   (*out).c13.re=(*z).re*(*in).c13.re-(*z).im*(*in).c13.im;
   (*out).c21.re=(*z).re*(*in).c21.re-(*z).im*(*in).c21.im;
   (*out).c22.re=(*z).re*(*in).c22.re-(*z).im*(*in).c22.im;
   (*out).c23.re=(*z).re*(*in).c23.re-(*z).im*(*in).c23.im;
   (*out).c31.re=(*z).re*(*in).c31.re-(*z).im*(*in).c31.im;
   (*out).c32.re=(*z).re*(*in).c32.re-(*z).im*(*in).c32.im;
   (*out).c33.re=(*z).re*(*in).c33.re-(*z).im*(*in).c33.im;

   (*out).c11.im=(*z).re*(*in).c11.im+(*z).im*(*in).c11.re;
   (*out).c12.im=(*z).re*(*in).c12.im+(*z).im*(*in).c12.re;
   (*out).c13.im=(*z).re*(*in).c13.im+(*z).im*(*in).c13.re;
   (*out).c21.im=(*z).re*(*in).c21.im+(*z).im*(*in).c21.re;
   (*out).c22.im=(*z).re*(*in).c22.im+(*z).im*(*in).c22.re;
   (*out).c23.im=(*z).re*(*in).c23.im+(*z).im*(*in).c23.re;
   (*out).c31.im=(*z).re*(*in).c31.im+(*z).im*(*in).c31.re;
   (*out).c32.im=(*z).re*(*in).c32.im+(*z).im*(*in).c32.re;
   (*out).c33.im=(*z).re*(*in).c33.im+(*z).im*(*in).c33.re;
}


extern void mul_cmplx_fsv_dble(complex_dble *phase,full_spinor_dble *in,full_spinor_dble *out)
{
      mul_cplx_su3_dble(phase,&((*in).c11),&((*out).c11));
      mul_cplx_su3_dble(phase,&((*in).c12),&((*out).c12));
      mul_cplx_su3_dble(phase,&((*in).c13),&((*out).c13));
      mul_cplx_su3_dble(phase,&((*in).c14),&((*out).c14));
      mul_cplx_su3_dble(phase,&((*in).c21),&((*out).c21));
      mul_cplx_su3_dble(phase,&((*in).c22),&((*out).c22));
      mul_cplx_su3_dble(phase,&((*in).c23),&((*out).c23));
      mul_cplx_su3_dble(phase,&((*in).c24),&((*out).c24));
      mul_cplx_su3_dble(phase,&((*in).c31),&((*out).c31));
      mul_cplx_su3_dble(phase,&((*in).c32),&((*out).c32));
      mul_cplx_su3_dble(phase,&((*in).c33),&((*out).c33));
      mul_cplx_su3_dble(phase,&((*in).c34),&((*out).c34));
      mul_cplx_su3_dble(phase,&((*in).c41),&((*out).c41));
      mul_cplx_su3_dble(phase,&((*in).c42),&((*out).c42));
      mul_cplx_su3_dble(phase,&((*in).c43),&((*out).c43));
      mul_cplx_su3_dble(phase,&((*in).c44),&((*out).c44));
}


extern void mul_gamma_r(int mu,full_spinor_dble *in, full_spinor_dble *out)
{
   complex_dble i,m_i,m_1,p_1;
   
   i.re=0.0;
   i.im=1.0;
   
   m_i.re=0.0;
   m_i.im=-1.0;
   
   m_1.re=-1.0;
   m_1.im=0.0;

   p_1.re=1.0;
   p_1.im=0.0;

   /* The following encoding for the gamma matrices has been used here: 
      0 -> gamma_0 (time-component)
      5 -> gamma_5
     */
   if (mu==0)
   {
      mul_cplx_su3_dble(&m_1,&((*in).c13),&((*out).c11));
      mul_cplx_su3_dble(&m_1,&((*in).c14),&((*out).c12));
      mul_cplx_su3_dble(&m_1,&((*in).c11),&((*out).c13));
      mul_cplx_su3_dble(&m_1,&((*in).c12),&((*out).c14));
      mul_cplx_su3_dble(&m_1,&((*in).c23),&((*out).c21));
      mul_cplx_su3_dble(&m_1,&((*in).c24),&((*out).c22));
      mul_cplx_su3_dble(&m_1,&((*in).c21),&((*out).c23));
      mul_cplx_su3_dble(&m_1,&((*in).c22),&((*out).c24));
      mul_cplx_su3_dble(&m_1,&((*in).c33),&((*out).c31));
      mul_cplx_su3_dble(&m_1,&((*in).c34),&((*out).c32));
      mul_cplx_su3_dble(&m_1,&((*in).c31),&((*out).c33));
      mul_cplx_su3_dble(&m_1,&((*in).c32),&((*out).c34));
      mul_cplx_su3_dble(&m_1,&((*in).c43),&((*out).c41));
      mul_cplx_su3_dble(&m_1,&((*in).c44),&((*out).c42));
      mul_cplx_su3_dble(&m_1,&((*in).c41),&((*out).c43));
      mul_cplx_su3_dble(&m_1,&((*in).c42),&((*out).c44));
   }
   else if (mu==5)
   {
      mul_cplx_su3_dble(&p_1,&((*in).c11),&((*out).c11));
      mul_cplx_su3_dble(&p_1,&((*in).c12),&((*out).c12));
      mul_cplx_su3_dble(&m_1,&((*in).c13),&((*out).c13));
      mul_cplx_su3_dble(&m_1,&((*in).c14),&((*out).c14));
      mul_cplx_su3_dble(&p_1,&((*in).c21),&((*out).c21));
      mul_cplx_su3_dble(&p_1,&((*in).c22),&((*out).c22));
      mul_cplx_su3_dble(&m_1,&((*in).c23),&((*out).c23));
      mul_cplx_su3_dble(&m_1,&((*in).c24),&((*out).c24));
      mul_cplx_su3_dble(&p_1,&((*in).c31),&((*out).c31));
      mul_cplx_su3_dble(&p_1,&((*in).c32),&((*out).c32));
      mul_cplx_su3_dble(&m_1,&((*in).c33),&((*out).c33));
      mul_cplx_su3_dble(&m_1,&((*in).c34),&((*out).c34));
      mul_cplx_su3_dble(&p_1,&((*in).c41),&((*out).c41));
      mul_cplx_su3_dble(&p_1,&((*in).c42),&((*out).c42));
      mul_cplx_su3_dble(&m_1,&((*in).c43),&((*out).c43));
      mul_cplx_su3_dble(&m_1,&((*in).c44),&((*out).c44));
   }
   else
   {
      error_root(1,1,"mul_gamma_l [algebra.c]",
                 "Gamma matrix %d not defined ",mu);
   }
}


extern void mul_gamma_l(int mu,full_spinor_dble *in, full_spinor_dble *out)
{
   complex_dble i,m_i,m_1,p_1;

   i.re=0.0;
   i.im=1.0;

   m_i.re=0.0;
   m_i.im=-1.0;

   m_1.re=-1.0;
   m_1.im=0.0;

   p_1.re=1.0;
   p_1.im=0.0;

   /* The following encoding for the gamma matrices has been used here: 
      0 -> gamma_0 (time-component)
      1 -> gamma_1
      2 -> gamma_2
      3 -> gamma_3
      4 -> unity
      5 -> gamma_5
      6 -> gamma_0 gamma_5
      7 -> gamma_1 gamma_5
      8 -> gamma_2 gamma_5
      9 -> gamma_3 gamma_5
      10-> gamma_0 gamma_1
      11-> gamma_0 gamma_2
      12-> gamma_0 gamma_3
      13-> gamma_1 gamma_2
      14-> gamma_1 gamma_3
      15-> gamma_2 gamma_3
     */
   if (mu==0)
   {
      mul_cplx_su3_dble(&m_1,&((*in).c31),&((*out).c11));
      mul_cplx_su3_dble(&m_1,&((*in).c32),&((*out).c12));
      mul_cplx_su3_dble(&m_1,&((*in).c33),&((*out).c13));
      mul_cplx_su3_dble(&m_1,&((*in).c34),&((*out).c14));
      mul_cplx_su3_dble(&m_1,&((*in).c41),&((*out).c21));
      mul_cplx_su3_dble(&m_1,&((*in).c42),&((*out).c22));
      mul_cplx_su3_dble(&m_1,&((*in).c43),&((*out).c23));
      mul_cplx_su3_dble(&m_1,&((*in).c44),&((*out).c24));
      mul_cplx_su3_dble(&m_1,&((*in).c11),&((*out).c31));
      mul_cplx_su3_dble(&m_1,&((*in).c12),&((*out).c32));
      mul_cplx_su3_dble(&m_1,&((*in).c13),&((*out).c33));
      mul_cplx_su3_dble(&m_1,&((*in).c14),&((*out).c34));
      mul_cplx_su3_dble(&m_1,&((*in).c21),&((*out).c41));
      mul_cplx_su3_dble(&m_1,&((*in).c22),&((*out).c42));
      mul_cplx_su3_dble(&m_1,&((*in).c23),&((*out).c43));
      mul_cplx_su3_dble(&m_1,&((*in).c24),&((*out).c44));
   }
   else if (mu==1)
   {
      mul_cplx_su3_dble(&m_i,&((*in).c41),&((*out).c11));
      mul_cplx_su3_dble(&m_i,&((*in).c42),&((*out).c12));
      mul_cplx_su3_dble(&m_i,&((*in).c43),&((*out).c13));
      mul_cplx_su3_dble(&m_i,&((*in).c44),&((*out).c14));
      mul_cplx_su3_dble(&m_i,&((*in).c31),&((*out).c21));
      mul_cplx_su3_dble(&m_i,&((*in).c32),&((*out).c22));
      mul_cplx_su3_dble(&m_i,&((*in).c33),&((*out).c23));
      mul_cplx_su3_dble(&m_i,&((*in).c34),&((*out).c24));
      mul_cplx_su3_dble(  &i,&((*in).c21),&((*out).c31));
      mul_cplx_su3_dble(  &i,&((*in).c22),&((*out).c32));
      mul_cplx_su3_dble(  &i,&((*in).c23),&((*out).c33));
      mul_cplx_su3_dble(  &i,&((*in).c24),&((*out).c34));
      mul_cplx_su3_dble(  &i,&((*in).c11),&((*out).c41));
      mul_cplx_su3_dble(  &i,&((*in).c12),&((*out).c42));
      mul_cplx_su3_dble(  &i,&((*in).c13),&((*out).c43));
      mul_cplx_su3_dble(  &i,&((*in).c14),&((*out).c44));
   }
   else if (mu==2) 
   {
      mul_cplx_su3_dble(&m_1,&((*in).c41),&((*out).c11));
      mul_cplx_su3_dble(&m_1,&((*in).c42),&((*out).c12));
      mul_cplx_su3_dble(&m_1,&((*in).c43),&((*out).c13));
      mul_cplx_su3_dble(&m_1,&((*in).c44),&((*out).c14));
      mul_cplx_su3_dble(&p_1,&((*in).c31),&((*out).c21));
      mul_cplx_su3_dble(&p_1,&((*in).c32),&((*out).c22));
      mul_cplx_su3_dble(&p_1,&((*in).c33),&((*out).c23));
      mul_cplx_su3_dble(&p_1,&((*in).c34),&((*out).c24));
      mul_cplx_su3_dble(&p_1,&((*in).c21),&((*out).c31));
      mul_cplx_su3_dble(&p_1,&((*in).c22),&((*out).c32));
      mul_cplx_su3_dble(&p_1,&((*in).c23),&((*out).c33));
      mul_cplx_su3_dble(&p_1,&((*in).c24),&((*out).c34));
      mul_cplx_su3_dble(&m_1,&((*in).c11),&((*out).c41));
      mul_cplx_su3_dble(&m_1,&((*in).c12),&((*out).c42));
      mul_cplx_su3_dble(&m_1,&((*in).c13),&((*out).c43));
      mul_cplx_su3_dble(&m_1,&((*in).c14),&((*out).c44));
   }
   else if (mu==3)
   {
     mul_cplx_su3_dble(&m_i,&((*in).c31),&((*out).c11));
     mul_cplx_su3_dble(&m_i,&((*in).c32),&((*out).c12));
     mul_cplx_su3_dble(&m_i,&((*in).c33),&((*out).c13));
     mul_cplx_su3_dble(&m_i,&((*in).c34),&((*out).c14));
     mul_cplx_su3_dble(  &i,&((*in).c41),&((*out).c21));
     mul_cplx_su3_dble(  &i,&((*in).c42),&((*out).c22));
     mul_cplx_su3_dble(  &i,&((*in).c43),&((*out).c23));
     mul_cplx_su3_dble(  &i,&((*in).c44),&((*out).c24));
     mul_cplx_su3_dble(  &i,&((*in).c11),&((*out).c31));
     mul_cplx_su3_dble(  &i,&((*in).c12),&((*out).c32));
     mul_cplx_su3_dble(  &i,&((*in).c13),&((*out).c33));
     mul_cplx_su3_dble(  &i,&((*in).c14),&((*out).c34));
     mul_cplx_su3_dble(&m_i,&((*in).c21),&((*out).c41));
     mul_cplx_su3_dble(&m_i,&((*in).c22),&((*out).c42));
     mul_cplx_su3_dble(&m_i,&((*in).c23),&((*out).c43));
     mul_cplx_su3_dble(&m_i,&((*in).c24),&((*out).c44));
   }
   else if (mu==4)
   {
      mul_cplx_su3_dble(&p_1,&((*in).c11),&((*out).c11));
      mul_cplx_su3_dble(&p_1,&((*in).c12),&((*out).c12));
      mul_cplx_su3_dble(&p_1,&((*in).c13),&((*out).c13));
      mul_cplx_su3_dble(&p_1,&((*in).c14),&((*out).c14));
      mul_cplx_su3_dble(&p_1,&((*in).c21),&((*out).c21));
      mul_cplx_su3_dble(&p_1,&((*in).c22),&((*out).c22));
      mul_cplx_su3_dble(&p_1,&((*in).c23),&((*out).c23));
      mul_cplx_su3_dble(&p_1,&((*in).c24),&((*out).c24));
      mul_cplx_su3_dble(&p_1,&((*in).c31),&((*out).c31));
      mul_cplx_su3_dble(&p_1,&((*in).c32),&((*out).c32));
      mul_cplx_su3_dble(&p_1,&((*in).c33),&((*out).c33));
      mul_cplx_su3_dble(&p_1,&((*in).c34),&((*out).c34));
      mul_cplx_su3_dble(&p_1,&((*in).c41),&((*out).c41));
      mul_cplx_su3_dble(&p_1,&((*in).c42),&((*out).c42));
      mul_cplx_su3_dble(&p_1,&((*in).c43),&((*out).c43));
      mul_cplx_su3_dble(&p_1,&((*in).c44),&((*out).c44));
   }
   else if (mu==5)
   {
      mul_cplx_su3_dble(&p_1,&((*in).c11),&((*out).c11));
      mul_cplx_su3_dble(&p_1,&((*in).c12),&((*out).c12));
      mul_cplx_su3_dble(&p_1,&((*in).c13),&((*out).c13));
      mul_cplx_su3_dble(&p_1,&((*in).c14),&((*out).c14));
      mul_cplx_su3_dble(&p_1,&((*in).c21),&((*out).c21));
      mul_cplx_su3_dble(&p_1,&((*in).c22),&((*out).c22));
      mul_cplx_su3_dble(&p_1,&((*in).c23),&((*out).c23));
      mul_cplx_su3_dble(&p_1,&((*in).c24),&((*out).c24));
      mul_cplx_su3_dble(&m_1,&((*in).c31),&((*out).c31));
      mul_cplx_su3_dble(&m_1,&((*in).c32),&((*out).c32));
      mul_cplx_su3_dble(&m_1,&((*in).c33),&((*out).c33));
      mul_cplx_su3_dble(&m_1,&((*in).c34),&((*out).c34));
      mul_cplx_su3_dble(&m_1,&((*in).c41),&((*out).c41));
      mul_cplx_su3_dble(&m_1,&((*in).c42),&((*out).c42));
      mul_cplx_su3_dble(&m_1,&((*in).c43),&((*out).c43));
      mul_cplx_su3_dble(&m_1,&((*in).c44),&((*out).c44));
   }
   else if (mu==6) /* gamma_0 gamma_5*/
   {
      mul_cplx_su3_dble(&m_1,&((*in).c31),&((*out).c11));
      mul_cplx_su3_dble(&m_1,&((*in).c32),&((*out).c12));
      mul_cplx_su3_dble(&m_1,&((*in).c33),&((*out).c13));
      mul_cplx_su3_dble(&m_1,&((*in).c34),&((*out).c14));
      mul_cplx_su3_dble(&m_1,&((*in).c41),&((*out).c21));
      mul_cplx_su3_dble(&m_1,&((*in).c42),&((*out).c22));
      mul_cplx_su3_dble(&m_1,&((*in).c43),&((*out).c23));
      mul_cplx_su3_dble(&m_1,&((*in).c44),&((*out).c24));
      mul_cplx_su3_dble(&p_1,&((*in).c11),&((*out).c31));
      mul_cplx_su3_dble(&p_1,&((*in).c12),&((*out).c32));
      mul_cplx_su3_dble(&p_1,&((*in).c13),&((*out).c33));
      mul_cplx_su3_dble(&p_1,&((*in).c14),&((*out).c34));
      mul_cplx_su3_dble(&p_1,&((*in).c21),&((*out).c41));
      mul_cplx_su3_dble(&p_1,&((*in).c22),&((*out).c42));
      mul_cplx_su3_dble(&p_1,&((*in).c23),&((*out).c43));
      mul_cplx_su3_dble(&p_1,&((*in).c24),&((*out).c44));
   }
   else if (mu==7) /* gamma_1 gamma_5 */
   {
      mul_cplx_su3_dble(&m_i,&((*in).c41),&((*out).c11));
      mul_cplx_su3_dble(&m_i,&((*in).c42),&((*out).c12));
      mul_cplx_su3_dble(&m_i,&((*in).c43),&((*out).c13));
      mul_cplx_su3_dble(&m_i,&((*in).c44),&((*out).c14));
      mul_cplx_su3_dble(&m_i,&((*in).c31),&((*out).c21));
      mul_cplx_su3_dble(&m_i,&((*in).c32),&((*out).c22));
      mul_cplx_su3_dble(&m_i,&((*in).c33),&((*out).c23));
      mul_cplx_su3_dble(&m_i,&((*in).c34),&((*out).c24));
      mul_cplx_su3_dble(&m_i,&((*in).c21),&((*out).c31));
      mul_cplx_su3_dble(&m_i,&((*in).c22),&((*out).c32));
      mul_cplx_su3_dble(&m_i,&((*in).c23),&((*out).c33));
      mul_cplx_su3_dble(&m_i,&((*in).c24),&((*out).c34));
      mul_cplx_su3_dble(&m_i,&((*in).c11),&((*out).c41));
      mul_cplx_su3_dble(&m_i,&((*in).c12),&((*out).c42));
      mul_cplx_su3_dble(&m_i,&((*in).c13),&((*out).c43));
      mul_cplx_su3_dble(&m_i,&((*in).c14),&((*out).c44));
   }
   else if (mu==8) /* gamma_2 gamma_5 */
   {
      mul_cplx_su3_dble(&m_1,&((*in).c41),&((*out).c11));
      mul_cplx_su3_dble(&m_1,&((*in).c42),&((*out).c12));
      mul_cplx_su3_dble(&m_1,&((*in).c43),&((*out).c13));
      mul_cplx_su3_dble(&m_1,&((*in).c44),&((*out).c14));
      mul_cplx_su3_dble(&p_1,&((*in).c31),&((*out).c21));
      mul_cplx_su3_dble(&p_1,&((*in).c32),&((*out).c22));
      mul_cplx_su3_dble(&p_1,&((*in).c33),&((*out).c23));
      mul_cplx_su3_dble(&p_1,&((*in).c34),&((*out).c24));
      mul_cplx_su3_dble(&m_1,&((*in).c21),&((*out).c31));
      mul_cplx_su3_dble(&m_1,&((*in).c22),&((*out).c32));
      mul_cplx_su3_dble(&m_1,&((*in).c23),&((*out).c33));
      mul_cplx_su3_dble(&m_1,&((*in).c24),&((*out).c34));
      mul_cplx_su3_dble(&p_1,&((*in).c11),&((*out).c41));
      mul_cplx_su3_dble(&p_1,&((*in).c12),&((*out).c42));
      mul_cplx_su3_dble(&p_1,&((*in).c13),&((*out).c43));
      mul_cplx_su3_dble(&p_1,&((*in).c14),&((*out).c44));
   }
   else if (mu==9) /* gamma_3 gamma_5 */
   {
      mul_cplx_su3_dble(&m_i,&((*in).c31),&((*out).c11));
      mul_cplx_su3_dble(&m_i,&((*in).c32),&((*out).c12));
      mul_cplx_su3_dble(&m_i,&((*in).c33),&((*out).c13));
      mul_cplx_su3_dble(&m_i,&((*in).c34),&((*out).c14));
      mul_cplx_su3_dble(  &i,&((*in).c41),&((*out).c21));
      mul_cplx_su3_dble(  &i,&((*in).c42),&((*out).c22));
      mul_cplx_su3_dble(  &i,&((*in).c43),&((*out).c23));
      mul_cplx_su3_dble(  &i,&((*in).c44),&((*out).c24));
      mul_cplx_su3_dble(&m_i,&((*in).c11),&((*out).c31));
      mul_cplx_su3_dble(&m_i,&((*in).c12),&((*out).c32));
      mul_cplx_su3_dble(&m_i,&((*in).c13),&((*out).c33));
      mul_cplx_su3_dble(&m_i,&((*in).c14),&((*out).c34));
      mul_cplx_su3_dble(  &i,&((*in).c21),&((*out).c41));
      mul_cplx_su3_dble(  &i,&((*in).c22),&((*out).c42));
      mul_cplx_su3_dble(  &i,&((*in).c23),&((*out).c43));
      mul_cplx_su3_dble(  &i,&((*in).c24),&((*out).c44));
   }
   else if (mu==10) /* gamma_0 gamma_1 */
   {
      mul_cplx_su3_dble(&m_i,&((*in).c21),&((*out).c11));
      mul_cplx_su3_dble(&m_i,&((*in).c22),&((*out).c12));
      mul_cplx_su3_dble(&m_i,&((*in).c23),&((*out).c13));
      mul_cplx_su3_dble(&m_i,&((*in).c24),&((*out).c14));
      mul_cplx_su3_dble(&m_i,&((*in).c11),&((*out).c21));
      mul_cplx_su3_dble(&m_i,&((*in).c12),&((*out).c22));
      mul_cplx_su3_dble(&m_i,&((*in).c13),&((*out).c23));
      mul_cplx_su3_dble(&m_i,&((*in).c14),&((*out).c24));
      mul_cplx_su3_dble(  &i,&((*in).c41),&((*out).c31));
      mul_cplx_su3_dble(  &i,&((*in).c42),&((*out).c32));
      mul_cplx_su3_dble(  &i,&((*in).c43),&((*out).c33));
      mul_cplx_su3_dble(  &i,&((*in).c44),&((*out).c34));
      mul_cplx_su3_dble(  &i,&((*in).c31),&((*out).c41));
      mul_cplx_su3_dble(  &i,&((*in).c32),&((*out).c42));
      mul_cplx_su3_dble(  &i,&((*in).c33),&((*out).c43));
      mul_cplx_su3_dble(  &i,&((*in).c34),&((*out).c44));
   }
   else if (mu==11) /* gamam_0 gamma_2 */
   {
      mul_cplx_su3_dble(&m_1,&((*in).c21),&((*out).c11));
      mul_cplx_su3_dble(&m_1,&((*in).c22),&((*out).c12));
      mul_cplx_su3_dble(&m_1,&((*in).c23),&((*out).c13));
      mul_cplx_su3_dble(&m_1,&((*in).c24),&((*out).c14));
      mul_cplx_su3_dble(&p_1,&((*in).c11),&((*out).c21));
      mul_cplx_su3_dble(&p_1,&((*in).c12),&((*out).c22));
      mul_cplx_su3_dble(&p_1,&((*in).c13),&((*out).c23));
      mul_cplx_su3_dble(&p_1,&((*in).c14),&((*out).c24));
      mul_cplx_su3_dble(&p_1,&((*in).c41),&((*out).c31));
      mul_cplx_su3_dble(&p_1,&((*in).c42),&((*out).c32));
      mul_cplx_su3_dble(&p_1,&((*in).c43),&((*out).c33));
      mul_cplx_su3_dble(&p_1,&((*in).c44),&((*out).c34));
      mul_cplx_su3_dble(&m_1,&((*in).c31),&((*out).c41));
      mul_cplx_su3_dble(&m_1,&((*in).c32),&((*out).c42));
      mul_cplx_su3_dble(&m_1,&((*in).c33),&((*out).c43));
      mul_cplx_su3_dble(&m_1,&((*in).c34),&((*out).c44));
   }
   else if (mu==12) /* gamma_0 gamma_3 */
   {
      mul_cplx_su3_dble(&m_i,&((*in).c11),&((*out).c11));
      mul_cplx_su3_dble(&m_i,&((*in).c12),&((*out).c12));
      mul_cplx_su3_dble(&m_i,&((*in).c13),&((*out).c13));
      mul_cplx_su3_dble(&m_i,&((*in).c14),&((*out).c14));
      mul_cplx_su3_dble(  &i,&((*in).c21),&((*out).c21));
      mul_cplx_su3_dble(  &i,&((*in).c22),&((*out).c22));
      mul_cplx_su3_dble(  &i,&((*in).c23),&((*out).c23));
      mul_cplx_su3_dble(  &i,&((*in).c24),&((*out).c24));
      mul_cplx_su3_dble(  &i,&((*in).c31),&((*out).c31));
      mul_cplx_su3_dble(  &i,&((*in).c32),&((*out).c32));
      mul_cplx_su3_dble(  &i,&((*in).c33),&((*out).c33));
      mul_cplx_su3_dble(  &i,&((*in).c34),&((*out).c34));
      mul_cplx_su3_dble(&m_i,&((*in).c41),&((*out).c41));
      mul_cplx_su3_dble(&m_i,&((*in).c42),&((*out).c42));
      mul_cplx_su3_dble(&m_i,&((*in).c43),&((*out).c43));
      mul_cplx_su3_dble(&m_i,&((*in).c44),&((*out).c44));
   }
   else if (mu==13) /* gamma_1 gamma_2 */
   {
      mul_cplx_su3_dble(  &i,&((*in).c11),&((*out).c11));
      mul_cplx_su3_dble(  &i,&((*in).c12),&((*out).c12));
      mul_cplx_su3_dble(  &i,&((*in).c13),&((*out).c13));
      mul_cplx_su3_dble(  &i,&((*in).c14),&((*out).c14));
      mul_cplx_su3_dble(&m_i,&((*in).c21),&((*out).c21));
      mul_cplx_su3_dble(&m_i,&((*in).c22),&((*out).c22));
      mul_cplx_su3_dble(&m_i,&((*in).c23),&((*out).c23));
      mul_cplx_su3_dble(&m_i,&((*in).c24),&((*out).c24));
      mul_cplx_su3_dble(  &i,&((*in).c31),&((*out).c31));
      mul_cplx_su3_dble(  &i,&((*in).c32),&((*out).c32));
      mul_cplx_su3_dble(  &i,&((*in).c33),&((*out).c33));
      mul_cplx_su3_dble(  &i,&((*in).c34),&((*out).c34));
      mul_cplx_su3_dble(&m_i,&((*in).c41),&((*out).c41));
      mul_cplx_su3_dble(&m_i,&((*in).c42),&((*out).c42));
      mul_cplx_su3_dble(&m_i,&((*in).c43),&((*out).c43));
      mul_cplx_su3_dble(&m_i,&((*in).c44),&((*out).c44));
   }
   else if (mu==14) /* gamma_1 gamma_3 */ 
   {
      mul_cplx_su3_dble(&m_1,&((*in).c21),&((*out).c11));
      mul_cplx_su3_dble(&m_1,&((*in).c22),&((*out).c12));
      mul_cplx_su3_dble(&m_1,&((*in).c23),&((*out).c13));
      mul_cplx_su3_dble(&m_1,&((*in).c24),&((*out).c14));
      mul_cplx_su3_dble(&p_1,&((*in).c11),&((*out).c21));
      mul_cplx_su3_dble(&p_1,&((*in).c12),&((*out).c22));
      mul_cplx_su3_dble(&p_1,&((*in).c13),&((*out).c23));
      mul_cplx_su3_dble(&p_1,&((*in).c14),&((*out).c24));
      mul_cplx_su3_dble(&m_1,&((*in).c41),&((*out).c31));
      mul_cplx_su3_dble(&m_1,&((*in).c42),&((*out).c32));
      mul_cplx_su3_dble(&m_1,&((*in).c43),&((*out).c33));
      mul_cplx_su3_dble(&m_1,&((*in).c44),&((*out).c34));
      mul_cplx_su3_dble(&p_1,&((*in).c31),&((*out).c41));
      mul_cplx_su3_dble(&p_1,&((*in).c32),&((*out).c42));
      mul_cplx_su3_dble(&p_1,&((*in).c33),&((*out).c43));
      mul_cplx_su3_dble(&p_1,&((*in).c34),&((*out).c44));
   }
   else if (mu==15) /* gamma_2 gamma_3 */
   {
      mul_cplx_su3_dble(&i,&((*in).c21),&((*out).c11));
      mul_cplx_su3_dble(&i,&((*in).c22),&((*out).c12));
      mul_cplx_su3_dble(&i,&((*in).c23),&((*out).c13));
      mul_cplx_su3_dble(&i,&((*in).c24),&((*out).c14));
      mul_cplx_su3_dble(&i,&((*in).c11),&((*out).c21));
      mul_cplx_su3_dble(&i,&((*in).c12),&((*out).c22));
      mul_cplx_su3_dble(&i,&((*in).c13),&((*out).c23));
      mul_cplx_su3_dble(&i,&((*in).c14),&((*out).c24));
      mul_cplx_su3_dble(&i,&((*in).c41),&((*out).c31));
      mul_cplx_su3_dble(&i,&((*in).c42),&((*out).c32));
      mul_cplx_su3_dble(&i,&((*in).c43),&((*out).c33));
      mul_cplx_su3_dble(&i,&((*in).c44),&((*out).c34));
      mul_cplx_su3_dble(&i,&((*in).c31),&((*out).c41));
      mul_cplx_su3_dble(&i,&((*in).c32),&((*out).c42));
      mul_cplx_su3_dble(&i,&((*in).c33),&((*out).c43));
      mul_cplx_su3_dble(&i,&((*in).c34),&((*out).c44));
   }
   else
   {
      error_root(1,1,"mul_gamma_l [algebra.c]",
                 "Gamma matrix %d not defined ",mu);
   }
}


extern void mul_s_gamma_r(int mu,spinor_dble *in,spinor_dble *out)
{

   complex_dble i,m_i,m_1;

   i.re=0.0;
   i.im=1.0;

   m_i.re=0.0;
   m_i.im=-1.0;

   m_1.re=-1.0;
   m_1.im=0.0;

   if (mu==0)
   {
      mul_cplx(&m_1,&((*in).c3),&((*out).c1));
      mul_cplx(&m_1,&((*in).c4),&((*out).c2));
      mul_cplx(&m_1,&((*in).c1),&((*out).c3));
      mul_cplx(&m_1,&((*in).c2),&((*out).c4));
   }
   else if (mu==1)
   {
      mul_cplx(  &i,&((*in).c4),&((*out).c1));
      mul_cplx(  &i,&((*in).c3),&((*out).c2));
      mul_cplx(&m_i,&((*in).c2),&((*out).c3));
      mul_cplx(&m_i,&((*in).c1),&((*out).c4));
   }
   else if (mu==2)
   {
      mul_cplx(&m_1,&((*in).c4),&((*out).c1));
      (*out).c2=(*in).c3;
      (*out).c3=(*in).c2;
      mul_cplx(&m_1,&((*in).c1),&((*out).c4));
   }
   else if (mu==3)
   {
      mul_cplx(  &i,&((*in).c3),&((*out).c1));
      mul_cplx(&m_i,&((*in).c4),&((*out).c2));
      mul_cplx(&m_i,&((*in).c1),&((*out).c3));
      mul_cplx(  &i,&((*in).c2),&((*out).c4));
   }
   else if (mu==4)
   {
      (*out).c1=(*in).c1;
      (*out).c2=(*in).c2;
      (*out).c3=(*in).c3;
      (*out).c4=(*in).c4;
   }
   else if (mu==5)
   {
      (*out).c1=(*in).c1;
      (*out).c2=(*in).c2;
      mul_cplx(&m_1,&((*in).c3),&((*out).c3));
      mul_cplx(&m_1,&((*in).c4),&((*out).c4));
   }
   else
   {
      error_root(1,1,"mul_s_gamma_r [algebra.c]",
                 "Gamma matrix %d not defined ",mu);
   }

}
extern void mul_s_gamma_l(int mu,spinor_dble *in,spinor_dble *out)
{
   complex_dble i,m_i,m_1;

   i.re=0.0;
   i.im=1.0;

   m_i.re=0.0;
   m_i.im=-1.0;

   m_1.re=-1.0;
   m_1.im=0.0;

   if (mu==0)
   {
      mul_cplx(&m_1,&((*in).c3),&((*out).c1));
      mul_cplx(&m_1,&((*in).c4),&((*out).c2));
      mul_cplx(&m_1,&((*in).c1),&((*out).c3));
      mul_cplx(&m_1,&((*in).c2),&((*out).c4));
   }
   else if (mu==1)
   {
      mul_cplx(&m_i,&((*in).c4),&((*out).c1));
      mul_cplx(&m_i,&((*in).c3),&((*out).c2));
      mul_cplx(  &i,&((*in).c2),&((*out).c3));
      mul_cplx(  &i,&((*in).c1),&((*out).c4));
   }
   else if (mu==2)
   {
      mul_cplx(&m_1,&((*in).c4),&((*out).c1));
      (*out).c2=(*in).c3;
      (*out).c3=(*in).c2;
      mul_cplx(&m_1,&((*in).c1),&((*out).c4));
   }
   else if (mu==3)
   {
      mul_cplx(&m_i,&((*in).c3),&((*out).c1));
      mul_cplx(  &i,&((*in).c4),&((*out).c2));
      mul_cplx(  &i,&((*in).c1),&((*out).c3));
      mul_cplx(&m_i,&((*in).c2),&((*out).c4));
   }
   else if (mu==4)
   {
      (*out).c1=(*in).c1;
      (*out).c2=(*in).c2;
      (*out).c3=(*in).c3;
      (*out).c4=(*in).c4;
   }
   else if (mu==5)
   {
      (*out).c1=(*in).c1;
      (*out).c2=(*in).c2;
      mul_cplx(&m_1,&((*in).c3),&((*out).c3));
      mul_cplx(&m_1,&((*in).c4),&((*out).c4));
   }
   else if (mu==6)
   {
      mul_cplx(&m_1,&((*in).c3),&((*out).c1));
      mul_cplx(&m_1,&((*in).c4),&((*out).c2));
      (*out).c3= (*in).c1;
      (*out).c4= (*in).c2;
   }
   else if (mu==7)
   {
      mul_cplx(&m_i,&((*in).c4),&((*out).c1));
      mul_cplx(&m_i,&((*in).c3),&((*out).c2));
      mul_cplx(&m_i,&((*in).c2),&((*out).c3));
      mul_cplx(&m_i,&((*in).c1),&((*out).c4));
   }
   else if (mu==8)
   {
      mul_cplx(&m_1,&((*in).c4),&((*out).c1));
      (*out).c2= (*in).c3;
      mul_cplx(&m_1,&((*in).c2),&((*out).c3));
      (*out).c4= (*in).c1;
   }
   else if (mu==9)
   {
      mul_cplx(&m_i  ,&((*in).c3),&((*out).c1));
      mul_cplx(&i,&((*in).c4),&((*out).c2));
      mul_cplx(&m_i  ,&((*in).c1),&((*out).c3));
      mul_cplx(&i,&((*in).c2),&((*out).c4));
   }
   else if (mu==10)
   {
      mul_cplx(&m_i,&((*in).c2),&((*out).c1));
      mul_cplx(&m_i,&((*in).c1),&((*out).c2));
      mul_cplx(  &i,&((*in).c4),&((*out).c3));
      mul_cplx(  &i,&((*in).c3),&((*out).c4));
   }
   else if (mu==11)
   {
      mul_cplx(&m_1,&((*in).c2),&((*out).c1));
      (*out).c2= (*in).c1;
      (*out).c3= (*in).c4;
      mul_cplx(&m_1,&((*in).c3),&((*out).c4));
   }
   else if (mu==12)
   {
      mul_cplx(&m_i, &((*in).c1),&((*out).c1));
      mul_cplx(&i  , &((*in).c2),&((*out).c2));
      mul_cplx(&i  , &((*in).c3),&((*out).c3));
      mul_cplx(&m_i, &((*in).c4),&((*out).c4));
   }
   else if (mu==13)
   {
      mul_cplx(&i  ,&((*in).c1),&((*out).c1));
      mul_cplx(&m_i,&((*in).c2),&((*out).c2));
      mul_cplx(&i  ,&((*in).c3),&((*out).c3));
      mul_cplx(&m_i,&((*in).c4),&((*out).c4));
   }
   else if (mu==14)
   {
      mul_cplx(&m_1,&((*in).c2),&((*out).c1));
      (*out).c2= (*in).c1;
      mul_cplx(&m_1,&((*in).c4),&((*out).c3));
      (*out).c4= (*in).c3;
   }
   else if (mu==15)
   {
      mul_cplx(&i,&((*in).c2),&((*out).c1));
      mul_cplx(&i,&((*in).c1),&((*out).c2));
      mul_cplx(&i,&((*in).c4),&((*out).c3));
      mul_cplx(&i,&((*in).c3),&((*out).c4));
   }
   else
   {
      error_root(1,1,"mul_s_gamma_l [algebra.c]",
                 "Gamma matrix %d not defined ",mu);
   }
}

extern complex_dble spinor_prod_dble_local(spinor_dble *k,spinor_dble *l)
{
 complex_dble z;
        z.re=(_vector_prod_re((*k).c1,(*l).c1)+
              _vector_prod_re((*k).c2,(*l).c2)+
              _vector_prod_re((*k).c3,(*l).c3)+
              _vector_prod_re((*k).c4,(*l).c4));

        z.im=(_vector_prod_im((*k).c1,(*l).c1)+
              _vector_prod_im((*k).c2,(*l).c2)+
              _vector_prod_im((*k).c3,(*l).c3)+
              _vector_prod_im((*k).c4,(*l).c4));
 return z;
}
extern full_spinor_dble su3_times_fsd(su3_dble u, full_spinor_dble s){

   full_spinor_dble su;
   _su3_times_su3(su.c11,u,s.c11);
   _su3_times_su3(su.c12,u,s.c12);
   _su3_times_su3(su.c13,u,s.c13);
   _su3_times_su3(su.c14,u,s.c14);
   _su3_times_su3(su.c21,u,s.c21);
   _su3_times_su3(su.c22,u,s.c22);
   _su3_times_su3(su.c23,u,s.c23);
   _su3_times_su3(su.c24,u,s.c24);
   _su3_times_su3(su.c31,u,s.c31);
   _su3_times_su3(su.c32,u,s.c32);
   _su3_times_su3(su.c33,u,s.c33);
   _su3_times_su3(su.c34,u,s.c34);
   _su3_times_su3(su.c41,u,s.c41);
   _su3_times_su3(su.c42,u,s.c42);
   _su3_times_su3(su.c43,u,s.c43);
   _su3_times_su3(su.c44,u,s.c44);
   
   return su;
}

extern full_spinor_dble fsd_times_su3(su3_dble u, full_spinor_dble s){

   full_spinor_dble su;
   _su3_times_su3(su.c11,s.c11,u);
   _su3_times_su3(su.c12,s.c12,u);
   _su3_times_su3(su.c13,s.c13,u);
   _su3_times_su3(su.c14,s.c14,u);
   _su3_times_su3(su.c21,s.c21,u);
   _su3_times_su3(su.c22,s.c22,u);
   _su3_times_su3(su.c23,s.c23,u);
   _su3_times_su3(su.c24,s.c24,u);
   _su3_times_su3(su.c31,s.c31,u);
   _su3_times_su3(su.c32,s.c32,u);
   _su3_times_su3(su.c33,s.c33,u);
   _su3_times_su3(su.c34,s.c34,u);
   _su3_times_su3(su.c41,s.c41,u);
   _su3_times_su3(su.c42,s.c42,u);
   _su3_times_su3(su.c43,s.c43,u);
   _su3_times_su3(su.c44,s.c44,u);

   return su;
}
extern su3_dble mul_cplx_su3_dble_out(complex_dble z,su3_dble s)
{
   su3_dble r;
   if((z.re==1.0)&&(z.im==0.0)){
    r.c11.re=s.c11.re;
    r.c12.re=s.c12.re;
    r.c13.re=s.c13.re;
    r.c21.re=s.c21.re;
    r.c22.re=s.c22.re;
    r.c23.re=s.c23.re;
    r.c31.re=s.c31.re;
    r.c32.re=s.c32.re;
    r.c33.re=s.c33.re;

    r.c11.im=s.c11.im;
    r.c12.im=s.c12.im;
    r.c13.im=s.c13.im;
    r.c21.im=s.c21.im;
    r.c22.im=s.c22.im;
    r.c23.im=s.c23.im;
    r.c31.im=s.c31.im;
    r.c32.im=s.c32.im;
    r.c33.im=s.c33.im;
   }
   else{
    r.c11.re=z.re*s.c11.re-z.im*s.c11.im;
    r.c12.re=z.re*s.c12.re-z.im*s.c12.im;
    r.c13.re=z.re*s.c13.re-z.im*s.c13.im;
    r.c21.re=z.re*s.c21.re-z.im*s.c21.im;
    r.c22.re=z.re*s.c22.re-z.im*s.c22.im;
    r.c23.re=z.re*s.c23.re-z.im*s.c23.im;
    r.c31.re=z.re*s.c31.re-z.im*s.c31.im;
    r.c32.re=z.re*s.c32.re-z.im*s.c32.im;
    r.c33.re=z.re*s.c33.re-z.im*s.c33.im;

    r.c11.im=z.re*s.c11.im+z.im*s.c11.re;
    r.c12.im=z.re*s.c12.im+z.im*s.c12.re;
    r.c13.im=z.re*s.c13.im+z.im*s.c13.re;
    r.c21.im=z.re*s.c21.im+z.im*s.c21.re;
    r.c22.im=z.re*s.c22.im+z.im*s.c22.re;
    r.c23.im=z.re*s.c23.im+z.im*s.c23.re;
    r.c31.im=z.re*s.c31.im+z.im*s.c31.re;
    r.c32.im=z.re*s.c32.im+z.im*s.c32.re;
    r.c33.im=z.re*s.c33.im+z.im*s.c33.re;
   }
   return r;
}
extern full_spinor_dble mul_one_pm_gamma_l(double sign, int mu,full_spinor_dble s)
{
   full_spinor_dble r;
   complex_dble i,m_i,m_1,p_1,csign;
   csign.re = sign;
   csign.im = 0.0;
   i.re=0.0;
   i.im=1.0;

   m_i.re=0.0;
   m_i.im=-1.0;

   m_1.re=-1.0;
   m_1.im=0.0;

   p_1.re=1.0;
   p_1.im=0.0;
   /* The following encoding for the gamma matrices has been used here:
      0 -> gamma_0         ± 1  (time-component)
      1 -> gamma_1         ± 1
      2 -> gamma_2         ± 1
      3 -> gamma_3         ± 1
      4 -> unity           ± 1
      5 -> gamma_5         ± 1
      6 -> gamma_0 gamma_5 ± 1
      7 -> gamma_1 gamma_5 ± 1
      8 -> gamma_2 gamma_5 ± 1
      9 -> gamma_3 gamma_5 ± 1
      10-> gamma_0 gamma_1 ± 1
      11-> gamma_0 gamma_2 ± 1
      12-> gamma_0 gamma_3 ± 1
      13-> gamma_1 gamma_2 ± 1
      14-> gamma_1 gamma_3 ± 1
      15-> gamma_2 gamma_3 ± 1
     */
   if (mu==0)
   {
    _su3_add(r.c11,mul_cplx_su3_dble_out(m_1,s.c31),mul_cplx_su3_dble_out(csign,s.c11));
    _su3_add(r.c12,mul_cplx_su3_dble_out(m_1,s.c32),mul_cplx_su3_dble_out(csign,s.c12));
    _su3_add(r.c13,mul_cplx_su3_dble_out(m_1,s.c33),mul_cplx_su3_dble_out(csign,s.c13));
    _su3_add(r.c14,mul_cplx_su3_dble_out(m_1,s.c34),mul_cplx_su3_dble_out(csign,s.c14));
    _su3_add(r.c21,mul_cplx_su3_dble_out(m_1,s.c41),mul_cplx_su3_dble_out(csign,s.c21));
    _su3_add(r.c22,mul_cplx_su3_dble_out(m_1,s.c42),mul_cplx_su3_dble_out(csign,s.c22));
    _su3_add(r.c23,mul_cplx_su3_dble_out(m_1,s.c43),mul_cplx_su3_dble_out(csign,s.c23));
    _su3_add(r.c24,mul_cplx_su3_dble_out(m_1,s.c44),mul_cplx_su3_dble_out(csign,s.c24));
    _su3_add(r.c31,mul_cplx_su3_dble_out(m_1,s.c11),mul_cplx_su3_dble_out(csign,s.c31));
    _su3_add(r.c32,mul_cplx_su3_dble_out(m_1,s.c12),mul_cplx_su3_dble_out(csign,s.c32));
    _su3_add(r.c33,mul_cplx_su3_dble_out(m_1,s.c13),mul_cplx_su3_dble_out(csign,s.c33));
    _su3_add(r.c34,mul_cplx_su3_dble_out(m_1,s.c14),mul_cplx_su3_dble_out(csign,s.c34));
    _su3_add(r.c41,mul_cplx_su3_dble_out(m_1,s.c21),mul_cplx_su3_dble_out(csign,s.c41));
    _su3_add(r.c42,mul_cplx_su3_dble_out(m_1,s.c22),mul_cplx_su3_dble_out(csign,s.c42));
    _su3_add(r.c43,mul_cplx_su3_dble_out(m_1,s.c23),mul_cplx_su3_dble_out(csign,s.c43));
    _su3_add(r.c44,mul_cplx_su3_dble_out(m_1,s.c24),mul_cplx_su3_dble_out(csign,s.c44));
  }
   else if (mu==1)
   {
    _su3_add(r.c11,mul_cplx_su3_dble_out(m_i,s.c41),mul_cplx_su3_dble_out(csign,s.c11));
    _su3_add(r.c12,mul_cplx_su3_dble_out(m_i,s.c42),mul_cplx_su3_dble_out(csign,s.c12));
    _su3_add(r.c13,mul_cplx_su3_dble_out(m_i,s.c43),mul_cplx_su3_dble_out(csign,s.c13));
    _su3_add(r.c14,mul_cplx_su3_dble_out(m_i,s.c44),mul_cplx_su3_dble_out(csign,s.c14));
    _su3_add(r.c21,mul_cplx_su3_dble_out(m_i,s.c31),mul_cplx_su3_dble_out(csign,s.c21));
    _su3_add(r.c22,mul_cplx_su3_dble_out(m_i,s.c32),mul_cplx_su3_dble_out(csign,s.c22));
    _su3_add(r.c23,mul_cplx_su3_dble_out(m_i,s.c33),mul_cplx_su3_dble_out(csign,s.c23));
    _su3_add(r.c24,mul_cplx_su3_dble_out(m_i,s.c34),mul_cplx_su3_dble_out(csign,s.c24));
    _su3_add(r.c31,mul_cplx_su3_dble_out(i,s.c21),mul_cplx_su3_dble_out(csign,s.c31));
    _su3_add(r.c32,mul_cplx_su3_dble_out(i,s.c22),mul_cplx_su3_dble_out(csign,s.c32));
    _su3_add(r.c33,mul_cplx_su3_dble_out(i,s.c23),mul_cplx_su3_dble_out(csign,s.c33));
    _su3_add(r.c34,mul_cplx_su3_dble_out(i,s.c24),mul_cplx_su3_dble_out(csign,s.c34));
    _su3_add(r.c41,mul_cplx_su3_dble_out(i,s.c11),mul_cplx_su3_dble_out(csign,s.c41));
    _su3_add(r.c42,mul_cplx_su3_dble_out(i,s.c12),mul_cplx_su3_dble_out(csign,s.c42));
    _su3_add(r.c43,mul_cplx_su3_dble_out(i,s.c13),mul_cplx_su3_dble_out(csign,s.c43));
    _su3_add(r.c44,mul_cplx_su3_dble_out(i,s.c14),mul_cplx_su3_dble_out(csign,s.c44));
   }
   else if (mu==2)
   {
    _su3_add(r.c11,mul_cplx_su3_dble_out(m_1,s.c41),mul_cplx_su3_dble_out(csign,s.c11));
    _su3_add(r.c12,mul_cplx_su3_dble_out(m_1,s.c42),mul_cplx_su3_dble_out(csign,s.c12));
    _su3_add(r.c13,mul_cplx_su3_dble_out(m_1,s.c43),mul_cplx_su3_dble_out(csign,s.c13));
    _su3_add(r.c14,mul_cplx_su3_dble_out(m_1,s.c44),mul_cplx_su3_dble_out(csign,s.c14));
    _su3_add(r.c21,mul_cplx_su3_dble_out(p_1,s.c31),mul_cplx_su3_dble_out(csign,s.c21));
    _su3_add(r.c22,mul_cplx_su3_dble_out(p_1,s.c32),mul_cplx_su3_dble_out(csign,s.c22));
    _su3_add(r.c23,mul_cplx_su3_dble_out(p_1,s.c33),mul_cplx_su3_dble_out(csign,s.c23));
    _su3_add(r.c24,mul_cplx_su3_dble_out(p_1,s.c34),mul_cplx_su3_dble_out(csign,s.c24));
    _su3_add(r.c31,mul_cplx_su3_dble_out(p_1,s.c21),mul_cplx_su3_dble_out(csign,s.c31));
    _su3_add(r.c32,mul_cplx_su3_dble_out(p_1,s.c22),mul_cplx_su3_dble_out(csign,s.c32));
    _su3_add(r.c33,mul_cplx_su3_dble_out(p_1,s.c23),mul_cplx_su3_dble_out(csign,s.c33));
    _su3_add(r.c34,mul_cplx_su3_dble_out(p_1,s.c24),mul_cplx_su3_dble_out(csign,s.c34));
    _su3_add(r.c41,mul_cplx_su3_dble_out(m_1,s.c11),mul_cplx_su3_dble_out(csign,s.c41));
    _su3_add(r.c42,mul_cplx_su3_dble_out(m_1,s.c12),mul_cplx_su3_dble_out(csign,s.c42));
    _su3_add(r.c43,mul_cplx_su3_dble_out(m_1,s.c13),mul_cplx_su3_dble_out(csign,s.c43));
    _su3_add(r.c44,mul_cplx_su3_dble_out(m_1,s.c14),mul_cplx_su3_dble_out(csign,s.c44));
   }
   else if (mu==3)
   {
    _su3_add(r.c11,mul_cplx_su3_dble_out(m_i,s.c31),mul_cplx_su3_dble_out(csign,s.c11));
    _su3_add(r.c12,mul_cplx_su3_dble_out(m_i,s.c32),mul_cplx_su3_dble_out(csign,s.c12));
    _su3_add(r.c13,mul_cplx_su3_dble_out(m_i,s.c33),mul_cplx_su3_dble_out(csign,s.c13));
    _su3_add(r.c14,mul_cplx_su3_dble_out(m_i,s.c34),mul_cplx_su3_dble_out(csign,s.c14));
    _su3_add(r.c21,mul_cplx_su3_dble_out(i,s.c41),mul_cplx_su3_dble_out(csign,s.c21));
    _su3_add(r.c22,mul_cplx_su3_dble_out(i,s.c42),mul_cplx_su3_dble_out(csign,s.c22));
    _su3_add(r.c23,mul_cplx_su3_dble_out(i,s.c43),mul_cplx_su3_dble_out(csign,s.c23));
    _su3_add(r.c24,mul_cplx_su3_dble_out(i,s.c44),mul_cplx_su3_dble_out(csign,s.c24));
    _su3_add(r.c31,mul_cplx_su3_dble_out(i,s.c11),mul_cplx_su3_dble_out(csign,s.c31));
    _su3_add(r.c32,mul_cplx_su3_dble_out(i,s.c12),mul_cplx_su3_dble_out(csign,s.c32));
    _su3_add(r.c33,mul_cplx_su3_dble_out(i,s.c13),mul_cplx_su3_dble_out(csign,s.c33));
    _su3_add(r.c34,mul_cplx_su3_dble_out(i,s.c14),mul_cplx_su3_dble_out(csign,s.c34));
    _su3_add(r.c41,mul_cplx_su3_dble_out(m_i,s.c21),mul_cplx_su3_dble_out(csign,s.c41));
    _su3_add(r.c42,mul_cplx_su3_dble_out(m_i,s.c22),mul_cplx_su3_dble_out(csign,s.c42));
    _su3_add(r.c43,mul_cplx_su3_dble_out(m_i,s.c23),mul_cplx_su3_dble_out(csign,s.c43));
    _su3_add(r.c44,mul_cplx_su3_dble_out(m_i,s.c24),mul_cplx_su3_dble_out(csign,s.c44));
   }
   else
   {
      r=fsd0;
      error_root(1,1,"mul_one_pm_gamma_l [algebra.c]",
                 "Gamma matrix %d not defined ",mu);
   }
   return r;
}

