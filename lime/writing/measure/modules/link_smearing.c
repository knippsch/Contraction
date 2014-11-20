/*******************************************************************************
*
* File link_smearing.c
*
* Copyright (C) 2008-9 Bastian Knippschild <knippsch@kph.uni-mainz.de>
*              					   Michele Della Morte <morte@kph.uni-mainz.de>
* 
* All externally accessible functions are explained in README.link_smearing
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/
#define LINK_SMEARING_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "start.h"
#include "global.h"
#include "measure.h"
#include "linalg.h"
#include "random.h"
#include "update.h"
#include "su3.h"
#include "misc.h"
#include "flags.h"


/*
* real number c times SU(3) matrix u
*/
#define _su3_mul(r, c, u)\
(r).c11.re = (u).c11.re * (c); \
(r).c11.im = (u).c11.im * (c); \
(r).c12.re = (u).c12.re * (c); \
(r).c12.im = (u).c12.im * (c); \
(r).c13.re = (u).c13.re * (c); \
(r).c13.im = (u).c13.im * (c); \
(r).c21.re = (u).c21.re * (c); \
(r).c21.im = (u).c21.im * (c); \
(r).c22.re = (u).c22.re * (c); \
(r).c22.im = (u).c22.im * (c); \
(r).c23.re = (u).c23.re * (c); \
(r).c23.im = (u).c23.im * (c); \
(r).c31.re = (u).c31.re * (c); \
(r).c31.im = (u).c31.im * (c); \
(r).c32.re = (u).c32.re * (c); \
(r).c32.im = (u).c32.im * (c); \
(r).c33.re = (u).c33.re * (c); \
(r).c33.im = (u).c33.im * (c)

/*
* complex number c times SU(3) matrix u
*/
#define _su3_mulc(z,c,u) \
  (z).c11.re=(c).re*(u).c11.re-(c).im*(u).c11.im;\
  (z).c11.im=(c).re*(u).c11.im+(c).im*(u).c11.re;\
  (z).c12.re=(c).re*(u).c12.re-(c).im*(u).c12.im;\
  (z).c12.im=(c).re*(u).c12.im+(c).im*(u).c12.re;\
  (z).c13.re=(c).re*(u).c13.re-(c).im*(u).c13.im;\
  (z).c13.im=(c).re*(u).c13.im+(c).im*(u).c13.re;\
  (z).c21.re=(c).re*(u).c21.re-(c).im*(u).c21.im;\
  (z).c21.im=(c).re*(u).c21.im+(c).im*(u).c21.re;\
  (z).c22.re=(c).re*(u).c22.re-(c).im*(u).c22.im;\
  (z).c22.im=(c).re*(u).c22.im+(c).im*(u).c22.re;\
  (z).c23.re=(c).re*(u).c23.re-(c).im*(u).c23.im;\
  (z).c23.im=(c).re*(u).c23.im+(c).im*(u).c23.re;\
  (z).c31.re=(c).re*(u).c31.re-(c).im*(u).c31.im;\
  (z).c31.im=(c).re*(u).c31.im+(c).im*(u).c31.re;\
  (z).c32.re=(c).re*(u).c32.re-(c).im*(u).c32.im;\
  (z).c32.im=(c).re*(u).c32.im+(c).im*(u).c32.re;\
  (z).c33.re=(c).re*(u).c33.re-(c).im*(u).c33.im;\
  (z).c33.im=(c).re*(u).c33.im+(c).im*(u).c33.re


/*
* i times SU(3) matrix u
*/
#define _su3_i_mul(r, u)\
(r).c11.re = -(u).c11.im; \
(r).c11.im = (u).c11.re; \
(r).c12.re = -(u).c12.im; \
(r).c12.im = (u).c12.re; \
(r).c13.re = -(u).c13.im; \
(r).c13.im = (u).c13.re; \
(r).c21.re = -(u).c21.im; \
(r).c21.im = (u).c21.re; \
(r).c22.re = -(u).c22.im; \
(r).c22.im = (u).c22.re; \
(r).c23.re = -(u).c23.im; \
(r).c23.im = (u).c23.re; \
(r).c31.re = -(u).c31.im; \
(r).c31.im = (u).c31.re; \
(r).c32.re = -(u).c32.im; \
(r).c32.im = (u).c32.re; \
(r).c33.re = -(u).c33.im; \
(r).c33.im = (u).c33.re


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


/*
* trace of a SU(3) matrix times the unity matrix
*/
#define _su3_trace(r, u)\
(r).c11.re = (u).c11.re+(u).c22.re+(u).c33.re; \
(r).c11.im = (u).c11.im+(u).c22.im+(u).c33.im; \
(r).c12.re = 0; \
(r).c12.im = 0; \
(r).c13.re = 0; \
(r).c13.im = 0; \
(r).c21.re = 0; \
(r).c21.im = 0; \
(r).c22.re = (u).c11.re+(u).c22.re+(u).c33.re; \
(r).c22.im = (u).c11.im+(u).c22.im+(u).c33.im; \
(r).c23.re = 0; \
(r).c23.im = 0; \
(r).c31.re = 0; \
(r).c31.im = 0; \
(r).c32.re = 0; \
(r).c32.im = 0; \
(r).c33.re = (u).c11.re+(u).c22.re+(u).c33.re; \
(r).c33.im = (u).c11.im+(u).c22.im+(u).c33.im


/*
* u=u+r*v (r real; u,v su3)
*/
#define _su3_mulr_add_assign(u,v,r) \
   (u).c11.re+=(r)*(v).c11.re; \
   (u).c11.im+=(r)*(v).c11.im; \
   (u).c12.re+=(r)*(v).c12.re; \
   (u).c12.im+=(r)*(v).c12.im; \
   (u).c13.re+=(r)*(v).c13.re; \
   (u).c13.im+=(r)*(v).c13.im; \
   (u).c21.re+=(r)*(v).c21.re; \
   (u).c21.im+=(r)*(v).c21.im; \
   (u).c22.re+=(r)*(v).c22.re; \
   (u).c22.im+=(r)*(v).c22.im; \
   (u).c23.re+=(r)*(v).c23.re; \
   (u).c23.im+=(r)*(v).c23.im; \
   (u).c31.re+=(r)*(v).c31.re; \
   (u).c31.im+=(r)*(v).c31.im; \
   (u).c32.re+=(r)*(v).c32.re; \
   (u).c32.im+=(r)*(v).c32.im; \
   (u).c33.re+=(r)*(v).c33.re; \
   (u).c33.im+=(r)*(v).c33.im


#define _su3_addr(z,a) \
  (z).c11.re+=a;\
  (z).c22.re+=a;\
  (z).c33.re+=a

#define SMEARLINKS  

typedef struct
{
   int nu0;
   int saddr,raddr;
   int iu0;
} comdat_t_reverse;


static su3_dble *C,/* *Omega, *Q,*/ *udz1=NULL, *udz2=NULL;
static const su3_dble uddd0={{0.0}};

static int nfc[4],ofs[4];
static int *idx_u0=NULL;

static su3_dble *sdbuf_u0=NULL,*rdbuf_u0;
static comdat_t_reverse comdat[4];


static double det_im_dble(su3_dble *u)
{
   double detuim;
   complex_dble det1,det2,det3;

   det1.re=
         ((*u).c22.re*(*u).c33.re-(*u).c22.im*(*u).c33.im)-
         ((*u).c23.re*(*u).c32.re-(*u).c23.im*(*u).c32.im);
   det1.im=
         ((*u).c22.re*(*u).c33.im+(*u).c22.im*(*u).c33.re)-
         ((*u).c23.re*(*u).c32.im+(*u).c23.im*(*u).c32.re);
   det2.re=
         ((*u).c21.re*(*u).c33.re-(*u).c21.im*(*u).c33.im)-
         ((*u).c23.re*(*u).c31.re-(*u).c23.im*(*u).c31.im);
   det2.im=
         ((*u).c21.re*(*u).c33.im+(*u).c21.im*(*u).c33.re)-
         ((*u).c23.re*(*u).c31.im+(*u).c23.im*(*u).c31.re);
   det3.re=
         ((*u).c21.re*(*u).c32.re-(*u).c21.im*(*u).c32.im)-
         ((*u).c22.re*(*u).c31.re-(*u).c22.im*(*u).c31.im);
   det3.im=
         ((*u).c21.re*(*u).c32.im+(*u).c21.im*(*u).c32.re)-
         ((*u).c22.re*(*u).c31.im+(*u).c22.im*(*u).c31.re);

   detuim=
         ((*u).c11.re*det1.im+(*u).c11.im*det1.re)-
         ((*u).c12.re*det2.im+(*u).c12.im*det2.re)+
         ((*u).c13.re*det3.im+(*u).c13.im*det3.re);

   return detuim;
}


/************************************/
/*	  a covariant way to project	*/
/*			back to SU(3)			*/
/************************************/
extern void approx_project_to_su3_dble(su3_dble *u, int iter){

   int i;
   double d;
   complex_dble z;
   su3_dble x;
   
   d= (*u).c11.re*(*u).c11.re+(*u).c11.im*(*u).c11.im  \
         +(*u).c12.re*(*u).c12.re+(*u).c12.im*(*u).c12.im \
         +(*u).c13.re*(*u).c13.re+(*u).c13.im*(*u).c13.im \
         +(*u).c21.re*(*u).c21.re+(*u).c21.im*(*u).c21.im \
         +(*u).c22.re*(*u).c22.re+(*u).c22.im*(*u).c22.im \
         +(*u).c23.re*(*u).c23.re+(*u).c23.im*(*u).c23.im \
         +(*u).c31.re*(*u).c31.re+(*u).c31.im*(*u).c31.im \
         +(*u).c32.re*(*u).c32.re+(*u).c32.im*(*u).c32.im \
         +(*u).c33.re*(*u).c33.re+(*u).c33.im*(*u).c33.im;
   d=1.0/sqrt(d/3.0);
   
   _su3_mul(*u,d,*u);
   
   z.re=1.0;
   for (i=0;i<iter;i++){
      su3dagxsu3(u,u,&x);
      _su3_mul(x,-0.5,x);
      _su3_addr(x,1.5);
      su3xsu3(u,&x,&x);

      z.im=det_im_dble(&x);
      z.im/=-3.0;
      _su3_mulc(*u,z,x);
   }
}


/************************************/
/*  allocates the memoryspcace for	*/
/*		  the lookuptable for		*/
/*		  copy_bnd_ud_reverse		*/
/************************************/
static void alloc_idx_reverse(void)
{
   int n,mu,iu0;
   comdat_t_reverse *c;

   nfc[0]=FACE0;
   nfc[1]=FACE1;
   nfc[2]=FACE2;
   nfc[3]=FACE3;

   ofs[0]=FACE0/2;
   ofs[1]=ofs[0]+(FACE0+FACE1)/2;
   ofs[2]=ofs[1]+(FACE1+FACE2)/2;
   ofs[3]=ofs[2]+(FACE2+FACE3)/2;

   n=BNDRY/4;
   idx_u0=amalloc(n*sizeof(int),3);
   error(idx_u0==NULL,1,"alloc_idx_diverse [ucom.c]","Unable to allocate index array");

   iu0=0;

   for (mu=0;mu<4;mu++)
   {
      c=comdat+mu;

      (*c).nu0=nfc[mu]/2;
      (*c).saddr=npr[2*mu];
      (*c).raddr=npr[2*mu+1];
      (*c).iu0=iu0;

      iu0+=(*c).nu0;
   }
}


/****************************************/
/* lookuptable for copy_bnd_ud reverse	*/
/****************************************/
static void set_idx_reverse(void)
{
   int mu;
   int ioe,ioo,ix,iyo;
   int nu0,*u0;

   u0=idx_u0;

   for (mu=0;mu<4;mu++)
   {
      nu0=nfc[mu]/2;
      ioe=ofs[mu];
      ioo=ioe+(BNDRY/2);

      for (ix=0;ix<nu0;ix++)
      {
         iyo=map[ioo+ix];

         *u0=8*(iyo-(VOLUME/2))+2*mu+1;
         u0+=1;
	  }
   }
}


/************************************/
/* allocates the needed memoryspace	*/
/*		for copy_bnd_ud reverse		*/
/************************************/
static void alloc_udbufs_reverse(su3_dble *u1[VOLUME][4])
{
   int n;

   error(pud[VOLUME/2][0]==NULL,1,"alloc_udbufs_diverse [ucom.c]",
         "Double-precision gauge field is not allocated");

   if (idx_u0==NULL)
   {
      alloc_idx_reverse();
      set_idx_reverse();
   }

   n=BNDRY/4;
   sdbuf_u0=amalloc(n*sizeof(su3_dble),ALIGN);
   error(sdbuf_u0==NULL,1,"alloc_udbufs_diverse [ucom.c]","Unable to allocate buffers");

   rdbuf_u0=u1[VOLUME/2][0]+4*VOLUME;
}


/************************************/
/*  sorts the needed linkvariables	*/
/*		for copy_bnd_ud reverse		*/
/************************************/
static void pack_ud0_reverse(su3_dble *u1[VOLUME][4])
{
   int *iu,*ium;
   su3_dble *u,*ub;

   u=sdbuf_u0;
   ub=u1[VOLUME/2][0];
   iu=idx_u0;
   ium=idx_u0+(BNDRY/4);

   for (;iu<ium;iu++)
   {
	  *(ub+(*iu)) = *u;
      u+=1;
   }
}


/************************************/
/*  sends the needed linkvariables	*/
/*		for copy_bnd_ud reverse		*/
/************************************/
static void send_ud0_reverse(void)
{
   int tag,nbf,saddr,raddr,np;
   double *sbuf,*rbuf;
   comdat_t_reverse *c,*cm;
   MPI_Status stat;

   np=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;
   c=comdat;
   cm=c+4;

   for (;c<cm;c++)
   {
      nbf=18*(*c).nu0;

      if (nbf>0)
      {
         tag=mpi_tag();
         saddr=(*c).saddr;
         raddr=(*c).raddr;
         sbuf=(double*)(sdbuf_u0+(*c).iu0);
         rbuf=(double*)(rdbuf_u0+(*c).iu0);

         if (np==0)
         {
            MPI_Send(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD);
            MPI_Recv(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD,&stat);
         }
         else
         {
            MPI_Recv(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD,&stat);
            MPI_Send(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD);
         }
      }
   }
}




/************************************/
/*   frees the variables which are	*/
/*	 needed for copy_bnd_ud_reverse */
/************************************/
extern void free_ucom_bufs_reverse(void)
{

   if ((sdbuf_u0!=NULL))
   {
      afree(sdbuf_u0);
      sdbuf_u0=NULL;
   }
}

/************************************/
/* copies the linkvariables from pud*/
/* in the correct boundary positions*/
/************************************/
extern void copy_bnd_ud_reverse(su3_dble *u1[VOLUME][4])
{
   if (NPROC>1)
   {
      if (sdbuf_u0==NULL)
         alloc_udbufs_reverse(u1);
		 
      send_ud0_reverse();
      pack_ud0_reverse(u1);

   }
free_ucom_bufs_reverse();
}

/************************************/
/*   Allocates the dummy variable	*/
/*			  pud_sm1				*/
/************************************/
extern void alloc_pud_sm1(void)
{
   int ix,iy,mu;
   int io[4],nu[4];
   su3_dble unity,*p;

   error(iup[0][0]==0,1,"alloc_ud [start.c]",
         "Geometry arrays are not initialized");

   error_root(sizeof(su3_dble)!=(18*sizeof(double)),1,"alloc_pud_APE [start.c]",
         "The su3_dble structures are not properly packed in the link smearing");
   
   if (udz1==NULL)
      udz1=amalloc((4*VOLUME+BNDRY/4)*sizeof(su3_dble),ALIGN);
   error(udz1==NULL,1,"alloc_pud_APE [start.c]",
         "Could not allocate memory space for the gauge fieldin the link smearing");

   unity=uddd0;
   unity.c11.re=1.0;
   unity.c22.re=1.0;
   unity.c33.re=1.0;
   p=udz1;

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      for (mu=0;mu<4;mu++)
      {
         *p=unity;
         pud_sm1[ix][mu]=p;
         p+=1;

         iy=idn[ix][mu];
         *p=unity;
         if (iy<VOLUME)
            pud_sm1[iy][mu]=p;
         p+=1;
      }
   }

   io[0]=      (BNDRY+FACE0)/2;
   io[1]=io[0]+(FACE0+FACE1)/2;
   io[2]=io[1]+(FACE1+FACE2)/2;
   io[3]=io[2]+(FACE2+FACE3)/2;

   nu[0]=FACE0/2;
   nu[1]=FACE1/2;
   nu[2]=FACE2/2;
   nu[3]=FACE3/2;

   for (mu=0;mu<4;mu++)
   {
      for (iy=0;iy<nu[mu];iy++)
      {
         ix=map[io[mu]+iy];
         ix=idn[ix][mu];
         ix=map[ix-VOLUME];

         *p=unity;
         pud_sm1[ix][mu]=p;
         p+=1;
      }
   }
}


/************************************/
/* Frees the dummy variable pud_sm1	*/
/************************************/
extern void free_pud_sm1(void)
{
   int ix,mu;

   if (udz1!=NULL)
   {
      afree(udz1);
      udz1=NULL;
   }

   for (ix=0;ix<VOLUME;ix++)
   {
      for (mu=0;mu<4;mu++)
         pud_sm1[ix][mu]=NULL;
   }

}


/************************************/
/*   Allocates the dummy variable	*/
/*			  pud_sm2				*/
/************************************/
extern void alloc_pud_sm2(void)
{
   int ix,iy,mu;
   int io[4],nu[4];
   su3_dble unity,*p;

   error(iup[0][0]==0,1,"alloc_ud [start.c]",
         "Geometry arrays are not initialized");

   error_root(sizeof(su3_dble)!=(18*sizeof(double)),1,"alloc_pud_APE [start.c]",
         "The su3_dble structures are not properly packed in the link smearing");
   
   if (udz2==NULL)
      udz2=amalloc((4*VOLUME+BNDRY/4)*sizeof(su3_dble),ALIGN);
   error(udz2==NULL,1,"alloc_pud_HYP [start.c]",
         "Could not allocate memory space for the gauge fieldin the link smearing");

   unity=uddd0;
   unity.c11.re=1.0;
   unity.c22.re=1.0;
   unity.c33.re=1.0;
   p=udz2;

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      for (mu=0;mu<4;mu++)
      {
         *p=unity;
         pud_sm2[ix][mu]=p;
         p+=1;

         iy=idn[ix][mu];
         *p=unity;
         if (iy<VOLUME)
            pud_sm2[iy][mu]=p;
         p+=1;
      }
   }

   io[0]=      (BNDRY+FACE0)/2;
   io[1]=io[0]+(FACE0+FACE1)/2;
   io[2]=io[1]+(FACE1+FACE2)/2;
   io[3]=io[2]+(FACE2+FACE3)/2;

   nu[0]=FACE0/2;
   nu[1]=FACE1/2;
   nu[2]=FACE2/2;
   nu[3]=FACE3/2;

   for (mu=0;mu<4;mu++)
   {
      for (iy=0;iy<nu[mu];iy++)
      {
         ix=map[io[mu]+iy];
         ix=idn[ix][mu];
         ix=map[ix-VOLUME];

         *p=unity;
         pud_sm2[ix][mu]=p;
         p+=1;
      }
   }
}


/************************************/
/* Frees the dummy variable pud_sm2	*/
/************************************/
extern void free_pud_sm2(void)
{
   int ix,mu;

   if (udz2!=NULL)
   {
      afree(udz2);
      udz2=NULL;
   }

   for (ix=0;ix<VOLUME;ix++)
   {
      for (mu=0;mu<4;mu++)
         pud_sm2[ix][mu]=NULL;
   }

}

/************************************/
/*	Gives the ips for sending and	*/
/*	 recieving the linkvariables	*/
/*	 for a L-shape communication	*/
/************************************/
static void L_shape_comm(int *ip1, int *ip2, int mu2, int nu2){

 int L_up[4], L_dn[4], i;
	
 for (i=0; i<4; i++){
	L_up[i] = cpr[i];
	L_dn[i] = cpr[i];
 }
 
 L_up[nu2] = cpr[nu2]-1;
 L_up[mu2] = cpr[mu2]+1;
 L_dn[nu2] = cpr[nu2]+1;
 L_dn[mu2] = cpr[mu2]-1;
 
 *ip1 = ipr_global(L_up);
 *ip2 = ipr_global(L_dn);
 
}


/************************************/
/*communication of the linkvariables*/
/************************************/
static void u_comm(char direction, int iy, int mu1, int nu1, su3_dble *u, su3_dble *u1[VOLUME][4]){

 int ip1=0, ip2=0, iw=0;
 double *zw;
 
 MPI_Status status1;

 zw = amalloc(18*sizeof(double), 4);

/* up direction */
 if(direction == 'u'){ 
	ip1 = npr[2*mu1+1];
	ip2 = npr[2*mu1];
	iw = map[iy-VOLUME];
 }
/* down direction */
 else if(direction == 'd'){
	ip1 = npr[2*mu1];
	ip2 = npr[2*mu1+1];
	iw = map[iy-VOLUME];
 }
/* L-shape */
 else if(direction == 'L'){
	L_shape_comm(&ip1, &ip2, mu1, nu1);
	iw = map[iy-VOLUME];		
 }
 else message("unknown character for the direction in the link-communication \n");
 
 *u = *u1[iw][nu1];
 
 MPI_Sendrecv((double*)(u), 18, MPI_DOUBLE, ip2, 37, zw, 18, MPI_DOUBLE, ip1, 37, MPI_COMM_WORLD, &status1);
 
 *u = *((su3_dble*)zw);

 afree(zw);
 
}


/************************************/
/*   Calculates the V_bar staples	*/
/*		 for the APE-Smearing		*/
/************************************/
static void C_mu(int ix, int mu, su3_dble *u1[VOLUME][4], double alpha1){

 int nu=0, i, iy=0, iz=0;
 su3_dble *C1, *C2, *C3, *u1_mu, *u2_mu, *u1_nu, *u2_nu, *u2_nu_dagger, *u3_nu, *u3_nu_dagger, *u4_nu;
 double *p;

 C1 = amalloc(sizeof(su3_dble), 4);
 C2 = amalloc(sizeof(su3_dble), 4);
 C3 = amalloc(sizeof(su3_dble), 4);
 u1_mu = amalloc(sizeof(su3_dble), 4);
 u2_mu = amalloc(sizeof(su3_dble), 4);
 u1_nu = amalloc(sizeof(su3_dble), 4);
 u2_nu = amalloc(sizeof(su3_dble), 4);
 u3_nu = amalloc(sizeof(su3_dble), 4);
 u4_nu = amalloc(sizeof(su3_dble), 4);
 u2_nu_dagger = amalloc(sizeof(su3_dble), 4);
 u3_nu_dagger = amalloc(sizeof(su3_dble), 4);

 p = (double*)C;
 for(i=0; i<18; i++){
	p[i] = 0.0;
 }

 for(nu=0; nu<4; nu++){
	
	if(nu == mu) continue;
					
	iy = iup[ix][nu];
	if((iy>=VOLUME) && (iy<(VOLUME+BNDRY))) u_comm('u', iy, nu, mu, u1_mu, u1);
	else *u1_mu = *u1[iy][mu];
		
	*u1_nu = *u1[ix][nu];
			
	iy = iup[ix][mu];
	if((iy>=VOLUME) && (iy<(VOLUME+BNDRY))) u_comm('u', iy, mu, nu, u2_nu, u1);
	else *u2_nu = *u1[iy][nu];
		
	_su3_times_su3((*C1), (*u1_nu), (*u1_mu));
	_su3_dagger((*u2_nu_dagger), (*u2_nu));
	_su3_times_su3((*C2), (*C1), (*u2_nu_dagger));

	
	iy = idn[ix][nu];
	if((iy>=VOLUME) && (iy<(VOLUME+BNDRY))){
		u_comm('d', iy, nu, mu, u2_mu, u1);
		u_comm('d', iy, nu, nu, u3_nu, u1);
	}
	else{
		*u2_mu = *u1[iy][mu];
		*u3_nu = *u1[iy][nu];
	}
	
	iy = idn[ix][nu];
	iz = iup[ix][mu];
	if((iy>=VOLUME) && (iy<(VOLUME+BNDRY)) && (iz<VOLUME)){
		iy = idn[iz][nu];
		u_comm('d', iy, nu, nu, u4_nu, u1);
	}
	else if((iz>=VOLUME) && (iz<(VOLUME+BNDRY)) && (iy<VOLUME)){
		iy = iup[iy][mu];
		u_comm('u', iy, mu, nu, u4_nu, u1);
	}
	else if(((iy>=VOLUME) && (iy<(VOLUME+BNDRY))) && ((iz>=VOLUME) && (iz<(VOLUME+BNDRY)))){
		iy = iup[map[iy-VOLUME]][mu];
		u_comm('L', iy, mu, nu, u4_nu, u1);
	}
	else{
		iy = iup[iy][mu];
		*u4_nu = *u1[iy][nu];
	} 
 
 	_su3_dagger((*u3_nu_dagger), (*u3_nu));
	_su3_times_su3((*C1), (*u3_nu_dagger), (*u2_mu));
	_su3_times_su3((*C3), (*C1), (*u4_nu));
	
	_su3_add((*C1), (*C2), (*C3));		
	_su3_add((*C), (*C), (*C1));
	
 }
 	
 *C1 = *u1[ix][mu];
 _su3_mul((*C1), 1.0-alpha1, (*C1));
 
 _su3_mul((*C), alpha1/6.0, (*C));
 _su3_add((*C), (*C), (*C1));

 
 afree(C1);
 afree(C2);
 afree(C3);
 afree(u1_mu);
 afree(u2_mu);
 afree(u1_nu);
 afree(u2_nu);
 afree(u3_nu);
 afree(u4_nu);
 afree(u2_nu_dagger);
 afree(u3_nu_dagger);

}


/************************************/
/*	    Calculates the staples		*/
/*		 for the 1st step in		*/
/*		  the HYP-Smearing			*/
/************************************/
static void C_HYP_1(int ix, int mu, int eta, su3_dble *u1[VOLUME][4], double alpha3){

 int i, iy=0, iz=0;
 su3_dble *C1, *C2, *C3, *u1_mu, *u2_mu, *u1_nu, *u2_nu, *u2_nu_dagger, *u3_nu, *u3_nu_dagger, *u4_nu;
 double *p;

 C1 = amalloc(sizeof(su3_dble), 4);
 C2 = amalloc(sizeof(su3_dble), 4);
 C3 = amalloc(sizeof(su3_dble), 4);
 u1_mu = amalloc(sizeof(su3_dble), 4);
 u2_mu = amalloc(sizeof(su3_dble), 4);
 u1_nu = amalloc(sizeof(su3_dble), 4);
 u2_nu = amalloc(sizeof(su3_dble), 4);
 u3_nu = amalloc(sizeof(su3_dble), 4);
 u4_nu = amalloc(sizeof(su3_dble), 4);
 u2_nu_dagger = amalloc(sizeof(su3_dble), 4);
 u3_nu_dagger = amalloc(sizeof(su3_dble), 4);

 p = (double*)C;
 for(i=0; i<18; i++){
	p[i] = 0.0;
 }
					
	iy = iup[ix][eta];
	if((iy>=VOLUME) && (iy<(VOLUME+BNDRY))) u_comm('u', iy, eta, mu, u1_mu, u1);
	else *u1_mu = *u1[iy][mu];
		
	*u1_nu = *u1[ix][eta];
			
	iy = iup[ix][mu];
	if((iy>=VOLUME) && (iy<(VOLUME+BNDRY))) u_comm('u', iy, mu, eta, u2_nu, u1);
	else *u2_nu = *u1[iy][eta];
		
	_su3_times_su3((*C1), (*u1_nu), (*u1_mu));
	_su3_dagger((*u2_nu_dagger), (*u2_nu));
	_su3_times_su3((*C2), (*C1), (*u2_nu_dagger));

	
	iy = idn[ix][eta];
	if((iy>=VOLUME) && (iy<(VOLUME+BNDRY))){
		u_comm('d', iy, eta, mu, u2_mu, u1);
		u_comm('d', iy, eta, eta, u3_nu, u1);
	}
	else{
		*u2_mu = *u1[iy][mu];
		*u3_nu = *u1[iy][eta];
	}
	
	iy = idn[ix][eta];
	iz = iup[ix][mu];
	if((iy>=VOLUME) && (iy<(VOLUME+BNDRY)) && (iz<VOLUME)){
		iy = idn[iz][eta];
		u_comm('d', iy, eta, eta, u4_nu, u1);
	}
	else if((iz>=VOLUME) && (iz<(VOLUME+BNDRY)) && (iy<VOLUME)){
		iy = iup[iy][mu];
		u_comm('u', iy, mu, eta, u4_nu, u1);
	}
	else if(((iy>=VOLUME) && (iy<(VOLUME+BNDRY))) && ((iz>=VOLUME) && (iz<(VOLUME+BNDRY)))){
		iy = iup[map[iy-VOLUME]][mu];
		u_comm('L', iy, mu, eta, u4_nu, u1);
	}
	else{
		iy = iup[iy][mu];
		*u4_nu = *u1[iy][eta];
	}
	
	_su3_dagger((*u3_nu_dagger), (*u3_nu));
	_su3_times_su3((*C1), (*u3_nu_dagger), (*u2_mu));
	_su3_times_su3((*C3), (*C1), (*u4_nu));

	
	_su3_add((*C1), (*C2), (*C3));	
	_su3_add((*C), (*C), (*C1));
	
 *C1 = *u1[ix][mu];
 _su3_mul((*C1), 1.0-alpha3, (*C1));
 
 _su3_mul((*C), alpha3/2.0, (*C));
 _su3_add((*C), (*C), (*C1));
 
 
 afree(C1);
 afree(C2);
 afree(C3);
 afree(u1_mu);
 afree(u2_mu);
 afree(u1_nu);
 afree(u2_nu);
 afree(u3_nu);
 afree(u4_nu);
 afree(u2_nu_dagger);
 afree(u3_nu_dagger);

}


/************************************/
/*	    lookuptable for C_HYP_2		*/
/************************************/
static int lookup_HYP_2(int mu, int nu, int roh, int *ord){
 
 int j=0, i=0;
 
 for(i=0; i<12; i++){
	if((*(ord+4*i) == mu) && (*(ord+4*i+1) == nu) && (*(ord+4*i+2) == roh)) j = *(ord+4*i+3);
	else if((*(ord+4*i) == mu) && (*(ord+4*i+1) == roh) && (*(ord+4*i+2) == nu)) j = *(ord+4*i+3);
 }
 
 return j;
}


/************************************/
/*communication of the linkvariables*/
/*		for C_HYP_2	and C_HYP_3		*/
/************************************/
static void u_comm_HYP(char direction, int iy, int mu1, int nu1, int j, su3_dble *u, su3_dble *V){

 int ip1=0, ip2=0, iw=0;
 double *zw;
 
 MPI_Status status1;

 zw = amalloc(18*sizeof(double), 4);

/* up direction */
 if(direction == 'u'){ 
	ip1 = npr[2*mu1+1];
	ip2 = npr[2*mu1];
	iw = map[iy-VOLUME];
 }
/* down direction */
 else if(direction == 'd'){
	ip1 = npr[2*mu1];
	ip2 = npr[2*mu1+1];
	iw = map[iy-VOLUME];
 }
/* L-shape */
 else if(direction == 'L'){
	L_shape_comm(&ip1, &ip2, mu1, nu1);
	iw = map[iy-VOLUME];		
 }
 else message("unknown character for the direction in the link-communication \n");
 
 *u = *(V+12*iw+j);
 
 MPI_Sendrecv((double*)(u), 18, MPI_DOUBLE, ip2, 37, zw, 18, MPI_DOUBLE, ip1, 37, MPI_COMM_WORLD, &status1);
 
 *u = *((su3_dble*)zw);

 afree(zw);
 
}


/************************************/
/*	    Calculates the staples		*/
/*		 for the 2nd step in		*/
/*		  the HYP-Smearing			*/
/************************************/
static void C_HYP_2(int ix, int mu, int nu, int *ord, su3_dble *V_bar, su3_dble *u1[VOLUME][4], double alpha2){

 int i, iy=0, iz=0, roh, j;
 su3_dble *C1, *C2, *C3, *u1_mu, *u2_mu, *u1_nu, *u2_nu, *u2_nu_dagger, *u3_nu, *u3_nu_dagger, *u4_nu;
 double *p;

 C1 = amalloc(sizeof(su3_dble), 4);
 C2 = amalloc(sizeof(su3_dble), 4);
 C3 = amalloc(sizeof(su3_dble), 4);
 u1_mu = amalloc(sizeof(su3_dble), 4);
 u2_mu = amalloc(sizeof(su3_dble), 4);
 u1_nu = amalloc(sizeof(su3_dble), 4);
 u2_nu = amalloc(sizeof(su3_dble), 4);
 u3_nu = amalloc(sizeof(su3_dble), 4);
 u4_nu = amalloc(sizeof(su3_dble), 4);
 u2_nu_dagger = amalloc(sizeof(su3_dble), 4);
 u3_nu_dagger = amalloc(sizeof(su3_dble), 4);

 p = (double*)C;
 for(i=0; i<18; i++){
	p[i] = 0.0;
 }
 
 for(roh=0; roh<4; roh++){
	if((roh == nu) || (roh == mu)) continue;
					
	iy = iup[ix][roh];
	j = lookup_HYP_2(mu, roh, nu, ord);
	if((iy>=VOLUME) && (iy<(VOLUME+BNDRY))) u_comm_HYP('u', iy, roh, mu, j, u1_mu, V_bar);
	else *u1_mu = *(V_bar+12*iy+j);
		
	j = lookup_HYP_2(roh, nu, mu, ord);
	*u1_nu = *(V_bar+12*ix+j);
			
	iy = iup[ix][mu];
	j = lookup_HYP_2(roh, nu, mu, ord);
	if((iy>=VOLUME) && (iy<(VOLUME+BNDRY))) u_comm_HYP('u', iy, mu, roh, j, u2_nu, V_bar);
	else *u2_nu = *(V_bar+12*iy+j);
		
	_su3_times_su3((*C1), (*u1_nu), (*u1_mu));
	_su3_dagger((*u2_nu_dagger), (*u2_nu));
	_su3_times_su3((*C2), (*C1), (*u2_nu_dagger));

	
	iy = idn[ix][roh];
	if((iy>=VOLUME) && (iy<(VOLUME+BNDRY))){
		j = lookup_HYP_2(mu, nu, roh, ord);
		u_comm_HYP('d', iy, roh, mu, j, u2_mu, V_bar);
		j = lookup_HYP_2(roh, nu, mu, ord);
		u_comm_HYP('d', iy, roh, roh, j, u3_nu, V_bar);
	}
	else{
		j = lookup_HYP_2(mu, nu, roh, ord);
		*u2_mu = *(V_bar+12*iy+j);
		j = lookup_HYP_2(roh, nu, mu, ord);
		*u3_nu = *(V_bar+12*iy+j);
	}
	
	iy = idn[ix][roh];
	iz = iup[ix][mu];
	j = lookup_HYP_2(roh, nu, mu, ord);
	if((iy>=VOLUME) && (iy<(VOLUME+BNDRY)) && (iz<VOLUME)){
		iy = idn[iz][roh];
		u_comm_HYP('d', iy, roh, roh, j, u4_nu, V_bar);
	}
	else if((iz>=VOLUME) && (iz<(VOLUME+BNDRY)) && (iy<VOLUME)){
		iy = iup[iy][mu];
		u_comm_HYP('u', iy, mu, roh, j, u4_nu, V_bar);
	}
	else if(((iy>=VOLUME) && (iy<(VOLUME+BNDRY))) && ((iz>=VOLUME) && (iz<(VOLUME+BNDRY)))){
		iy = iup[map[iy-VOLUME]][mu];
		u_comm_HYP('L', iy, mu, roh, j, u4_nu, V_bar);
	}
	else{
		iy = iup[iy][mu];
		*u4_nu = *(V_bar+12*iy+j);
	}
	
	_su3_dagger((*u3_nu_dagger), (*u3_nu));
	_su3_times_su3((*C1), (*u3_nu_dagger), (*u2_mu));
	_su3_times_su3((*C3), (*C1), (*u4_nu));

	
	_su3_add((*C1), (*C2), (*C3));		
	_su3_add((*C), (*C), (*C1));
 }
 	
 *C1 = *u1[ix][mu];
 _su3_mul((*C1), 1.0-alpha2, (*C1));
 
 _su3_mul((*C), alpha2/4.0, (*C));
 _su3_add((*C), (*C), (*C1));
 
 
 afree(C1);
 afree(C2);
 afree(C3);
 afree(u1_mu);
 afree(u2_mu);
 afree(u1_nu);
 afree(u2_nu);
 afree(u3_nu);
 afree(u4_nu);
 afree(u2_nu_dagger);
 afree(u3_nu_dagger);

}


/************************************/
/*	    lookuptable for C_HYP_3		*/
/************************************/
static int lookup_HYP_3(int mu, int nu, int *ord){
 
 int j=0, i=0;
 
 for(i=0; i<12; i++){
	if((*(ord+3*i) == mu) && (*(ord+3*i+1) == nu)) 
		j = *(ord+3*i+2);
 }
 
 return j;
}


/************************************/
/*	    Calculates the staples		*/
/*		 for the 3rd step in		*/
/*		  the HYP-Smearing			*/
/************************************/
static void C_HYP_3(int ix, int mu, int *ord, su3_dble *V_tilde, su3_dble *u1[VOLUME][4], double alpha1){

 int i, iy=0, iz=0, nu, j;
 su3_dble *C1, *C2, *C3, *u1_mu, *u2_mu, *u1_nu, *u2_nu, *u2_nu_dagger, *u3_nu, *u3_nu_dagger, *u4_nu;
 double *p;

 C1 = amalloc(sizeof(su3_dble), 4);
 C2 = amalloc(sizeof(su3_dble), 4);
 C3 = amalloc(sizeof(su3_dble), 4);
 u1_mu = amalloc(sizeof(su3_dble), 4);
 u2_mu = amalloc(sizeof(su3_dble), 4);
 u1_nu = amalloc(sizeof(su3_dble), 4);
 u2_nu = amalloc(sizeof(su3_dble), 4);
 u3_nu = amalloc(sizeof(su3_dble), 4);
 u4_nu = amalloc(sizeof(su3_dble), 4);
 u2_nu_dagger = amalloc(sizeof(su3_dble), 4);
 u3_nu_dagger = amalloc(sizeof(su3_dble), 4);

 p = (double*)C;
 for(i=0; i<18; i++){
	p[i] = 0.0;
 }
 
 for(nu=0; nu<4; nu++){
	if(nu == mu) continue;
					
	iy = iup[ix][nu];
	j = lookup_HYP_3(mu, nu, ord);
	if((iy>=VOLUME) && (iy<(VOLUME+BNDRY))) u_comm_HYP('u', iy, nu, mu, j, u1_mu, V_tilde);
	else *u1_mu = *(V_tilde+12*iy+j);
		
	j = lookup_HYP_3(nu, mu, ord);
	*u1_nu = *(V_tilde+12*ix+j);
			
	iy = iup[ix][mu];
	j = lookup_HYP_3(nu, mu, ord);
	if((iy>=VOLUME) && (iy<(VOLUME+BNDRY))) u_comm_HYP('u', iy, mu, nu, j, u2_nu, V_tilde);
	else *u2_nu = *(V_tilde+12*iy+j);
		
	_su3_times_su3((*C1), (*u1_nu), (*u1_mu));
	_su3_dagger((*u2_nu_dagger), (*u2_nu));
	_su3_times_su3((*C2), (*C1), (*u2_nu_dagger));

	
	iy = idn[ix][nu];
	if((iy>=VOLUME) && (iy<(VOLUME+BNDRY))){
		j = lookup_HYP_3(mu, nu, ord);
		u_comm_HYP('d', iy, nu, mu, j, u2_mu, V_tilde);
		j = lookup_HYP_3(nu, mu, ord);
		u_comm_HYP('d', iy, nu, nu, j, u3_nu, V_tilde);
	}
	else{
		j = lookup_HYP_3(mu, nu, ord);
		*u2_mu = *(V_tilde+12*iy+j);
		j = lookup_HYP_3(nu, mu, ord);
		*u3_nu = *(V_tilde+12*iy+j);
	}
	
	iy = idn[ix][nu];
	iz = iup[ix][mu];
	j = lookup_HYP_3(nu, mu, ord);
	if((iy>=VOLUME) && (iy<(VOLUME+BNDRY)) && (iz<VOLUME)){
		iy = idn[iz][nu];
		u_comm_HYP('d', iy, nu, nu, j, u4_nu, V_tilde);
	}
	else if((iz>=VOLUME) && (iz<(VOLUME+BNDRY)) && (iy<VOLUME)){
		iy = iup[iy][mu];
		u_comm_HYP('u', iy, mu, nu, j, u4_nu, V_tilde);
	}
	else if(((iy>=VOLUME) && (iy<(VOLUME+BNDRY))) && ((iz>=VOLUME) && (iz<(VOLUME+BNDRY)))){
		iy = iup[map[iy-VOLUME]][mu];
		u_comm_HYP('L', iy, mu, nu, j, u4_nu, V_tilde);
	}
	else{
		iy = iup[iy][mu];
		*u4_nu = *(V_tilde+12*iy+j);
	}
	
	_su3_dagger((*u3_nu_dagger), (*u3_nu));
	_su3_times_su3((*C1), (*u3_nu_dagger), (*u2_mu));
	_su3_times_su3((*C3), (*C1), (*u4_nu));

	
	_su3_add((*C1), (*C2), (*C3));		
	_su3_add((*C), (*C), (*C1));
 }
 	
 *C1 = *u1[ix][mu];
 _su3_mul((*C1), 1.0-alpha1, (*C1));
 
 _su3_mul((*C), alpha1/6.0, (*C));
 _su3_add((*C), (*C), (*C1));
 
 
 afree(C1);
 afree(C2);
 afree(C3);
 afree(u1_mu);
 afree(u2_mu);
 afree(u1_nu);
 afree(u2_nu);
 afree(u3_nu);
 afree(u4_nu);
 afree(u2_nu_dagger);
 afree(u3_nu_dagger);

}


/************************************/
/*   is needed for stout-smearing	*/
/************************************/
/*static void Omega_mu(int ix, int mu){

 su3_dble *u1, *u2;
 
 u1 = pud[ix][mu];
 u2 = amalloc(sizeof(su3_dble), 4);
 
 _su3_dagger((*u2), (*u1));
 _su3_times_su3((*Omega), (*C), (*u2));
 
 afree(u2);
}*/


/************************************/
/*   is needed for stout-smearing	*/
/************************************/
/*static void Q_mu(void){
 
 const double c0 = -1, c1 = -1/3, c2 = 1/2;
 su3_dble *u1, *u2, *u3;
 
 u1 = amalloc(sizeof(su3_dble), 4);
 u2 = amalloc(sizeof(su3_dble), 4);
 u3 = amalloc(sizeof(su3_dble), 4);

 _su3_dagger((*u1), (*Omega));
 _su3_mul((*u2), c0, (*Omega));
 
 _su3_add((*u3), (*u1), (*u2));
 
 _su3_trace((*u1), (*u3));
 _su3_mul((*u1), c1, (*u1));
 
 _su3_add((*u2), (*u3), (*u1));
 _su3_i_mul((*u3), (*u2));
 _su3_mul((*Q), c2, (*u3));
 
 afree(u1);
 afree(u2);
 afree(u3);

}*/


/****************************/
/*	  converts a su3 to		*/
/*		  a su3_dble		*/
/****************************/
extern void su3_dble_to_su3(su3_dble *um1, su3 *um2){

 double *p_d;
 float *p_f;
 int i=0;
 p_d = (double*)um1;
 p_f = (float*)um2;

 for(i=0; i<18; i++){
	p_f[i] = (float)(p_d[i]);
 }
}


/****************************/
/*	converts a su3_dble to	*/
/*			a su3			*/
/****************************/
extern void su3_to_su3_dble(su3_dble *um1, su3 *um2){

 double *p_d;
 float *p_f;
 int i=0;
 p_d = (double*)um1;
 p_f = (float*)um2;

 for(i=0; i<18; i++){
	p_d[i] = (double)(p_f[i]);
 }
}


/********************************/
/* writes any double link fields*/
/*	 in the log file. This is	*/
/*		indipendent of the		*/
/*		 parallelisation		*/
/********************************/
extern void pud_output(su3_dble *doublelinks[VOLUME][4]){

 su3_dble *u1, *u2;
 double *p;
 int i, x, y, z, t, mu, ort[4], ix, ip, my_rank, send[5], recv[5];
 
 MPI_Status status;
 MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
 
 u1 = amalloc(sizeof(su3_dble), 4);
 u2 = amalloc(sizeof(su3_dble), 4);
 p = amalloc(sizeof(double), 4);
 
 for(x=0; x<NPROC1*L1; x++){
	for(y=0; y<NPROC2*L2; y++){
		for(z=0; z<NPROC3*L3; z++){
			for(t=0; t<NPROC0*L0; t++){
				ort[0] = t;
				ort[1] = x;
				ort[2] = y;
				ort[3] = z;
				ipt_global(ort, &ip, &ix);
				
				for(mu=0; mu<4; mu++){
					if(my_rank == ip){
						*u1 = *doublelinks[ix][mu];
						MPI_Send((double*)u1, 18, MPI_DOUBLE, 0, 37, MPI_COMM_WORLD);
					}
					if(my_rank == 0) MPI_Recv((double*)u2, 18, MPI_DOUBLE, ip, 37, MPI_COMM_WORLD, &status);
					if(my_rank == ip){
						send[0] = t;
						send[1] = x;
						send[2] = y;
						send[3] = z;
						send[4] = mu;
						MPI_Send(send, 5, MPI_INT, 0, 37, MPI_COMM_WORLD);
					}
					if(my_rank == 0) MPI_Recv(recv, 5, MPI_INT, ip, 37, MPI_COMM_WORLD, &status);

					if(my_rank == 0){
						
						p = (double*)u2;
						
						for(i=0; i<18; i++){
							message("t=%d, x=%d, y=%d, z=%d, mu=%d, link=%.10lf\n", recv[0], recv[1], recv[2], recv[3], recv[4], p[i]);
						}
					}
				}
			}
		}
	}
 }

 afree(u1);
 afree(u2);
 afree(p);

}


/********************************/
/* Compares two links to each	*/
/*		  other. This is		*/
/*		indipendent of the		*/
/*		 parallelisation		*/
/********************************/
extern void cmp_2pud(su3_dble *doublelinks1[VOLUME][4], su3_dble *doublelinks2[VOLUME][4]){

 su3_dble *u1_s, *u1_r, *u2_s, *u2_r;
 double *p1, *p2;
 int i, x, y, z, t, mu, ort[4], ix, ip, my_rank, send[5], recv[5];
 
 MPI_Status status;
 MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
 
 u1_s = amalloc(sizeof(su3_dble), 4);
 u1_r = amalloc(sizeof(su3_dble), 4);
 u2_s = amalloc(sizeof(su3_dble), 4);
 u2_r = amalloc(sizeof(su3_dble), 4);
 p1 = amalloc(sizeof(double), 4);
 p2 = amalloc(sizeof(double), 4);
 
 for(x=0; x<NPROC1*L1; x++){
	for(y=0; y<NPROC2*L2; y++){
		for(z=0; z<NPROC3*L3; z++){
			for(t=0; t<NPROC0*L0; t++){
				ort[0] = t;
				ort[1] = x;
				ort[2] = y;
				ort[3] = z;
				ipt_global(ort, &ip, &ix);
				
				for(mu=0; mu<4; mu++){
					if(my_rank == ip){
						*u1_s = *doublelinks1[ix][mu];
						MPI_Send((double*)u1_s, 18, MPI_DOUBLE, 0, 37, MPI_COMM_WORLD);
					}
					if(my_rank == 0) MPI_Recv((double*)u1_r, 18, MPI_DOUBLE, ip, 37, MPI_COMM_WORLD, &status);
					
					if(my_rank == ip){
						*u2_s = *doublelinks2[ix][mu];
						MPI_Send((double*)u2_s, 18, MPI_DOUBLE, 0, 37, MPI_COMM_WORLD);
					}
					if(my_rank == 0) MPI_Recv((double*)u2_r, 18, MPI_DOUBLE, ip, 37, MPI_COMM_WORLD, &status);
					
					if(my_rank == ip){
						send[0] = t;
						send[1] = x;
						send[2] = y;
						send[3] = z;
						send[4] = mu;
						MPI_Send(send, 5, MPI_INT, 0, 37, MPI_COMM_WORLD);
					}
					if(my_rank == 0) MPI_Recv(recv, 5, MPI_INT, ip, 37, MPI_COMM_WORLD, &status);

					if(my_rank == 0){
						
						p1 = (double*)u1_r;
						p2 = (double*)u2_r;
						
						for(i=0; i<18; i++){
							if(sqrt((p1[i]-p2[i])*(p1[i]-p2[i]))>0.001)
							message("t=%d, x=%d, y=%d, z=%d, mu=%d, link_diff=%.10lf\n", 
								recv[0], recv[1], recv[2], recv[3], recv[4], sqrt((p1[i]-p2[i])*(p1[i]-p2[i])));
						}
					}
				}
			}
		}
	}
 }

 afree(u1_s);
 afree(u1_r);
 afree(u2_s);
 afree(u2_r);
 afree(p1);
 afree(p2);

}


/********************************/
/* writes any single link fields*/
/*	 the output file. This is	*/
/*		indipendent of the		*/
/*		 parallelisation		*/
/********************************/
extern void pu_output(su3 *singlelinks[VOLUME][4]){

 su3 *u1, *u2;
 float *p;
 int i, x, y, z, t, mu, ort[4], ix, ip, my_rank, send[5], recv[5];
 
 MPI_Status status;
 MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
 
 u1 = amalloc(sizeof(su3), 4);
 u2 = amalloc(sizeof(su3), 4);
 p = amalloc(sizeof(float), 4);
 
 for(x=0; x<NPROC1*L1; x++){
	for(y=0; y<NPROC2*L2; y++){
		for(z=0; z<NPROC3*L3; z++){
			for(t=0; t<NPROC0*L0; t++){
				ort[0] = t;
				ort[1] = x;
				ort[2] = y;
				ort[3] = z;
				ipt_global(ort, &ip, &ix);
				
				for(mu=0; mu<4; mu++){
					if(my_rank == ip){
						*u1 = *singlelinks[ix][mu];
						MPI_Send((float*)u1, 18, MPI_FLOAT, 0, 37, MPI_COMM_WORLD);
					}
					if(my_rank == 0) MPI_Recv((float*)u2, 18, MPI_FLOAT, ip, 37, MPI_COMM_WORLD, &status);
					if(my_rank == ip){
						send[0] = t;
						send[1] = x;
						send[2] = y;
						send[3] = z;
						send[4] = mu;
						MPI_Send(send, 5, MPI_INT, 0, 37, MPI_COMM_WORLD);
					}
					if(my_rank == 0) MPI_Recv(recv, 5, MPI_INT, ip, 37, MPI_COMM_WORLD, &status);

					if(my_rank == 0){
						
						p = (float*)u2;
						
						for(i=0; i<18; i++){
							message("t=%d, x=%d, y=%d, z=%d, mu=%d, link=%.10lf\n", recv[0], recv[1], recv[2], recv[3], recv[4], p[i]);
						}
					}
				}
			}
		}
	}
 }

 afree(u1);
 afree(u2);
 afree(p);

}


/************************************/
/*	checks if the staples function	*/
/*	and C_mu do the same thing on	*/
/*		    one processor			*/
/************************************/
/*static void staples_check(int ix, int mu){

 su3_dble *um1;
 su3 *um1_f;
 double *p1, *p2;
 float *p3;
 int i;
 
 um1_f = amalloc(sizeof(su3), 4);
 um1 = amalloc(sizeof(su3_dble), 4);

 staples(ix, mu, um1_f);
 su3_to_su3_dble(um1, um1_f);
 project_to_su3_dble(um1);
 
 p1 = (double*)um1;
 p2 = (double*)C;
 p3 = (float*)um1_f;
 
 for(i=0; i<18; i++){ 
	if((p1[i]-p2[i]) > 0.000001) message("The difference between staples and C_mu is: %.7lf \n", p1[i]-p2[i]);
 }
 
 afree(um1_f);
 afree(um1);

}
*/

/************************************/
/* Copies the dummy variable pud_APE*/
/*	  into the pud link variables	*/
/************************************/
extern void pud_copy(su3_dble *u1[VOLUME][4], su3_dble *u2[VOLUME][4]){

 int ix, mu;
 
 for(ix=0; ix<VOLUME; ix++){
	for(mu=0; mu<4; mu++){
		*u2[ix][mu] = *u1[ix][mu];
	}
 }
  
 copy_bnd_ud_reverse(u2);
 free_ucom_bufs_reverse();
}


/************************************/
/*		 APE-smearing routine		*/
/************************************/
extern void APE_smearing(su3_dble *u1[VOLUME][4], su3_dble *u2[VOLUME][4], smearparm parm){

 int ibc, ix=0, mu=0;
 double alpha1 = parm.alpha[0];

 ibc=query_flags(BCD_FLIPPED);
 if (ibc==1)
   flipbcd();

 copy_bnd_ud();

 for(ix=0; ix<VOLUME; ix++){
	for(mu=0; mu<4; mu++){ 

		C = u2[ix][mu];
		
		C_mu(ix, mu, u1, alpha1);
		
		approx_project_to_su3_dble(C, 7);

	}
 }

 if (ibc==1)
   flipbcd();

}


/************************************/
/*	HYP-smearing routine	    */
/************************************/
extern void HYP_smearing(su3_dble *u1[VOLUME][4], su3_dble *u2[VOLUME][4], smearparm parm){

 int ibc, mu, nu, roh, eta, ix, i, j, *ord;
 su3_dble *V_bar, *V_tilde;
 double alpha1 = parm.alpha[0], alpha2 = parm.alpha[1], alpha3 = parm.alpha[2];

 ibc=query_flags(BCD_FLIPPED);
 if (ibc==1)
   flipbcd();

 V_bar = amalloc(12*VOLUME*sizeof(su3_dble), ALIGN);
 ord = amalloc(48*sizeof(int), ALIGN);

 copy_bnd_ud();
 
 /*lookuptable for V_bar in C_HYP_2*/
 j=0;
 for(mu=0; mu<4; mu++){
   for(i=1; i<4; i++){
     nu=(mu+i)&0x3;
     for(roh=(nu+1); roh<4; roh++){
	if(roh==mu) continue;
	  *(ord+4*j) = mu; 
	  *(ord+4*j+1) = nu; 
	  *(ord+4*j+2) = roh; 
	  *(ord+4*j+3) = j; 
	  j++;
	}
     }
 }
 
 for(ix=0; ix<VOLUME; ix++){
	j=0;
	for(mu=0; mu<4; mu++){
		for(i=1; i<4; i++){
			nu=(mu+i)&0x3; /*the 'bitwise AND' takes care that mu and nu are different*/
			for(roh=(nu+1); roh<4; roh++){
				if(roh==mu) continue;
				eta=6-(mu+nu+roh); /*now eta is the complement to mu, nu and roh*/
				C = (V_bar+j+ix*12);
				C_HYP_1(ix, mu, eta, u1, alpha3);
				approx_project_to_su3_dble(C, 4);
				j++;
			}
		}
	}
 }
 
 V_tilde = amalloc(12*VOLUME*sizeof(su3_dble), ALIGN);
 
 for(ix=0; ix<VOLUME; ix++){
	j=0;
	for(mu=0; mu<4; mu++){
		for(i=1; i<4; i++){
			nu=(mu+i)&0x3;
			C = (V_tilde+j+ix*12);
			C_HYP_2(ix, mu, nu, ord, V_bar, u1, alpha2);
			approx_project_to_su3_dble(C, 4);
			j++;
		}
	}
 }

 afree(V_bar);
 afree(ord);
 
 ord = amalloc(36*sizeof(int), ALIGN);
 
 /*lookuptable for V_tilde in C_HYP_3*/
 j=0;
 for(mu=0; mu<4; mu++){
	for(i=1; i<4; i++){
		nu=(mu+i)&0x3;
		*(ord+3*j) = mu; 
		*(ord+3*j+1) = nu; 
		*(ord+3*j+2) = j; 
		j++;
	}
 }

 for(ix=0; ix<VOLUME; ix++){
	for(mu=0; mu<4; mu++){
		C = u2[ix][mu];
		C_HYP_3(ix, mu, ord, V_tilde, u1, alpha1);
		approx_project_to_su3_dble(C, 4);
	}
 }
 
 afree(V_tilde);
 afree(ord);

 if (ibc==1)
   flipbcd();

}

static complex_dble su3_determinant(su3_dble u){

  complex_dble det;
  su3_vector_dble *v1, *v2, *v3, v4;

  v1=(su3_vector_dble*)(&u);
  v2=v1+1;
  v3=v1+2;

  _vector_cross_prod(v4,*v1,*v2);

  det.re = _vector_prod_re(v4,*v3);
  det.im = _vector_prod_im(v4,*v3);

  return det;

}



extern void check_SU3(su3_dble *u1[VOLUME][4]){

  int ibc,my_rank,id, ixx, mumu, sendbuffer[7], recvbuffer[7], testdet, testtrace;
  double sendbuffer_d[4], recvbuffer_d[4];
  complex_dble trace,det;
  su3_dble u2, u3;
  MPI_Status status;

  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

 ibc=query_flags(BCD_FLIPPED);
 if (ibc==1)
   flipbcd();

  message("Testing su3 matrices\n");

  for(ixx=0;ixx<VOLUME;ixx++){
  for(mumu=0;mumu<4;mumu++){


  _su3_dagger(u2,*(u1[ixx][mumu]));
  _su3_times_su3(u3,u2,*(u1[ixx][mumu]));
  trace.re = u3.c11.re+u3.c22.re+u3.c33.re +
             u3.c12.re+u3.c23.re+u3.c31.re +
             u3.c13.re+u3.c21.re+u3.c32.re;
  trace.im = u3.c11.im+u3.c22.im+u3.c33.im +
             u3.c12.im+u3.c23.im+u3.c31.im +
             u3.c13.im+u3.c21.im+u3.c32.im;

  det = su3_determinant(*(u1[ixx][mumu]));

  testtrace=0;
  if( (fabs(trace.re-3.0)>0.000001) || (fabs(trace.im)>0.000001) ) testtrace=1;
  testdet=0;
  if( (fabs(det.re-1.0)>0.000001) || (fabs(det.im)>0.000001) ) testdet=1;

  for(id=0;id<NPROC;id++){
    MPI_Barrier(MPI_COMM_WORLD);
    if( (my_rank == id) ){
      sendbuffer[0] = coords[ixx].t;
      sendbuffer[1] = coords[ixx].x;
      sendbuffer[2] = coords[ixx].y;
      sendbuffer[3] = coords[ixx].z;
      sendbuffer[4] = testdet;  
      sendbuffer[5] = testtrace;
      sendbuffer[6] = mumu;

      sendbuffer_d[0] = trace.re ;
      sendbuffer_d[1] = trace.im ;
      sendbuffer_d[2] = det.re;
      sendbuffer_d[3] = det.im;
      MPI_Send( sendbuffer, 7, MPI_INT, 0, 37, MPI_COMM_WORLD );
    }
    if(my_rank==0)
      MPI_Recv( recvbuffer, 7, MPI_INT, id, 37, MPI_COMM_WORLD, &status);

    if(my_rank==id)
      MPI_Send( sendbuffer_d, 4, MPI_DOUBLE, 0, 38, MPI_COMM_WORLD );
    if(my_rank==0)
      MPI_Recv( recvbuffer_d, 4, MPI_DOUBLE, id, 38, MPI_COMM_WORLD, &status);
    
    if(my_rank==0){
      if( (recvbuffer[4]==1)  || (recvbuffer[5]==1) ){
        message("\nThe matrix at (%d,%d,%d,%d) in %d-direction is not in su(3), dettest=%d, tracetest=%d\n",
            recvbuffer[0],recvbuffer[1],recvbuffer[2],recvbuffer[3],recvbuffer[6],recvbuffer[4],recvbuffer[5]);
        message("det.re=%lf, det.im=%lf, trace.re=%lf, trace.im=%lf\n",
                 recvbuffer_d[2],recvbuffer_d[3],recvbuffer_d[0],recvbuffer_d[1]);
        fflush(stdout);
      }
    }
  }


  }}

 if (ibc==1)
   flipbcd();

}






