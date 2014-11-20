/*******************************************************************************
*
* File source.c
*
* Copyright (C) 2008-9 Andreas Juettner, Eric Endress, Bastian Knippschild
* 
* The externally accessible functions are:
* 
* extern void random_Z4(int vol,spinor_dble *pk)
*	generate a Z2\timesZ2 random vector
*
* extern void SourceRadius(srcdef src, spinor_dble *sd)
* 	compute width of gaussian distributed source
*
* extern void jacobi(spinor_dble *sd, srcdef src, int idirac, int icolor)
*	compute the jacobi smeared source
*
* extern void srcfld(spinor_dble *rnd, spinor_dble *sd, srcdef src, int idirac,
*	int icolor)
*	generates quark sources of type
*		pt	= point
*		J_thin	= jacobi
*		J_APE = jacobi with APE-smearing
*		J_HYP = jacobi with HYP-smearing
*		Z4	= Z4
*		Z4spin	= Z4 spin diluted
*		Z2	= Z2
*
* extern void extended_src(int igamma, int *mom, int t,
*	full_spinor_dble *k, full_spinor_dble *l)
* 	compute extended source (sequential source)
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "mpi.h"
#include "start.h"
#include "global.h"
#include "measure.h"
#include "linalg.h"
#include "update.h"
#include "random.h"

#define SOURCE_C

/*  c = a * b */
#define CMUL(a,b,c) { (c).re = (a).re*(b).re - (a).im*(b).im; \
                      (c).im = (a).re*(b).im + (a).im*(b).re; }
static void pointsource(spinor_dble *sd, srcdef src, int idirac, int icolor)
{
   int x[4],ix,ip;
   int my_rank;
   double *s;
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   if(my_rank==0) 
	message("Creating point source\n");
   x[0]=src.pos[0];
   x[1]=src.pos[1];
   x[2]=src.pos[2];
   x[3]=src.pos[3];
   
   ipt_global(x,&ip,&ix);
   if (my_rank==ip)
   {
      s=(double*)(sd+ix);
      s[6*(idirac-1)+2*(icolor-1)]=1.0;
   }

}
extern void random_Z4(int vol,spinor_dble *pk)
{
   int i;
   double r[24],*s,norm,neg,pos;
   spinor_dble *rpk;
   norm = sqrt(L1*NPROC1*L2*NPROC2*L3*NPROC3);
   neg  = -1.0/sqrt(2)/norm;
   pos  =  1.0/sqrt(2)/norm;
   for (rpk=pk;rpk<(pk+vol);rpk++)
   {
      ranlxd(r,24);
      s=(double*)(rpk);
      for (i=0;i<24;i++)
      {
       if(r[i]<0.5)
         s[i]=neg; 
       else
         s[i]=pos; 
      }
   }
}

extern void SourceRadius(srcdef src, spinor_dble *sd)
{
  /*****************************************************************************
  *                                                                     
  *                    Radius calculation                               
  *                                                                     
  *       To calculate the radius of the source the                     
  *       sourceorigin is shifted to the point (0,0,0).        
  *       Only points for (x1, x2, x3)<(L1/2, L2/2, L3/2) hold 
  *	are used because of the periodic boundary conditions.
  *                                                            
  *                       ++++++++++++++                       
  *                       +            +                       
  *                       +            +                       
  *                   L/2 +++++++      +                       
  *                       +     +      +                       
  *                       +     +      +                       
  *                 (0,0) ++++++++++++++                       
  *                             L/2                            
  *****************************************************************************/
  
  int my_rank, t, x, y, z, i, ix;
  int x_send[NPROC], y_send[NPROC], z_send[NPROC];
  double S2x2_send[NPROC], S2_send[NPROC], F_send[NPROC];
  double F=0, S2=0, S2x2=0, rmsradius=0, r;
  spinor_dble *psi;
 
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  
  for(i=0; i<NPROC; i++)
  {
    x_send[i]=0;
    y_send[i]=0;
    z_send[i]=0;
    S2x2_send[i]=0.0;
    S2_send[i]=0.0;
    F_send[i]=0.0;
  }
  
  for(ix=0; ix<VOLUME; ix++)
  {
    t=coords[ix].t;
    x=coords[ix].x;
    y=coords[ix].y;
    z=coords[ix].z;
  
    F=0;
    psi=sd+ix;
  
    F=_vector_prod_re((*psi).c1, (*psi).c1)+
      _vector_prod_re((*psi).c2, (*psi).c2)+
      _vector_prod_re((*psi).c3, (*psi).c3)+
      _vector_prod_re((*psi).c4, (*psi).c4);
  
    x -= src.pos[1];
    y -= src.pos[2];
    z -= src.pos[3];

    if(x<0) x += NPROC1*L1;
    if(y<0) y += NPROC2*L2;
    if(z<0) z += NPROC3*L3;
  
  if((x<=(NPROC1*L1/2)) && (y<=(NPROC2*L2/2)) && (z<=(NPROC3*L3/2))){ 
      r = (x*x+y*y+z*z);
      S2x2 += r*F;
      S2 += F;
    }
  
  }
  MPI_Gather(&S2x2, 1, MPI_DOUBLE, S2x2_send, 1, MPI_DOUBLE, 0,MPI_COMM_WORLD);
  MPI_Gather(&S2, 1, MPI_DOUBLE, S2_send, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  S2x2=0;
  S2=0;
  
  if(my_rank==0)
  {
    for(i=0; i<NPROC; i++)
    {
       S2x2 += S2x2_send[i];
       S2 += S2_send[i];
    }
  }
  rmsradius=sqrt(S2x2/S2);

  if(my_rank == 0)
    message("the rms radius of the source is %lf \n", rmsradius);
}
  

extern void jacobi_sink(spinor_dble *sd, srcdef src, int idirac, int icolor){

  int i=0, ix=0, iy=0, n=0, mu=0, my_rank, ip1=0, ip2=0, iw=0,shift[4];
  su3_dble *up, *um;
  spinor_dble *chi, *psi, *phi, **wsd0, **wsd1, **wsd2, *p, *(*p0)[NSPIN];
  const spinor_dble sd0={{{0.0}}};
  double zw[18]={0},d;
  MPI_Status status1;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  double norm = 1/(1.+6.*src.kappa); 

  p0=amalloc(3*sizeof(*p0),3);
  p=amalloc(3*NSPIN*sizeof(spinor_dble),ALIGN);

  shift[0] = -src.pos[0];
  shift[1] = -src.pos[1];
  shift[2] = -src.pos[2];
  shift[3] = -src.pos[3];
  shift_ud(shift);
  alloc_pud_sm1();

  /*defines the way of linksmearing*/
  if(strcmp(src.type,"J_thin")==0) 
      	 pud_copy(pud, pud_sm1);
  else if (strcmp(src.type,"J_APE")==0) 
      	 APE_smearing(pud, pud_sm1, src.smear);
  else if (strcmp(src.type,"J_HYP")==0)
      	 HYP_smearing(pud, pud_sm1, src.smear);
      	
  for(i=0; i<3; i++)
  {
    for (ix=0;ix<NSPIN;ix++)
    {
       p0[i][ix]=p;
       *p=sd0;
       p+=1;
    }
  }
  
  wsd0=&(p0[0][0]);
  wsd1=&(p0[1][0]);
  wsd2=&(p0[2][0]);

  if(my_rank==0) { 
      	message("Jacobi-smearing with n=%d and kappa=%.2f\n",
      			  src.n,src.kappa);
      if (strcmp(src.type,"J_APE")==0) 
      	message("and APE-smearing with alpha = %.2f\n",
      			  src.smear.alpha[0]);
      else if (strcmp(src.type,"J_HYP")==0)
      	message("and HYP-smearing with alpha1 = %.2f, alpha2 = %.2f and alpha3 = %.2f \n",
      			  src.smear.alpha[0], src.smear.alpha[1],src.smear.alpha[2]);
      }
      
  copy_bnd_ud();

  for(n=0; n<src.n; n++)
  {
  if(n>0){
  
    for(ix=0; ix<VOLUME; ix++) *wsd0[ix]=*wsd2[ix];
  
  }
  
  for(ix=0; ix<VOLUME; ix++)
  {
   /* if (coords[ix].t==src.pos[0])
      {*/
        phi=wsd1[ix];
        if(n==0) *wsd2[ix]=sd[ix];
        psi=wsd2[ix];
        for(mu=1; mu<4; mu++)
        {
          iy=idn[ix][mu];
          if((iy>=VOLUME) && (iy<(VOLUME+BNDRY)))
      		{
            ip1=npr[2*mu];
            ip2=npr[2*mu+1];
            iw=map[iy-VOLUME];
  
            if(n==0) *wsd0[iw]=sd[iw];
            *phi=*wsd0[iw];
            chi=wsd1[iw];
  
            um=pud_sm1[iw][mu];
  
            MPI_Sendrecv((double*)(um), 18, MPI_DOUBLE, ip2, 37,
               zw, 18, MPI_DOUBLE, ip1, 37, MPI_COMM_WORLD, &status1);
  
            um=(su3_dble*)zw;
            MPI_Sendrecv((double*)(phi), 24, MPI_DOUBLE, ip2, 37,
               (double*)(chi), 24, MPI_DOUBLE, ip1, 37,
               MPI_COMM_WORLD, &status1);
  
          }
          else 
          {
            um=pud_sm1[iy][mu];
            if(n==0) *wsd0[iy]=sd[iy];
            chi=wsd0[iy];
          }
          _su3_inverse_multiply((*phi).c1,(*um),(*chi).c1);
          _su3_inverse_multiply((*phi).c2,(*um),(*chi).c2);
          _su3_inverse_multiply((*phi).c3,(*um),(*chi).c3);
          _su3_inverse_multiply((*phi).c4,(*um),(*chi).c4);
          _vector_add_assign((*psi).c1, (*phi).c1);
          _vector_add_assign((*psi).c2, (*phi).c2);
          _vector_add_assign((*psi).c3, (*phi).c3);
          _vector_add_assign((*psi).c4, (*phi).c4);
  
          iy=iup[ix][mu];
  
          if((iy>=VOLUME) && (iy<(VOLUME+BNDRY)))
      		{
            ip1=npr[2*mu];
            ip2=npr[2*mu+1];
            iw=map[iy-VOLUME];
  
            if(n==0) *wsd0[iw]=sd[iw];
            *phi=*wsd0[iw];
            chi=wsd1[iw];
  
            MPI_Sendrecv((double*)(phi), 24, MPI_DOUBLE, ip1, 37,
                    (double*)(chi), 24, MPI_DOUBLE, ip2, 37,
                    MPI_COMM_WORLD, &status1);
          }
          else {
            if(n==0) *wsd0[iy]=sd[iy];
            chi=wsd0[iy];
          }
  
          up=pud_sm1[ix][mu];
  
          _su3_multiply((*phi).c1,(*up),(*chi).c1);
          _su3_multiply((*phi).c2,(*up),(*chi).c2);
          _su3_multiply((*phi).c3,(*up),(*chi).c3);
          _su3_multiply((*phi).c4,(*up),(*chi).c4);
          _vector_add_assign((*psi).c1, (*phi).c1);
          _vector_add_assign((*psi).c2, (*phi).c2);
          _vector_add_assign((*psi).c3, (*phi).c3);
          _vector_add_assign((*psi).c4, (*phi).c4);
        }
        _vector_mul((*psi).c1, src.kappa, (*psi).c1);
        _vector_mul((*psi).c2, src.kappa, (*psi).c2);
        _vector_mul((*psi).c3, src.kappa, (*psi).c3);
        _vector_mul((*psi).c4, src.kappa, (*psi).c4);
  
        *wsd1[ix]=sd[ix];
        phi=wsd1[ix];
        _vector_add_assign((*psi).c1, (*phi).c1);
        _vector_add_assign((*psi).c2, (*phi).c2);
        _vector_add_assign((*psi).c3, (*phi).c3);
        _vector_add_assign((*psi).c4, (*phi).c4);


        _spinor_mul(*wsd2[ix],norm,*wsd2[ix]);
  
   /*   }*/
    }
    if(n==(src.n-1)){ 

      d=norm_square_dble(VOLUME,1,wsd2[0]);
      d=(1/sqrt(d));
      message("jacobi-smeared vector normalized with %.2e\n",d);
      for(ix=0; ix<VOLUME; ix++){
      		sd[ix]=*wsd2[ix];
      /* _spinor_mul(sd[ix],d,sd[ix]);*/ 
     }
      d=norm_square_dble(VOLUME,1,sd);
      message("normalized jacobi-smeared vector %.2e\n",d);
    }
  }
  SourceRadius(src,sd);
  shift[0] = +src.pos[0];
  shift[1] = +src.pos[1];
  shift[2] = +src.pos[2];
  shift[3] = +src.pos[3];
  shift_ud(shift);
  fflush(stdout);
  afree(p);
  afree(p0[0][0]);
  afree(p0);
  p0=NULL;
  free_pud_sm1();
}


extern void jacobi_source(spinor_dble *sd, srcdef src, int idirac, int icolor){

 int i=0, ix=0, iy=0, n=0, mu=0, my_rank, ip1=0, ip2=0, iw=0;
    su3_dble *up, *um;
    spinor_dble *chi, *psi, *phi, **wsd0, **wsd1, **wsd2, *p, *(*p0)[NSPIN];
    const spinor_dble sd0={{{0.0}}};
    double zw[18]={0},d;
    MPI_Status status1;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    double norm = 1/(1.+6.*src.kappa); 

    copy_bnd_ud();

    p0=amalloc(3*sizeof(*p0),3);
    p=amalloc(3*NSPIN*sizeof(spinor_dble),ALIGN);

    alloc_pud_sm1();

    /*defines the way of linksmearing*/
    if(strcmp(src.type,"J_thin")==0) 
		 pud_copy(pud, pud_sm1);
    else if (strcmp(src.type,"J_APE")==0) 
		 APE_smearing(pud, pud_sm1, src.smear);
    else if (strcmp(src.type,"J_HYP")==0)
		 HYP_smearing(pud, pud_sm1, src.smear);

    copy_bnd_ud_reverse(pud_sm1);
		
    for(i=0; i<3; i++)
    {
      for (ix=0;ix<NSPIN;ix++)
      {
         *p=sd0;
         p0[i][ix]=p;
         p+=1;
      }
    }
    
    wsd0=&(p0[0][0]);
    wsd1=&(p0[1][0]);
    wsd2=&(p0[2][0]);

    if(my_rank==0) { 
		message("Jacobi-smearing with n=%d and kappa=%.2f\n",
				  src.n,src.kappa);
	if (strcmp(src.type,"J_APE")==0) 
		message("and APE-smearing with alpha = %.2f\n",
				  src.smear.alpha[0]);
	else if (strcmp(src.type,"J_HYP")==0)
		message("and HYP-smearing with alpha1 = %.2f, alpha2 = %.2f and alpha3 = %.2f \n",
				  src.smear.alpha[0], src.smear.alpha[1],src.smear.alpha[2]);
	}

    message("sourceposition = (%d, %d, %d, %d)\n", src.pos[0], src.pos[1], src.pos[2], src.pos[3]);

    for(n=0; n<src.n; n++)
    {
    if(n>0){
    
      for(ix=0; ix<VOLUME; ix++) *wsd0[ix]=*wsd2[ix];
    
    }
    
    for(ix=0; ix<VOLUME; ix++)
    {
      if (coords[ix].t==src.pos[0]){

          phi=wsd1[ix];
          if(n==0) *wsd2[ix]=sd[ix];
          psi=wsd2[ix];
          for(mu=1; mu<4; mu++){

            iy=idn[ix][mu];
            if((iy>=VOLUME) && (iy<(VOLUME+BNDRY))){

              ip1=npr[2*mu];
              ip2=npr[2*mu+1];
              iw=map[iy-VOLUME];
    
              if(n==0) *wsd0[iw]=sd[iw];
              *phi=*wsd0[iw];
              chi=wsd1[iw];
    
              um=pud_sm1[iw][mu];
    
              MPI_Sendrecv((double*)(um), 18, MPI_DOUBLE, ip2, 37,
                 zw, 18, MPI_DOUBLE, ip1, 37, MPI_COMM_WORLD, &status1);
    
              um=(su3_dble*)zw;
              MPI_Sendrecv((double*)(phi), 24, MPI_DOUBLE, ip2, 37,
                 (double*)(chi), 24, MPI_DOUBLE, ip1, 37,
                 MPI_COMM_WORLD, &status1);
    
            }
            else {
              um=pud_sm1[iy][mu];
              if(n==0) *wsd0[iy]=sd[iy];
              chi=wsd0[iy];
            }
            _su3_inverse_multiply((*phi).c1,(*um),(*chi).c1);
            _su3_inverse_multiply((*phi).c2,(*um),(*chi).c2);
            _su3_inverse_multiply((*phi).c3,(*um),(*chi).c3);
            _su3_inverse_multiply((*phi).c4,(*um),(*chi).c4);
            _vector_add_assign((*psi).c1, (*phi).c1);
            _vector_add_assign((*psi).c2, (*phi).c2);
            _vector_add_assign((*psi).c3, (*phi).c3);
            _vector_add_assign((*psi).c4, (*phi).c4);
    
            iy=iup[ix][mu];
    
            if((iy>=VOLUME) && (iy<(VOLUME+BNDRY)))
			{
              ip1=npr[2*mu];
              ip2=npr[2*mu+1];
              iw=map[iy-VOLUME];
    
              if(n==0) *wsd0[iw]=sd[iw];
              *phi=*wsd0[iw];
              chi=wsd1[iw];
    
              MPI_Sendrecv((double*)(phi), 24, MPI_DOUBLE, ip1, 37,
                      (double*)(chi), 24, MPI_DOUBLE, ip2, 37,
                      MPI_COMM_WORLD, &status1);
            }
            else {
              if(n==0) *wsd0[iy]=sd[iy];
              chi=wsd0[iy];
            }
    
            up=pud_sm1[ix][mu];
    
            _su3_multiply((*phi).c1,(*up),(*chi).c1);
            _su3_multiply((*phi).c2,(*up),(*chi).c2);
            _su3_multiply((*phi).c3,(*up),(*chi).c3);
            _su3_multiply((*phi).c4,(*up),(*chi).c4);
            _vector_add_assign((*psi).c1, (*phi).c1);
            _vector_add_assign((*psi).c2, (*phi).c2);
            _vector_add_assign((*psi).c3, (*phi).c3);
            _vector_add_assign((*psi).c4, (*phi).c4);
          }
          _vector_mul((*psi).c1, src.kappa, (*psi).c1);
          _vector_mul((*psi).c2, src.kappa, (*psi).c2);
          _vector_mul((*psi).c3, src.kappa, (*psi).c3);
          _vector_mul((*psi).c4, src.kappa, (*psi).c4);
    
          *wsd1[ix]=sd[ix];
          phi=wsd1[ix];
          _vector_add_assign((*psi).c1, (*phi).c1);
          _vector_add_assign((*psi).c2, (*phi).c2);
          _vector_add_assign((*psi).c3, (*phi).c3);
          _vector_add_assign((*psi).c4, (*phi).c4);

          _spinor_mul(*wsd2[ix],norm,*wsd2[ix]);
    
        }
      }
      if(n==(src.n-1)){ 
        
        d=norm_square_dble(VOLUME,1,wsd2[0]);
        d=1/sqrt(d);
        message("jacobi-smeared vector normalized with %.2e\n",d);
        for(ix=0; ix<VOLUME; ix++){
         sd[ix]=*wsd2[ix];
     /*    _spinor_mul(sd[ix],d,sd[ix]);*/
       }
        d=norm_square_dble(VOLUME,1,sd);
        message("normalized jacobi-smeared vector %.2e\n",d);
      }
    }
    SourceRadius(src,sd);
    afree(p);
    afree(p0[0][0]);
    afree(p0);
    p0=NULL;
    free_pud_sm1();
}


extern void srcfld(spinor_dble *rnd, spinor_dble *sd, srcdef src, int idirac,
	int icolor)
{
  int my_rank;
  int i;
  set_sd2zero(VOLUME,sd);
  /***************************************************************************/
  if(strcmp(src.type,"pt")==0)
  {
  /***************************************************************************/
   pointsource(sd,src,idirac,icolor);
  
  }
  /***************************************************************************/
  else if(strcmp(src.type,"Z4")==0){
  /***************************************************************************/
    int ix = 0;
    double *in,*out;
    spinor_dble *rpk;
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    if(my_rank==0) 
	message("Creating Z4 source on time-slice %d\n",src.pos[0]);

    if(icolor==1)
    {

            for( rpk = sd ; rpk < (sd+VOLUME) ; rpk ++)
            {
             if(coords[ix].t==src.pos[0])
             {
              in =(double *)(rnd+ix);
              out=(double *)rpk;
 	      for (i=0;i<24;i++){ 
               *(out+i)=*(in+i);
	      }
             }
             ix++;
            }
  }
  }
  else if(strcmp(src.type,"Z2")==0){
  /***************************************************************************/
    int ix = 0;
    double *in,*out;
    spinor_dble *rpk;
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    if(my_rank==0) 
	message("Creating Z2 source on time-slice %d\n",src.pos[0]);

    if(icolor==1)
    {

            for( rpk = sd ; rpk < (sd+VOLUME) ; rpk ++)
            {
             if(coords[ix].t==src.pos[0])
             {
              in =(double *)(rnd+ix);
              out=(double *)rpk;
 	      for (i=0;i<24;i=i+2){ 
               *(out+i)=*(in+i);
	      }
             }
             ix++;
            }
  }
  }
  /***************************************************************************/
  else if(strcmp(src.type,"Z4spin")==0){
  /***************************************************************************/
    int ix = 0;
    double *in,*out;
    spinor_dble *rpk;
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    if(my_rank==0) 
	message("Creating Z4spin source on time-slice %d\n",src.pos[0]);

    if(icolor==1)
    {

            for( rpk = sd ; rpk < (sd+VOLUME) ; rpk ++)
            {
             if(coords[ix].t==src.pos[0])
             {
              in =(double *)(rnd+ix);
              out=(double *)rpk;
 	      for (i=0;i<6;i++){ 
               *(out+(idirac-1)*6+i)=*(in+i);
	      }
             }
             ix++;
            }
  }
  }
  /***************************************************************************/

  /***************************************************************************/
  else if((strcmp(src.type,"J_thin")==0)||(strcmp(src.type,"J_HYP")==0)||(strcmp(src.type,"J_APE")==0))
  {
  /***************************************************************************/
	  

	  pointsource(sd, src, idirac, icolor);
       /*   transform_sd(sd,sd);
*/	  jacobi_source(sd, src, idirac, icolor);

  }		 
  /**************************************************************************/
  else{
    error_root(1,1,"srcfld [source.c]",
	"Couldn't find specified source type %s\n",src.type);   
  }
}

extern void extended_src(int igamma, int *mom, int t,
		full_spinor_dble *k, full_spinor_dble *l)
{  
   full_spinor_dble s2;
   int ix=0;
   int x0,x1,x2,x3;
   double PIdouble,phi;
   complex_dble phase;
   PIdouble=0.5*6.283185307179586476925286;
   set_fsd2zero(VOLUME,l);
   for (x0=0;x0<L0;x0++)
   { 
   for (x1=0;x1<L1;x1++)
    {
    for (x2=0;x2<L2;x2++)
     {
     for (x3=0;x3<L3;x3++)
      {
      if (cpr[0]*L0+x0==t){
      ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];
      phi=((2.0*PIdouble)*(
      ((double)((*(mom  ))*(cpr[1]*L1+x1))/(double)(NPROC1*L1)) +
      ((double)((*(mom+1))*(cpr[2]*L2+x2))/(double)(NPROC2*L2)) +
      ((double)((*(mom+2))*(cpr[3]*L3+x3))/(double)(NPROC3*L3))));
      phase.re = cos(phi);
      phase.im = sin(phi);
      mul_gamma_l(igamma,(k+ix),&s2);
      mul_cmplx_fsv_dble(&phase,&s2,(l+ix));
      }

      }
     }
    }
   }
}


