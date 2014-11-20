/*******************************************************************************
*
* File correlators.c
*
* Copyright (C) 2008-11 Andreas Juettner
*                       Bastian Knippschild
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* 2- and 3-point functions for mesons
*
* void C2(full_spinor_dble *k, full_spinor_dble *l, int num_c, int *mu,
*                int *momenta, complex_dble *sst);
*    Function that contracts the full spinors k and l for all mutual 
*    combinatinos of num_c gammas given in mu and for momenta
*    specified in momenta. The output is handed over to sst.
*
* void seq_src_test(full_spinor_dble *k, int gamma_insertion, int *pos)
*    Funtcion that does tests on the sequential propagator
* 
*
* 2- and 3-point functions for baryons
*
* void C2_b(full_spinor_dble *svk, full_spinor_dble *svl, full_spinor_dble *svm, 
*           complex_dble **s_pointer, int *momenta)
* 
* void C3_b(full_spinor_dble *sv1, full_spinor_dble *sv3, spinor_dble **wsd, 
*           complex_dble **s_pointer, int *momenta, int iw1, int iw3) 
*
*******************************************************************************/
#define CORRELATORS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "start.h"
#include "global.h"
#include "measure.h"
#include "linalg.h"
#include "mpi.h"
#include "random.h"
#include "misc.h"
#include "sw_term.h"
#include "dirac.h"
#include "version.h"
#include "version_qcd-measure.h"
#include "update.h"




static const complex_dble cplxd0={0};
static const full_spinor_dble fsd0={{{0}}};


/* generate momentum 3-vectors for given maximal component comp_max */
int gen_momenta(int comp_max,int *momenta)
{
  int i,ix,iy,iz;
  i=0;
  for(ix=-comp_max; ix<=comp_max; ix++)
  for(iy=-comp_max; iy<=comp_max; iy++)
  for(iz=-comp_max; iz<=comp_max; iz++){
   (*(momenta+i))   = ix;
   (*(momenta+i+1)) = iy;
   (*(momenta+i+2)) = iz;
   i=i+3;
  }
  return i/3;
}
static void flush_fs_double(full_spinor_dble *fs)
{
  int i;
  double *dum;
  dum=(double*)fs;
  for (i=0;i<288;i++) {
    *(dum+i)=0.0;
    *(dum+i)=0.0;
    }
}
static void flush_corr(int num_c,int num_mom, complex_dble *s)
{
  int i;
  for (i=0;i<num_c*num_c*num_mom;i++) {
    *(s+i)=cplxd0;
    }
}
static void flush_corr2(int max, complex_dble *s)
{
  int i;
  for (i=0;i<max;i++) {
    *(s+i)=cplxd0;
    }
}
static void add_corr(int num_c, int num_mom, complex_dble *si, complex_dble *so){
  int i;
  for (i=0;i<num_c*num_c*num_mom;i++){
   so[i].re += si[i].re;
   so[i].im += si[i].im;
   }
}
static void add_corr2(int max, complex_dble *si, complex_dble *so){
  int i;
  for (i=0;i<max;i++){
   so[i].re += si[i].re;
   so[i].im += si[i].im;
   }
}


extern void C2(full_spinor_dble *k, full_spinor_dble *l, int num_c, int *mu, 
		             int nmom_max, int *momenta, complex_dble *sst)
{
  int x0,x1,x2,x3,ix,T,pt;
  int i,ip,mu1,mu2,imu1,imu2,max,max2;
  int num_mom=nmom_max;

  double *phase,phase0,phase_sin,phase_cos,PI,wt1,wt2;
  double X,Y,Z,px,py,pz;
 
  complex_dble *st=NULL,*xt=NULL,*xhelp,*xhelp1,*xhelp2, corr;
  full_spinor_dble sk,sl,shelp1,shelp2;
    
  PI=0.5*6.283185307179586476925286;
  
  wt1=MPI_Wtime();
  
  max   = num_c*num_c*num_mom*NPROC0*L0;
  max2  = num_c*num_c*num_mom;

  st    = amalloc(max*sizeof(complex_dble),3);
  xt    = amalloc(max2*L0*sizeof(complex_dble),3);
  phase = amalloc(num_mom*sizeof(double),3);

  flush_corr2(max,st);
  flush_corr2(max2*L0,xt);
  flush_corr2(max,sst);  
  
  /**** shortcuts to save compute-time ****/
  T    =          (NPROC0*L0);
  X    = (double) (NPROC1*L1);
  Y    = (double) (NPROC2*L2);
  Z    = (double) (NPROC3*L3);
  pt   =          cpr[0]*L0;
  px   = (double) cpr[1]*L1;
  py   = (double) cpr[2]*L2;
  pz   = (double) cpr[3]*L3;
  
  for(x1=0;x1<L1;x1++){
    for(x2=0;x2<L2;x2++){
      for(x3=0;x3<L3;x3++){
      
        /**** calculating the phase for fouriertransformation ****/
        for(ip=0;ip<num_mom;ip++){
          phase[ip] = 2.0*PI*(
             ((double) (*(momenta+ip*3  ))*(px+x1))/X +
             ((double) (*(momenta+ip*3+1))*(py+x2))/Y +
             ((double) (*(momenta+ip*3+2))*(pz+x3))/Z);
        }
          
        for(x0=0;x0<L0;x0++){

          ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];
          i=0;

          xhelp = (xt + x0*max2);         
 
          /****  S'=gamma_0*gamma_5*adj(S)*gamma_5 ****/
          adj_full_spinor((k+ix), &sk);
          mul_gamma_l(5,&sk,&shelp1);
          mul_gamma_l(0,&shelp1,&sk);
          mul_gamma_r(5,&sk,&shelp1);

          /****  S''=S*gamma_0  ****/
          mul_gamma_r(0,(l+ix),&shelp2);

          for(imu1=0;imu1<num_c;imu1++){
            for(imu2=0;imu2<num_c;imu2++){
            
              mu1 = mu[imu1];
              mu2 = mu[imu2];
              
              /****  S=Gamma*S' ****/
              mul_gamma_l(mu2,&shelp1,&sk);
 
              /****  S=Gamma*S''  ****/
              mul_gamma_l(mu1,&shelp2,&sl);

              /****  Contraction  ****/
              meson_trace(&sk,&sl,&corr);
              
              /****  multiplication by  exp{ipy}  ****/
              for(ip=0;ip<num_mom;ip++){
                /***** performance ****/  
                phase0 = phase[ip];
                phase_sin = sin(phase0);
                phase_cos = cos(phase0);
                xhelp2 = (xhelp + i);
                /**********************/     
                (*xhelp2).re += corr.re*phase_cos - corr.im*phase_sin;
                (*xhelp2).im += corr.im*phase_cos + corr.re*phase_sin; 
                i=i+1;
              }
              
            }
          }
          
        }
      }
    }
  }

  for(x0=0;x0<L0;x0++){
    xhelp  = (xt + x0*max2);
    xhelp1 = (st + pt + x0);
    for(i=0;i<max2;i++){
      (*(xhelp1 + i*T)).re = (*(xhelp + i)).re;
      (*(xhelp1 + i*T)).im = (*(xhelp + i)).im;
    }
  }

  /**** combine the results from all processes in process 0 ****/
  MPI_Reduce(st,sst,2*max,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  wt2=MPI_Wtime();
  message("Time for meson contractions %.2e sec\n",wt2-wt1);
  fflush(stdout);

  afree(st);
  afree(phase);
  afree(xt); 
}


extern void seq_src_test(full_spinor_dble *k, int gamma_insertion, int *pos)
{
  int x0,x1,x2,x3,ix;
  int i,t,my_rank;

  complex_dble st1,sst1,*st2=NULL,*sst2=NULL;
  complex_dble corr1,corr2,dum;
 
  full_spinor_dble sk, shelp;
  
  st2  = amalloc(NPROC0*L0*sizeof(complex_dble),3);
  sst2 = amalloc(NPROC0*L0*sizeof(complex_dble),3);

  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

  st1.re=0.0;
  st1.im=0.0;
  sst1.re=0.0;
  sst1.im=0.0;

  for (i=0;i<NPROC0*L0;i++)
   {
    (*(st2+i)).re  = 0.0;
    (*(st2+i)).im  = 0.0;
    (*(sst2+i)).re = 0.0;
    (*(sst2+i)).im = 0.0;
   }

  for(x0=0;x0<L0;x0++){
  
    corr1.re=0.0;
    corr1.im=0.0;
    corr2.re=0.0;
    corr2.im=0.0;

    for(x1=0;x1<L1;x1++){
      for(x2=0;x2<L2;x2++){
        for(x3=0;x3<L3;x3++){
        
          ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];

          if(((cpr[0]*L0+x0)==pos[0])&&
             ((cpr[1]*L1+x1)==pos[1])&&
	            ((cpr[2]*L2+x2)==pos[2])&&
	            ((cpr[3]*L3+x3)==pos[3])){
      
            adj_full_spinor((k+ix), &shelp);
            mul_gamma_l(gamma_insertion,&shelp,&sk);
            fsv_trace(&sk,&dum);
            st1.re += dum.re; 
            st1.im += dum.im;
          }

          adj_full_spinor((k+ix), &sk);
          meson_trace(&sk,(k+ix),&dum);
          corr2.re += dum.re; 
          corr2.im += dum.im;

        }
      }
    }
    (*(st2+cpr[0]*L0 + x0)).re+=corr2.re;
    (*(st2+cpr[0]*L0 + x0)).im+=corr2.im;
  }
   
  MPI_Reduce(&st1,&sst1,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(st2,sst2,2*NPROC0*L0,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  if (my_rank==0){ 
  
    message("Tr(gamma_%d S_%d)|_(%d,%d,%d,%d) = %.6e + i %.6e, where S is\n",
           		gamma_insertion,gamma_insertion,
		           *pos,*(pos+1),*(pos+2),*(pos+3),sst1.re,sst1.im);
    message("the extended prop with insertion gamma_%d, should agree with\n",
		           gamma_insertion);
    message("<j_%d j_%d^+> on t_snk when using a point source for the first leg.\n",
		           gamma_insertion,gamma_insertion);
    message("Further test:\n");
    message("self-contraction of the extended prop: Tr(tilde S * tilde S^+):\n");

    for(t=0;t<L0*NPROC0;t++){
      message("t=%d %.6e + i %.6e \n",t,(*(sst2+t)).re,(*(sst2+t)).im); 
    }
    message("\n");
  }

  afree(st2);
  afree(sst2);

}


extern void C2_b(
	full_spinor_dble *svk, full_spinor_dble *svl, full_spinor_dble *svm, 
	int numb,int *mub, int nmom_max, int *momenta, complex_dble **s_pointer, int pol)
{
   int my_rank, max;
   double wt3,wt4;

   complex_dble *st21=NULL; 
   complex_dble *st22=NULL; 
   complex_dble *st23=NULL; 
   complex_dble *sst21=NULL; 
   complex_dble *sst22=NULL; 
   complex_dble *sst23=NULL; 

   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   wt3=MPI_Wtime();

   st21        = *(s_pointer + 0);
   st22        = *(s_pointer + 1);
   st23        = *(s_pointer + 2);
   sst21       = *(s_pointer + 3);
   sst22       = *(s_pointer + 4);
   sst23       = *(s_pointer + 5);

   max = NPROC0*L0*numb*nmom_max;

   flush_corr2(max,  st21);
   flush_corr2(max,  st22);
   flush_corr2(max,  st23);
   flush_corr2(max, sst21);
   flush_corr2(max, sst22);
   flush_corr2(max, sst23);

   /**** different polarisations ****/
   if(pol==0){
     twopt_baryons_dble(1,1,numb,mub,nmom_max,momenta,svk,svl,0.5,st21,st22,st23);
     twopt_baryons_dble(1,3,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
     twopt_baryons_dble(3,1,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
     twopt_baryons_dble(3,3,numb,mub,nmom_max,momenta,svk,svl,0.5,st21,st22,st23);
   }
   else if (pol==1){
     twopt_baryons_dble(1,1,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
     twopt_baryons_dble(1,3,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
     twopt_baryons_dble(3,1,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
     twopt_baryons_dble(3,3,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
   }
   else if(pol==2){
     twopt_baryons_dble(0,0,numb,mub,nmom_max,momenta,svk,svl,0.5,st21,st22,st23);
     twopt_baryons_dble(0,2,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
     twopt_baryons_dble(2,0,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
     twopt_baryons_dble(2,2,numb,mub,nmom_max,momenta,svk,svl,0.5,st21,st22,st23);
   }
   else if (pol==3){
     twopt_baryons_dble(0,0,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
     twopt_baryons_dble(0,2,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
     twopt_baryons_dble(2,0,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
     twopt_baryons_dble(2,2,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
   }
   else if (pol==4){
     twopt_baryons_dble(0,0,numb,mub,nmom_max,momenta,svk,svl,0.5,st21,st22,st23);
     twopt_baryons_dble(1,1,numb,mub,nmom_max,momenta,svk,svl,0.5,st21,st22,st23);
     twopt_baryons_dble(2,2,numb,mub,nmom_max,momenta,svk,svl,0.5,st21,st22,st23);
     twopt_baryons_dble(3,3,numb,mub,nmom_max,momenta,svk,svl,0.5,st21,st22,st23);
     twopt_baryons_dble(0,2,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
     twopt_baryons_dble(1,3,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
     twopt_baryons_dble(2,0,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
     twopt_baryons_dble(3,1,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
   }
   else if (pol==5){
     twopt_baryons_dble(0,0,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
     twopt_baryons_dble(1,1,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
     twopt_baryons_dble(2,2,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
     twopt_baryons_dble(3,3,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
     twopt_baryons_dble(0,2,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
     twopt_baryons_dble(1,3,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
     twopt_baryons_dble(2,0,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
     twopt_baryons_dble(3,1,numb,mub,nmom_max,momenta,svk,svl,-0.5,st21,st22,st23);
   }

   /**** combine the results from all processes in process 0 ****/
   MPI_Reduce(st21,sst21,2*max,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   MPI_Reduce(st22,sst22,2*max,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   MPI_Reduce(st23,sst23,2*max,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

   wt4=MPI_Wtime();
   message("\tTime for baryon 2pt contractions %.2e sec\n",wt4-wt3);
   fflush(stdout);

}


static void conserved(full_spinor_dble *sv0, full_spinor_dble *sv1,
                      int mu, int ix, complex_dble *corr)
{

  int y, ip1, ip2, z;
  complex_dble term1;
  full_spinor_dble shelp1, shelp2, shelp_send;
  su3_dble u;
  MPI_Status status;

  y = iup[ix][mu];

  if( y < VOLUME ){
    /**** Term 1 ****/
    shelp1 = su3_times_fsd(*(pud[ix][mu]),*(sv1+y));
    shelp2 = mul_one_pm_gamma_l(1.0, mu, shelp1);
    adj_full_spinor(sv0+ix,&shelp1);
    meson_trace(&shelp1,&shelp2,&term1);

    /**** Term 2 ****/
    _su3_dagger(u,*(pud[ix][mu]));
    shelp1 = su3_times_fsd(u,*(sv1+ix));
    shelp2 = mul_one_pm_gamma_l(-1.0, mu, shelp1);
    adj_full_spinor(sv0+y,&shelp1);
    meson_trace(&shelp1,&shelp2,corr);
  }
  /*************************/
  /**** boundary points ****/
  else{
    ip1 = npr[2 * mu];
    ip2 = npr[2 * mu + 1];
      z = map[y - VOLUME];

    /**** Term 1 ****/
    MPI_Sendrecv((double*) (sv1+z)    , 288, MPI_DOUBLE, ip1, 37,
                 (double*) &shelp_send, 288, MPI_DOUBLE, ip2, 37,
                                        MPI_COMM_WORLD, &status);

    shelp1 = su3_times_fsd(*(pud[ix][mu]),shelp_send);
    shelp2 = mul_one_pm_gamma_l(1.0, mu, shelp1);
    adj_full_spinor(sv0+ix,&shelp1);
    meson_trace(&shelp1,&shelp2,&term1);

    /**** Term 2 ****/
    MPI_Sendrecv((double*) (sv0+z)    , 288, MPI_DOUBLE, ip1, 37,
                 (double*) &shelp_send, 288, MPI_DOUBLE, ip2, 37,
                                        MPI_COMM_WORLD, &status);

    _su3_dagger(u,*(pud[ix][mu]));
    shelp1 = su3_times_fsd(u,*(sv1+ix));
    shelp2 = mul_one_pm_gamma_l(-1.0, mu, shelp1);
    adj_full_spinor(&shelp_send,&shelp1);
    meson_trace(&shelp1,&shelp2,corr);

  }

  /**** Term1 + Term2 ***/
  (*corr).re = 0.5*((*corr).re + term1.re); 
  (*corr).im = 0.5*((*corr).im + term1.im); 
  
}


static void derivative(full_spinor_dble *sv0, full_spinor_dble *sv1,
                       int mu, int ix, complex_dble *corr)
{

  int ix_dn, ix_up, ip1, ip2, z, nu;
  complex_dble term1, term2, term3, term4;
  full_spinor_dble fsd_up[2], fsd_dn[2], fsd_send[2];
  full_spinor_dble shelp, sv0_ix, shelp_term1, shelp_term2;
  su3_dble u_help;
  MPI_Status status; 

  /**** communication part ****/
  ix_dn = idn[ix][mu];
  ix_up = iup[ix][mu];
  ip1 = npr[2 * mu];
  ip2 = npr[2 * mu + 1];

  if(ix_up < VOLUME){  
    adj_full_spinor(sv0+ix_up, &(fsd_up[0]) );
    fsd_up[1] = *(sv1+ix_up);
  }
  else{
    z = map[ix_up - VOLUME];
    adj_full_spinor(sv0+z, &(fsd_send[0]) );
    fsd_send[1] = *(sv1+z);
    MPI_Sendrecv((double*) fsd_send, 576, MPI_DOUBLE, ip1, 37,
                 (double*) fsd_up  , 576, MPI_DOUBLE, ip2, 37,
                                     MPI_COMM_WORLD, &status);
  }
  if(ix_dn < VOLUME){
    _su3_dagger(u_help, *(pud[ix_dn][mu]));
    adj_full_spinor(sv0+ix_dn, &shelp);
    fsd_dn[0] = fsd_times_su3(*(pud[ix_dn][mu]),shelp);
    fsd_dn[1] = su3_times_fsd(u_help,*(sv1+ix_dn));
  }
  else{
    z = map[ix_dn - VOLUME];
    _su3_dagger(u_help, *(pud[z][mu]));
    adj_full_spinor(sv0+z, &shelp);
    fsd_send[0] = fsd_times_su3(*(pud[z][mu]),shelp);
    fsd_send[1] = su3_times_fsd(u_help,*(sv1+z));
    MPI_Sendrecv((double*) fsd_send, 576, MPI_DOUBLE, ip2, 37,
                 (double*) fsd_dn  , 576, MPI_DOUBLE, ip1, 37,
                                     MPI_COMM_WORLD, &status);
  }

  /**** some other preparations ****/
  adj_full_spinor(sv0+ix, &(sv0_ix));
  shelp_term1 = su3_times_fsd(*(pud[ix][mu]), fsd_up[1]);
  _su3_dagger(u_help, *(pud[ix][mu]));
  shelp_term2 = su3_times_fsd(u_help, *(sv1+ix));

  /**** combining everything ****/
  for(nu=0;nu<16;nu++){
    mul_gamma_l(nu, &shelp_term1, &shelp);
    meson_trace(&(sv0_ix), &shelp, &term1);

    mul_gamma_l(nu, &(fsd_dn[1]), &shelp);
    meson_trace(&(sv0_ix), &shelp, &term2);

    mul_gamma_l(nu, &shelp_term2, &shelp);
    meson_trace(&(fsd_up[0]), &shelp, &term3);

    mul_gamma_l(nu, sv1+ix, &shelp);
    meson_trace(&(fsd_dn[0]), &shelp, &term4);

    (*(corr+nu)).re = 0.25*(term1.re - term2.re + term3.re - term4.re);
    (*(corr+nu)).im = 0.25*(term1.im - term2.im + term3.im - term4.im);
  }
}


static void derivative2(full_spinor_dble *sv0, full_spinor_dble *sv1,
                       int mu, int ix, complex_dble *corr)
{

  int ix_dn, ix_up, ip1, ip2, z, nu;
  complex_dble term1, term2, term3, term4;
  full_spinor_dble fsd_up[2], fsd_dn[2], fsd_send[2];
  full_spinor_dble shelp, sv1_ix, shelp_term1, shelp_term2;
  su3_dble u_help;
  MPI_Status status; 

  /**** communication part ****/
  ix_dn = idn[ix][mu];
  ix_up = iup[ix][mu];
  ip1 = npr[2 * mu];
  ip2 = npr[2 * mu + 1];

  if(ix_up < VOLUME){  
    adj_full_spinor(sv1+ix_up, &(fsd_up[0]) );
    fsd_up[1] = *(sv0+ix_up);
  }
  else{
    z = map[ix_up - VOLUME];
    adj_full_spinor(sv1+z, &(fsd_send[0]) );
    fsd_send[1] = *(sv0+z);
    MPI_Sendrecv((double*) fsd_send, 576, MPI_DOUBLE, ip1, 37,
                 (double*) fsd_up  , 576, MPI_DOUBLE, ip2, 37,
                                     MPI_COMM_WORLD, &status);
  }
  if(ix_dn < VOLUME){
    _su3_dagger(u_help, *(pud[ix_dn][mu]));
    adj_full_spinor(sv1+ix_dn, &shelp);
    fsd_dn[0] = fsd_times_su3(*(pud[ix_dn][mu]),shelp);
    fsd_dn[1] = su3_times_fsd(u_help,*(sv0+ix_dn));
  }
  else{
    z = map[ix_dn - VOLUME];
    _su3_dagger(u_help, *(pud[z][mu]));
    adj_full_spinor(sv1+z, &shelp);
    fsd_send[0] = fsd_times_su3(*(pud[z][mu]),shelp);
    fsd_send[1] = su3_times_fsd(u_help,*(sv0+z));
    MPI_Sendrecv((double*) fsd_send, 576, MPI_DOUBLE, ip2, 37,
                 (double*) fsd_dn  , 576, MPI_DOUBLE, ip1, 37,
                                     MPI_COMM_WORLD, &status);
  }

  /**** some other preparations ****/
  adj_full_spinor(sv1+ix, &(sv1_ix));
  shelp_term1 = su3_times_fsd(*(pud[ix][mu]), fsd_up[1]);
  _su3_dagger(u_help, *(pud[ix][mu]));
  shelp_term2 = su3_times_fsd(u_help, *(sv0+ix));

  /**** combining everything ****/
  for(nu=0;nu<16;nu++){
    mul_gamma_l(nu, &shelp_term1, &shelp);
    meson_trace(&(sv1_ix), &shelp, &term1);

    mul_gamma_l(nu, &(fsd_dn[1]), &shelp);
    meson_trace(&(sv1_ix), &shelp, &term2);

    mul_gamma_l(nu, &shelp_term2, &shelp);
    meson_trace(&(fsd_up[0]), &shelp, &term3);

    mul_gamma_l(nu, sv1+ix, &shelp);
    meson_trace(&(fsd_dn[0]), &shelp, &term4);

    (*(corr+nu)).re = 0.25*(term1.re - term2.re + term3.re - term4.re);
    (*(corr+nu)).im = 0.25*(term1.im - term2.im + term3.im - term4.im);
  }
}
static void fourier_trafo(int *i, int nmom_max, double *phase, 
                          complex_dble *in, complex_dble *out)
{
  int ip;
  double phase0, phase_sin, phase_cos;
  complex_dble *xhelp2, corr_loc;

  corr_loc=*(in);  
  for(ip=0;ip<nmom_max;ip++){
    phase0 = phase[ip];
    phase_sin = sin(phase0);
    phase_cos = cos(phase0);
    xhelp2 = (out + *i);     
    (*xhelp2).re += corr_loc.re*phase_cos - corr_loc.im*phase_sin;
    (*xhelp2).im += corr_loc.im*phase_cos + corr_loc.re*phase_sin; 
    *i=*i+1;
  }

}

extern void C3_b(full_spinor_dble *sv0, full_spinor_dble *sv1, 
 	               int numb, int *mub, int nmom_max, int *momenta, 
 	               complex_dble *s_pointer, propinfo prop) 
{
  int i,my_rank,imu,mu,max,max2,x0,x1,x2,x3,ip,ix,T,pt,nu,j;
  double wt1,wt2,*phase,phase0,phase_sin,phase_cos;
  double X,Y,Z,px,py,pz,momx,momy,momz;
  double PI;
  int num_mom=nmom_max;
 
  complex_dble *st=NULL,*xt=NULL,*xhelp,*xhelp1,*xhelp2,corr[16],corr_loc;
  full_spinor_dble shelp1,shelp2;

  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  wt1 = MPI_Wtime();
   
  PI=0.5*6.283185307179586476925286;

  max   = 16*numb*num_mom*NPROC0*L0;
  max2  = 16*numb*num_mom;

  st    = amalloc(max*sizeof(complex_dble),3);
  xt    = amalloc(max2*L0*sizeof(complex_dble),3);
  phase = amalloc(num_mom*sizeof(double),3);

  flush_corr2(max, s_pointer);
  flush_corr2(max2*L0, xt);
  flush_corr2(max, st); 

  /**** shortcuts to save compute-time ****/
  T    =          (NPROC0*L0);
  X    = (double) (NPROC1*L1);
  Y    = (double) (NPROC2*L2);
  Z    = (double) (NPROC3*L3);
  pt   =          cpr[0]*L0;
  px   = (double) cpr[1]*L1;
  py   = (double) cpr[2]*L2;
  pz   = (double) cpr[3]*L3;
  momx = (double) *(prop.mom_ins+0);
  momy = (double) *(prop.mom_ins+1);
  momz = (double) *(prop.mom_ins+2);
 
  for(x1=0;x1<L1;x1++){
    for(x2=0;x2<L2;x2++){
      for(x3=0;x3<L3;x3++){
      
        /**** calculating the phase for fourier transformation ****/
        for(ip=0;ip<nmom_max;ip++){
          phase[ip] = 2.0*PI*(
             (double)( momx - (*(momenta+ip*3  )))*(px+x1)/X +
             (double)( momy - (*(momenta+ip*3+1)))*(py+x2)/Y +
             (double)( momz - (*(momenta+ip*3+2)))*(pz+x3)/Z);
        }
        
        for(x0=0;x0<L0;x0++){
      
          ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];
          i=0;

          xhelp = (xt + x0*max2);         

          for(imu=0;imu<numb;imu++){
            
            /****   computing the operator insertion   ****/
            mu = mub[imu];

            /**** conserved currents ****/
            if(mu==0){
              for(nu=0;nu<4;nu++){
                conserved(sv0, sv1, nu, ix, corr+nu);
                fourier_trafo(&i, nmom_max, phase, corr+nu, xhelp);
              }
            }
            /**** local currents ****/
            else if(mu==1){
              adj_full_spinor(sv0+ix,&shelp2);
              for(nu=0;nu<16;nu++){
                mul_gamma_l(nu,sv1+ix,&shelp1);
                meson_trace(&shelp1,&shelp2,&corr_loc);
                fourier_trafo(&i, nmom_max, phase, &corr_loc, xhelp);
              }
            }
            /**** derivative operators ****/
            else{
              derivative2(sv0, sv1, mu-2, ix, corr);
              for(j=0;j<16;j++){
                fourier_trafo(&i, nmom_max, phase, corr+j, xhelp);
              }
            }

          }

        }
      }
    }
  }


  for(x0=0;x0<L0;x0++){
    xhelp  = (xt + x0*max2);
    xhelp1 = (st + pt + x0);
    for(i=0;i<max2;i++){
      (*(xhelp1 + i*T)).re = (*(xhelp + i)).re;
      (*(xhelp1 + i*T)).im = (*(xhelp + i)).im;
    }
  }
  
  /**** combine the results from all processes in process 0 ****/
  MPI_Reduce(st,s_pointer,2*max,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
 
  afree(st);
  afree(phase); 
  afree(xt);

  wt2 = MPI_Wtime();
  message("\tTime for baryon 3pt contraction: %.2e sec\n",wt2-wt1);
  
}

















