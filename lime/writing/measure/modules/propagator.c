/*******************************************************************************
*
* File propagator.c
*
* Copyright (C) 2008-11 Andreas Juettner
*                       Bastian Knippschild
*
* Externally available functions:
*
* extern void propagator(full_spinor_dble *sv, spinor_dble *rnd, propinfo prop)
* 	Computes 12-colum (full_spinor_dble) propagator for the 
*	propagator description provided in prop (including source type).
*	the field rnd is a spinor_dble and contains a random gauge
*  	field of the desired type which will be used as the basis
*	for the construction of the noise source
*
* extern void seq_propagator(full_spinor_dble *sv, full_spinor_dble *src, 
*		propinfo prop)
*	Given the full_spinor_dble *sv this function constructs the
*	sequential source specified further in prop
*
* extern void propagator_col(spinor_dble *sv, spinor_dble *rnd, propinfo prop,
*		int idirac, int icolour)
* 	this function works in the same way as propagator(...) but returns
*	only the propagator column specified by idirca and icolour
*
* extern void seq_propagator_b(full_spinor_dble *sv_out, full_spinor_dble *svk, full_spinor_dble *svl, 
*                              full_spinor_dble *svm, spinor_dble **wsd, spinor_dble *grand_source, 
*                              propinfo prop, int iw0, int iw3)
* 
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*
*
*******************************************************************************/
#define CORRELATORS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#include "random.h"
#include "misc.h"
#include "start.h"
#include "linalg.h"
#include "sw_term.h"
#include "dirac.h"
#include "dfl.h"
#include "sap_gcr.h"
#include "version.h"
#include "version_qcd-measure.h"
#include "global.h"
#include "update.h"
#include "measure.h"

static double residue(spinor_dble *s,spinor_dble *r)
{
   int ir,iw;
   double d;
   complex_dble z;
   spinor_dble **wsd;
   wsd=reserve_wsd(1);

   ir=(r-psd[0][0])/NSPIN;   
   iw=(wsd[0]-psd[0][0])/NSPIN;

   Qnohat_dble(ir,iw);
   mulg5_dble(VOLUME,wsd[0]);
   z.re=-1.0;
   z.im=0.0;
   mulc_spinor_add_dble(VOLUME,wsd[0],s,z);
   d=norm_square_dble(VOLUME,1,wsd[0])/norm_square_dble(VOLUME,1,s);

   release_wsd();

   return sqrt(d);
}



extern void propagator(full_spinor_dble *sv, spinor_dble *rnd, propinfo prop)
{
   int idirac, icolour, status[3], ifail, shift[4];
   int iw0,iw1,iw2;
   double wt1,wt2;
   double d1,d2,cres;
   lat_parms_t lat; 
   spinor_dble **wsd;

   wsd=reserve_wsd(3);
   iw0=(wsd[0]-psd[0][0])/NSPIN;
   iw1=(wsd[1]-psd[0][0])/NSPIN;
   iw2=(wsd[2]-psd[0][0])/NSPIN;
   alloc_pud_sm1();

   lat=set_lat_parms(prop.beta,prop.kappa,prop.csw);

   wt1=MPI_Wtime();  

   twistbc(&(prop.theta[0]),+1.);   
   message("Applied twist angle (%f,%f,%f,%f)\n",
      	prop.theta[0],prop.theta[1],prop.theta[2],prop.theta[3]);
/*
   shift[0] = -prop.src.pos[0];
   shift[1] = -prop.src.pos[1];
   shift[2] = -prop.src.pos[2];
   shift[3] = -prop.src.pos[3];
   shift_ud(shift);
   sw_term();
   copy_bnd_ud();
   message("Shifted gauge field by (%d,%d,%d,%d)\n",
      	shift[0], shift[1], shift[2], shift[3]);
*/
   if(strcmp(prop.src.type, "J_thin") == 0){
     pud_copy(pud, pud_sm1);
     copy_bnd_ud_reverse(pud_sm1);
     message("gaussian smearing with thin links\n");
   }
   else if (strcmp(prop.src.type, "J_APE") == 0){
     APE_smearing(pud, pud_sm1, prop.src.smear);
     copy_bnd_ud_reverse(pud_sm1);
     message("gaussian smearing with APE links\n");
   }
   else if (strcmp(prop.src.type, "J_HYP") == 0){
     HYP_smearing(pud, pud_sm1, prop.src.smear);
     copy_bnd_ud_reverse(pud_sm1);
     message("gaussian smearing with HYP links\n");fflush(stdout);
   }
/*
   dfl_modes(0,status,&ifail);
   message("Deflation subspace generation: status = %d, ifail = %d\n\n",
           status[0],ifail);
*/
   message("Inverting propagator with kappa=%lf, precision %e\n",prop.kappa,
      		prop.res);
   fflush(stdout);

   for (idirac=1;idirac<=4;idirac++)
   {
    for (icolour=1;icolour<=3;icolour++)
    {
     prop.src.pos[0]=0;
     prop.src.pos[1]=0;
     prop.src.pos[2]=0;
     prop.src.pos[3]=0;
     srcfld(rnd,wsd[0],prop.src,idirac,icolour);         
     dfl_sap_gcr(prop.nmx,prop.res,iw0,iw2,status,&ifail);
     cres=residue(wsd[0],wsd[2]);
     d1=norm_square_dble(VOLUME,1,wsd[0]);
     d2=norm_square_dble(VOLUME,1,wsd[2]);
     copy_sv_fs(VOLUME,wsd[2],sv,idirac,icolour);
     fflush(stdout);
     message("(%d%d)-component: status = %3d (little: %3d,%3d), "
           "ifail = %d, check = %.1e\n",idirac,icolour,status[0],
           status[1],status[2],ifail,cres);
     error_root((status[0]<0)||(ifail!=0),1,"propagator [propagator.c]",
           "Propagator computation failed");
     fflush(stdout);
    }
   }
   wt2=MPI_Wtime();
   message("Time for inversion %.2e sec\n",wt2-wt1);
     fflush(stdout);
   release_wsd();
/*
   shift[0] = -shift[0];
   shift[1] = -shift[1];
   shift[2] = -shift[2];
   shift[3] = -shift[3];
   shift_ud(shift);
*/   twistbc(&(prop.theta[0]),-1.);
   fflush(stdout);
   free_pud_sm1();
}



extern void seq_propagator(full_spinor_dble *sv, full_spinor_dble *src, 
		propinfo prop)
{
   int idirac, icolour, status[3], ifail, shift[4];
   int iw0, iw1;
   double d1,d2,cres,wt1,wt2;
   lat_parms_t lat; 
   spinor_dble **wsd;
   full_spinor_dble *svaux;
   
   wsd=reserve_wsd(2);
   iw0=(wsd[0]-psd[0][0])/NSPIN;
   iw1=(wsd[1]-psd[0][0])/NSPIN;
   svaux   = amalloc(VOLUME*sizeof(full_spinor_dble),3); 

   set_fsd2zero(VOLUME,svaux);
   message("Creating ext. src with gamma_snk=%d, kappa=%.6f and t_snk=%d\n",
   	prop.gam_ins, prop.kappa,prop.t_snk);
   error_root((prop.t_snk>NPROC0*L0),0,"main [measure6.c]",
     "t_snk=%d is larger than NPROC0*L0=%d",prop.t_snk,NPROC0*L0);
   extended_src(prop.gam_ins,prop.mom_ins,prop.t_snk,src,svaux);   
   message("mom_insert p=(%d %d %d)/2pi\n",
     *(prop.mom_ins), *(prop.mom_ins+1), *(prop.mom_ins+2));
   twistbc(&(prop.theta[0]),+1.);
   lat=set_lat_parms(prop.beta,prop.kappa,prop.csw);
   message("Applied twist angle (%f,%f,%f)\n",
     prop.theta[1],prop.theta[2],prop.theta[3]);
   wt1=MPI_Wtime();
   shift[0] = -prop.src.pos[0];
   shift[1] = -prop.src.pos[1];
   shift[2] = -prop.src.pos[2];
   shift[3] = -prop.src.pos[3];
   shift_ud(shift);
   message("Shifted gauge field by (%d,%d,%d,%d)\n",
      	shift[0], shift[1], shift[2], shift[3]);
   sw_term();
   copy_bnd_ud();

   for (idirac=1;idirac<=4;idirac++)
   {
    for (icolour=1;icolour<=3;icolour++)
    {
     set_sd2zero(VOLUME,wsd[0]);
     set_sd2zero(VOLUME,wsd[1]);
     copy_fs_sv(VOLUME,svaux,wsd[0],idirac,icolour);
     d1=norm_square_dble(VOLUME,1,wsd[0]);
     sap_gcr(prop.nmx,prop.res,iw0,iw1,status,&ifail); /* fix dfl_sap_gcr!!! */
     copy_sv_fs(VOLUME,wsd[1],sv,idirac,icolour);
     cres=residue(wsd[0],wsd[1]);
     d2=norm_square_dble(VOLUME,1,wsd[1]);
     message("(%d%d)-component: status = %3d , "
         "ifail = %d, check = %.1e\n",idirac,icolour,status[0],ifail,cres);
     error_root((status[0]<0)||(ifail!=0),1,"main [measure6.c]",
          "Propagator computation failed with status[0] %d and ifail %d",
      	  status[0],ifail);
     fflush(stdout);
    }
   }
   shift[0] = +prop.src.pos[0];
   shift[1] = +prop.src.pos[1];
   shift[2] = +prop.src.pos[2];
   shift[3] = +prop.src.pos[3];
   shift_ud(shift);
   wt2=MPI_Wtime();
   message("Time for inversion %.2e sec\n",wt2-wt1);
   twistbc(&(prop.theta[0]),-1.);
   afree(svaux);
   release_wsd();
}



extern void seq_propagator_b(full_spinor_dble *sv_out,full_spinor_dble *svk,
			     full_spinor_dble *svl,   full_spinor_dble *svm,
			     propinfo prop, int snk_sm, int pol)
{

   int icolour, idirac, ifail, shift[4],status[3],time[4];
   int my_rank, ix, i, x0, x1, x2, x3, iw0, iw1;
   double cres, wt3,wt4, phase;
   const double PI=0.5*6.283185307179586476925286;
   lat_parms_t lat;
   complex_dble *s=NULL, *ss=NULL;
   spinor_dble shelp, *grand_source,**wsd;
   static const spinor_dble sd0={{{0.0}}};

   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

   lat=set_lat_parms(prop.beta,prop.kappa,prop.csw);

   alloc_pud_sm1();
   wsd=reserve_wsd(2);
   iw0=(wsd[0]-psd[0][0])/NSPIN;
   iw1=(wsd[1]-psd[0][0])/NSPIN;

   wt3=MPI_Wtime();

   twistbc(&prop.theta[0],+1.);
   message("Applied twist angle (%f,%f,%f,%f)\n",
      	prop.theta[0],prop.theta[1],prop.theta[2],prop.theta[3]);
/*
   shift[0] = -prop.src.pos[0];
   shift[1] = -prop.src.pos[1];
   shift[2] = -prop.src.pos[2];
   shift[3] = -prop.src.pos[3];
   shift_ud(shift);
   sw_term();
   copy_bnd_ud();   
   message("gauge field shifted by (%d, %d, %d, %d)\n",
           -prop.src.pos[0], -prop.src.pos[1],
           -prop.src.pos[2], -prop.src.pos[3] );
*/
   if(strcmp(prop.src.type, "J_thin") == 0){
     pud_copy(pud, pud_sm1);
     copy_bnd_ud_reverse(pud_sm1);
     message("gaussian smearing with thin links\n");
   }
   else if (strcmp(prop.src.type, "J_APE") == 0){
     APE_smearing(pud, pud_sm1, prop.src.smear);
     copy_bnd_ud_reverse(pud_sm1);
     message("gaussian smearing with APE links\n");
   }
   else if (strcmp(prop.src.type, "J_HYP") == 0){
     HYP_smearing(pud, pud_sm1, prop.src.smear);
     copy_bnd_ud_reverse(pud_sm1);
     message("gaussian smearing with HYP links\n");
   }
   fflush(stdout);

   grand_source = amalloc(12*VOLUME*sizeof(spinor_dble),3);	
   for (i=0;i<(12*VOLUME);i++){
     *(grand_source+i) = sd0;
   }

   if(pol==0){
     message("\t polarisation matrix no 0\n");
     d_source_dble(prop.t_snk,1,1,svk,0.5,grand_source);
     d_source_dble(prop.t_snk,1,3,svk,-0.5,grand_source);
     d_source_dble(prop.t_snk,3,1,svk,-0.5,grand_source);
     d_source_dble(prop.t_snk,3,3,svk,0.5,grand_source);

     u_source_dble(prop.t_snk,1,1,svk,svl,-0.5,grand_source);
     u_source_dble(prop.t_snk,1,3,svk,svl,0.5,grand_source);
     u_source_dble(prop.t_snk,3,1,svk,svl,0.5,grand_source);
     u_source_dble(prop.t_snk,3,3,svk,svl,-0.5,grand_source);
   }
   else if(pol==1){
     message("\t polarisation matrix no 1\n");
     d_source_dble(prop.t_snk,1,1,svk,-0.5,grand_source);
     d_source_dble(prop.t_snk,1,3,svk,-0.5,grand_source);
     d_source_dble(prop.t_snk,3,1,svk,-0.5,grand_source);
     d_source_dble(prop.t_snk,3,3,svk,-0.5,grand_source);

     u_source_dble(prop.t_snk,1,1,svk,svl,0.5,grand_source);
     u_source_dble(prop.t_snk,1,3,svk,svl,0.5,grand_source);
     u_source_dble(prop.t_snk,3,1,svk,svl,0.5,grand_source);
     u_source_dble(prop.t_snk,3,3,svk,svl,0.5,grand_source);
   }
   else if(pol==2){
     message("\t polarisation matrix no 2\n");
     d_source_dble(prop.t_snk,0,0,svk,0.5,grand_source);
     d_source_dble(prop.t_snk,0,2,svk,-0.5,grand_source);
     d_source_dble(prop.t_snk,2,0,svk,-0.5,grand_source);
     d_source_dble(prop.t_snk,2,2,svk,0.5,grand_source);

     u_source_dble(prop.t_snk,0,0,svk,svl,-0.5,grand_source);
     u_source_dble(prop.t_snk,0,2,svk,svl,0.5,grand_source);
     u_source_dble(prop.t_snk,2,0,svk,svl,0.5,grand_source);
     u_source_dble(prop.t_snk,2,2,svk,svl,-0.5,grand_source);
   }
   else if(pol==3){
     message("\t polarisation matrix no 3\n");
     d_source_dble(prop.t_snk,0,0,svk,-0.5,grand_source);
     d_source_dble(prop.t_snk,0,2,svk,-0.5,grand_source);
     d_source_dble(prop.t_snk,2,0,svk,-0.5,grand_source);
     d_source_dble(prop.t_snk,2,2,svk,-0.5,grand_source);

     u_source_dble(prop.t_snk,0,0,svk,svl,0.5,grand_source);
     u_source_dble(prop.t_snk,0,2,svk,svl,0.5,grand_source);
     u_source_dble(prop.t_snk,2,0,svk,svl,0.5,grand_source);
     u_source_dble(prop.t_snk,2,2,svk,svl,0.5,grand_source);
   }
   else if(pol==4){
     message("\t unpolarised nucleon\n");
     d_source_dble(prop.t_snk,0,0,svk,0.5,grand_source);
     d_source_dble(prop.t_snk,1,1,svk,0.5,grand_source);
     d_source_dble(prop.t_snk,2,2,svk,0.5,grand_source);
     d_source_dble(prop.t_snk,3,3,svk,0.5,grand_source);
     d_source_dble(prop.t_snk,0,2,svk,-0.5,grand_source);
     d_source_dble(prop.t_snk,1,3,svk,-0.5,grand_source);
     d_source_dble(prop.t_snk,2,0,svk,-0.5,grand_source);
     d_source_dble(prop.t_snk,3,1,svk,-0.5,grand_source);

     u_source_dble(prop.t_snk,0,0,svk,svl,-0.5,grand_source);
     u_source_dble(prop.t_snk,1,1,svk,svl,-0.5,grand_source);
     u_source_dble(prop.t_snk,2,2,svk,svl,-0.5,grand_source);
     u_source_dble(prop.t_snk,3,3,svk,svl,-0.5,grand_source);
     u_source_dble(prop.t_snk,0,2,svk,svl,0.5,grand_source);
     u_source_dble(prop.t_snk,1,3,svk,svl,0.5,grand_source);
     u_source_dble(prop.t_snk,2,0,svk,svl,0.5,grand_source);
     u_source_dble(prop.t_snk,3,1,svk,svl,0.5,grand_source);
   }
   else if(pol==5){
     message("\t unpolarised anti-nucleon\n");
     d_source_dble(prop.t_snk,0,0,svk,-0.5,grand_source);
     d_source_dble(prop.t_snk,1,1,svk,-0.5,grand_source);
     d_source_dble(prop.t_snk,2,2,svk,-0.5,grand_source);
     d_source_dble(prop.t_snk,3,3,svk,-0.5,grand_source);
     d_source_dble(prop.t_snk,0,2,svk,-0.5,grand_source);
     d_source_dble(prop.t_snk,1,3,svk,-0.5,grand_source);
     d_source_dble(prop.t_snk,2,0,svk,-0.5,grand_source);
     d_source_dble(prop.t_snk,3,1,svk,-0.5,grand_source);

     u_source_dble(prop.t_snk,0,0,svk,svl,0.5,grand_source);
     u_source_dble(prop.t_snk,1,1,svk,svl,0.5,grand_source);
     u_source_dble(prop.t_snk,2,2,svk,svl,0.5,grand_source);
     u_source_dble(prop.t_snk,3,3,svk,svl,0.5,grand_source);
     u_source_dble(prop.t_snk,0,2,svk,svl,0.5,grand_source);
     u_source_dble(prop.t_snk,1,3,svk,svl,0.5,grand_source);
     u_source_dble(prop.t_snk,2,0,svk,svl,0.5,grand_source);
     u_source_dble(prop.t_snk,3,1,svk,svl,0.5,grand_source);
   }

/*
   dfl_modes(0,status,&ifail);
   message("Deflation subspace generation: status = %d, ifail = %d\n\n",
           status[0],ifail);
*/
   message("Inverting propagator with kappa=%lf, precision %e\n",prop.kappa,
	   prop.res);
		
   for (idirac=1;idirac<=4;idirac++)
   {
    for (icolour=1;icolour<=3;icolour++)
    {
     set_sd2zero(VOLUME,wsd[0]);
     set_sd2zero(VOLUME,wsd[1]);
     for (x0=0;x0<L0;x0++)
     {
      for (x1=0;x1<L1;x1++)
      {
       for (x2=0;x2<L2;x2++)
       {
        for (x3=0;x3<L3;x3++)
        {
         ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];
         if (coords[ix].t==prop.t_snk)
         {
          s  = (complex_dble*)(wsd[0]+ix);
          ss = (complex_dble*)(grand_source+ix
                                 +VOLUME*(icolour-1+3*(idirac-1)));

          /****  multiplication by  exp{-ipx}  (momentum at the sink) ****/
          phase = 2.0*PI*(
                (double)((*(prop.mom_ins))*(cpr[1]*L1+x1))/
                                (double)(NPROC1*L1) +
                (double)((*(prop.mom_ins+1))*(cpr[2]*L2+x2))/
                                (double)(NPROC2*L2) +
                (double)((*(prop.mom_ins+2))*(cpr[3]*L3+x3))/
                                (double)(NPROC3*L3));
 
          for(i=0;i<12;i++)
          {
           (*(s+i)).re  = (*(ss+i)).re*cos(phase) + (*(ss+i)).im*sin(phase);
           (*(s+i)).im  = (*(ss+i)).im*cos(phase) - (*(ss+i)).re*sin(phase);
          }
         }
        }
       }
      }
     }
 
     /* sinksmearing */
     if( snk_sm >= 0 ){
       for(i=0;i<4;i++){
         time[i] = prop.src.pos[i];
         prop.src.pos[i] = 0;
       }
       prop.src.pos[0] = prop.t_snk;
       jacobi_source(wsd[0],prop.src,idirac,icolour);
       for(i=0;i<4;i++)
         prop.src.pos[i] = time[i];
     }

     dfl_sap_gcr(prop.nmx,prop.res,iw0,iw1,status,&ifail);
     cres=residue(wsd[0],wsd[1]);
  
     message("(%d%d)-component: status = %3d (little: %3d,%3d), "
             "ifail = %d, check = %.1e\n",idirac,icolour,status[0],
             status[1],status[2],ifail,cres);
     error_root((status[0]<0)||(ifail!=0),1,"seq_propagator_b [propagators.c]",
                 "Propagator computation failed");
     fflush(stdout);
 
     for (ix=0;ix<VOLUME;ix++){
       mul_s_gamma_l(5, (wsd[1]+ix), &shelp);
       *(wsd[1]+ix) = shelp;
     }
     copy_sv_fs(VOLUME,wsd[1],sv_out,idirac,icolour);
    }
   }

   /* check routine */
/*   check_closure_dble(grand_source, svk, sv_out, prop); 
*/
/*
   shift[0] = +prop.src.pos[0];
   shift[1] = +prop.src.pos[1];
   shift[2] = +prop.src.pos[2];
   shift[3] = +prop.src.pos[3];
   shift_ud(shift);
*/
   wt4=MPI_Wtime();
   message("Time for baryon extended source: %.2e sec\n", wt4-wt3);   

   twistbc(&prop.theta[0],-1.);
   afree(grand_source);
   release_wsd();
   free_pud_sm1();
}








/************ DON'T USE THIS ONE - UNDER DEVELOPMENT * AJ ********************/
extern void propagator_col(spinor_dble *sv, spinor_dble *rnd, propinfo prop,
		int idirac, int icolour)
{
   int i,status[3], ifail;
   int iw0,iw1,iw2;
   double wt1,wt2;
   double d1,d2,cres;
   lat_parms_t lat; 
   spinor_dble **wsd;
   wsd=reserve_wsd(3);
   iw0=(wsd[0]-psd[0][0])/NSPIN;
   iw1=(wsd[1]-psd[0][0])/NSPIN;
   iw2=(wsd[2]-psd[0][0])/NSPIN;
     lat=set_lat_parms(prop.beta,prop.kappa,prop.csw);
     twistbc(&(prop.theta[0]),+1.);
     sw_term();
     message("Applied twist angle (%f,%f,%f,%f)\n",
		prop.theta[0],prop.theta[1],prop.theta[2],prop.theta[3]);
     message("inverting propagator with kappa=%lf, precision %e\n",prop.kappa,
			prop.res);
     wt1=MPI_Wtime();  
/*     shift[0] = prop.src.pos[0];
     shift[1] = prop.src.pos[1];
     shift[2] = prop.src.pos[2];
     shift[3] = prop.src.pos[3];
     shift_ud(shift);
*/

       srcfld(rnd,wsd[0],prop.src,idirac,icolour);         
       dfl_modes(0,status,&ifail); /* Check whether they really have to be
				    generated every time */
       message("Deflation subspace generation: status = %d, ifail = %d\n\n",
             status[0],ifail);
       SourceRadius(prop.src,wsd[0]);
       for (i=0;i<VOLUME;i++){
        *(rnd+i) = *(wsd[0]+i);
       }
       dfl_sap_gcr(prop.nmx,prop.res,iw0,iw2,status,&ifail);
       cres=residue(wsd[0],wsd[2]);
       d1=norm_square_dble(VOLUME,1,wsd[0]);
       for (i=0;i<VOLUME;i++){
        *(sv+i)  = *(wsd[2]+i);
       }
       d2=norm_square_dble(VOLUME,1,sv);
       message("(%d%d)-component: status = %3d (little: %3d,%3d), "
             "ifail = %d, check = %.1e\n",idirac,icolour,status[0],
             status[1],status[2],ifail,cres);
       message("norm squared of source %.2e and of solution %.2e\n",d1,d2);
       error_root((status[0]<0)||(ifail!=0),1,"propagator [propagator.c]",
             "Propagator computation failed");
       fflush(stdout);

     twistbc(&(prop.theta[0]),-1.);
     wt2=MPI_Wtime();
     message("Time for inversion %.2e sec\n",wt2-wt1);
     release_wsd();

 /*    shift[0] = -prop.src.pos[0];
     shift[1] = -prop.src.pos[1];
     shift[2] = -prop.src.pos[2];
     shift[3] = -prop.src.pos[3];
     shift_ud(shift);
*/
}
