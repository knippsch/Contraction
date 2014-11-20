
/*******************************************************************************
*
* File twistbc.c
*
* Copyright (C) 2008 Andreas Juettner
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Programs to switch the boundary conditions in the time direction
*
*   void twistedbc(void)
*     multiplies the phase exp(i(\theta_1/L1+\theta_2/L2+\theta_3/L3))
*     to each single precision link variable on the lattice
*
*   void twistedbcd(void)
*     multiplies the phase exp(i(\theta_1/L1+\theta_2/L2+\theta_3/L3))
*     to each double precision link variable on the lattice
*
*******************************************************************************/

#define TWISTBC_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "start.h"
#include "global.h"
#include "measure.h"

static void su3_appl_twist(su3_dble *u,double phase){
   double c,s;
   su3_dble *r;
   r=amalloc(sizeof(su3_dble),3);
   c = (double)cos(phase);
   s = (double)sin(phase);
   (*(r)).c11.re=c*(*(u)).c11.re-s*(*(u)).c11.im; 
   (*(r)).c11.im=s*(*(u)).c11.re+c*(*(u)).c11.im; 
   (*(r)).c12.re=c*(*(u)).c12.re-s*(*(u)).c12.im; 
   (*(r)).c12.im=s*(*(u)).c12.re+c*(*(u)).c12.im; 
   (*(r)).c13.re=c*(*(u)).c13.re-s*(*(u)).c13.im; 
   (*(r)).c13.im=s*(*(u)).c13.re+c*(*(u)).c13.im; 
   (*(r)).c21.re=c*(*(u)).c21.re-s*(*(u)).c21.im; 
   (*(r)).c21.im=s*(*(u)).c21.re+c*(*(u)).c21.im; 
   (*(r)).c22.re=c*(*(u)).c22.re-s*(*(u)).c22.im; 
   (*(r)).c22.im=s*(*(u)).c22.re+c*(*(u)).c22.im; 
   (*(r)).c23.re=c*(*(u)).c23.re-s*(*(u)).c23.im; 
   (*(r)).c23.im=s*(*(u)).c23.re+c*(*(u)).c23.im; 
   (*(r)).c31.re=c*(*(u)).c31.re-s*(*(u)).c31.im; 
   (*(r)).c31.im=s*(*(u)).c31.re+c*(*(u)).c31.im; 
   (*(r)).c32.re=c*(*(u)).c32.re-s*(*(u)).c32.im; 
   (*(r)).c32.im=s*(*(u)).c32.re+c*(*(u)).c32.im; 
   (*(r)).c33.re=c*(*(u)).c33.re-s*(*(u)).c33.im; 
   (*(r)).c33.im=s*(*(u)).c33.re+c*(*(u)).c33.im;
   *u=*r;
   afree(r);
}

extern void twistbc(double *theta,double sign)
{
   int ibc,i,mu;
   su3_dble *u;
   int Lmu[4];

   if( (fabs(*(theta+0))<1e-10) && (fabs(*(theta+1))<1e-10) && 
       (fabs(*(theta+2))<1e-10) && (fabs(*(theta+3))<1e-10) )
     return ;

   message("Twisting the boundaries\n");

   ibc=query_flags(BCD_FLIPPED);
   if (ibc==1)
     flipbcd();

   Lmu[0]=NPROC0*L0;
   Lmu[1]=NPROC1*L1;
   Lmu[2]=NPROC2*L2;
   Lmu[3]=NPROC3*L3;
   u=pud[VOLUME/2][0];
   for (i=0;i<VOLUME/2;i++)
   {
      for(mu=0;mu<4;mu++){
      su3_appl_twist(u+8*i+2*mu  ,(double)(sign*(*(theta+mu))/Lmu[mu]));
      su3_appl_twist(u+8*i+2*mu+1,(double)(sign*(*(theta+mu))/Lmu[mu]));
      }
   }
   assign_ud2u();
   set_flags(NEW_UD);

   if (ibc==1)
     flipbcd();
}

