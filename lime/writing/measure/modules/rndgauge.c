/*******************************************************************************
*
* File rndgauge.c
*
* Copyright (C) 20010 Bastian Knippschild
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* the random gauge links are stored in su3_dble *g
*
* externally accassible functions are:
* 
* extern void transform_sd(spinor_dble *pk, spinor_dble *pl)
*   This functions allows a rnd gaguge transformation of a spinor. It is
*   used e.g. for the transformation of the pointsource befor applaying 
*   a smearing procedure on it!
*
* extern void transform_ud(su3_dble *input[VOLUME][4])
*   This function does the global random gauge transformation of the
*   gauge link array pud.
*
********************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "su3.h"
#include "start.h"
#include "global.h"
#include "random.h"
#include "update.h"
#include "measure.h"
#include "mpi.h"


static int nfc[8],ofs[8];
static su3_dble *g,*gbuf;

static void pack_gbuf(void)
{
   int n,ix,iy,io;

   nfc[0]=FACE0/2;
   nfc[1]=FACE0/2;
   nfc[2]=FACE1/2;
   nfc[3]=FACE1/2;
   nfc[4]=FACE2/2;
   nfc[5]=FACE2/2;
   nfc[6]=FACE3/2;
   nfc[7]=FACE3/2;

   ofs[0]=0;
   ofs[1]=ofs[0]+nfc[0];
   ofs[2]=ofs[1]+nfc[1];
   ofs[3]=ofs[2]+nfc[2];
   ofs[4]=ofs[3]+nfc[3];
   ofs[5]=ofs[4]+nfc[4];
   ofs[6]=ofs[5]+nfc[5];
   ofs[7]=ofs[6]+nfc[6];

   for (n=0;n<8;n++)
   {
      io=ofs[n];

      for (ix=0;ix<nfc[n];ix++)
      {
         iy=map[io+ix];
         gbuf[io+ix]=g[iy];
      }
   }
}

static void send_gbuf(void)
{
   int n,mu,np,saddr,raddr;
   int nbf;
   double *sbuf,*rbuf;
   MPI_Status stat;

   for (n=0;n<8;n++)
   {
      nbf=18*nfc[n];

      if (nbf>0)
      {
         mu=n/2;
         np=cpr[mu];

         if (n==(2*mu))
         {
            saddr=npr[n+1];
            raddr=npr[n];
         }
         else
         {
            saddr=npr[n-1];
            raddr=npr[n];
         }

         sbuf=(double*)(gbuf+ofs[n]);
         rbuf=(double*)(g+ofs[n]+VOLUME);

         if ((np|0x1)!=np)
         {
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,37,MPI_COMM_WORLD);
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,37,MPI_COMM_WORLD,&stat);
         }
         else
         {
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,37,MPI_COMM_WORLD,&stat);
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,37,MPI_COMM_WORLD);
         }
      }
   }
}


extern void transform_sd(spinor_dble *pk, spinor_dble *pl) {
	int ix;
	su3_dble gx;
	spinor_dble r, s;

	for (ix = 0; ix < VOLUME; ix++) {
		s = pk[ix];
		gx = g[ix];

		_su3_multiply(r.c1,gx,s.c1);
		_su3_multiply(r.c2,gx,s.c2);
		_su3_multiply(r.c3,gx,s.c3);
		_su3_multiply(r.c4,gx,s.c4);

		pl[ix] = r;
	}
}

static void random_g(void)
{
   su3_dble *gx,*gm;

   gm=g+VOLUME;

   for (gx=g;gx<gm;gx++)
      random_su3_dble(gx);

   if (BNDRY>0)
   {
      pack_gbuf();
      send_gbuf();
   }
}


extern void transform_ud(su3_dble *input[VOLUME][4])
{
   int ix,iy,mu;
   su3_dble *ub,u,v,w,gx,gxi,gy,gyi;

   g=amalloc(NSPIN*sizeof(su3_dble),4);

   if (BNDRY>0)
      gbuf=amalloc((BNDRY/2)*sizeof(su3_dble),4);

   random_g();

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      ub=input[ix][0];
      gx=g[ix];

      for (mu=0;mu<4;mu++)
      {
         iy=iup[ix][mu];
         gy=g[iy];
         u=ub[2*mu];
         _su3_dagger(gyi,gy);
         _su3_times_su3(v,u,gyi);
         _su3_times_su3(w,gx,v);
         ub[2*mu]=w;

         iy=idn[ix][mu];
         gy=g[iy];
         u=ub[2*mu+1];
         _su3_dagger(gxi,gx);
         _su3_times_su3(v,u,gxi);
         _su3_times_su3(w,gy,v);
         ub[2*mu+1]=w;
      }
   }

   set_flags(NEW_UD);
}

