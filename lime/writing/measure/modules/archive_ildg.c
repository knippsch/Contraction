
/*******************************************************************************
*
* File archive_ildg.c
*
* Copyright (C) 2008 Andreas Juettner
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Programs to read field configurations in the ILDG fomrat
*
* The externally accessible functions are
*
*   void import_cnfg_ildg(char *in)
*     Reads the global double-precision gauge field from the file "in".
* 
* The ILDG format uses the LIME file format. The LIME library is used
* in order to access the gauge configuration. Details on the ILDG format
* and the LIME library can be found here:
* http://www-zeuthen.desy.de/~pleiter/ildg/ildg-file-format-1.1.pdf
* http://usqcd.jlab.org/usqcd-docs/c-lime/lime_1p2.pdf
*
*******************************************************************************/

#define ARCH_ILDG_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "start.h"
#include "flags.h"
#include "update.h"
#include "misc.h"
#include "global.h"
#include <lime_config.h>
#include <lime.h>
#include <lime_fixed_types.h>
#include "measure.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static su3_dble *ubuf=NULL,*vbuf;
static su3 *vbuf_s;




static int check_machine(void)
{
   int ie;
   
   error_root(sizeof(stdint_t)!=4,1,"check_machine [archive_ildg.c]",
              "Size of a stdint_t integer is not 4");
   error_root(sizeof(double)!=8,1,"check_machine [archive_ildg.c]",
              "Size of a double is not 8");   

   ie=endianness();
   error_root(ie==UNKNOWN_ENDIAN,1,"check_machine [archive_ildg.c]",
              "Unkown endianness");

   return ie;
}


static void alloc_ubuf_ildg(int my_rank)
{
   if (my_rank==0)
   {
      ubuf=amalloc(8*sizeof(su3_dble),4);
      vbuf=ubuf+4;
      vbuf_s=amalloc(4*sizeof(su3),4);
   }
   else
      ubuf=amalloc(4*sizeof(su3_dble),4);

   error(ubuf==NULL,1,"alloc_ubuf [archive_ildg.c]",
         "Unable to allocate auxiliary array");
}

static void bswap_float(int n,float *a)
{
   unsigned char *ba,*bam,bas;

   ba=(unsigned char*)(a);
   bam=(unsigned char*)(a+n);

   for (;ba<bam;ba+=4)
   {
      bas=ba[3];
      ba[3]=ba[0];
      ba[0]=bas;

      bas=ba[2];
      ba[2]=ba[1];
      ba[1]=bas;      
   }
}



static void set_links_ildg(int iy,int y3)
{
   int iz[4],mu,ieo;
   int map[4]={3,0,1,2};
   su3_dble *u,*v;
   v=ubuf;
   ieo=0;

   if(ipt[iy+y3]<VOLUME/2){
    for (mu=0;mu<4;mu++){
     iz[mu] = iup[ipt[iy+y3]][mu];
     ieo=1;
    }
   }
   else {
    for (mu=0;mu<4;mu++){
     iz[mu] = ipt[iy+y3];
    } 
   }
   for (mu=0;mu<4;mu++){
    u = pud[iz[mu]][0];
    *(u+2*mu+ieo) = *(v+map[mu]);
   }
}
static void su3_copy_singletodouble(su3 *u,su3_dble *v){
   (*v).c11.re = (double)(*u).c11.re;
   (*v).c11.im = (double)(*u).c11.im;
   (*v).c12.re = (double)(*u).c12.re;
   (*v).c12.im = (double)(*u).c12.im;
   (*v).c13.re = (double)(*u).c13.re;
   (*v).c13.im = (double)(*u).c13.im;
   (*v).c21.re = (double)(*u).c21.re;
   (*v).c21.im = (double)(*u).c21.im;
   (*v).c22.re = (double)(*u).c22.re;
   (*v).c22.im = (double)(*u).c22.im;
   (*v).c23.re = (double)(*u).c23.re;
   (*v).c23.im = (double)(*u).c23.im;
   (*v).c31.re = (double)(*u).c31.re;
   (*v).c31.im = (double)(*u).c31.im;
   (*v).c32.re = (double)(*u).c32.re;
   (*v).c32.im = (double)(*u).c32.im;
   (*v).c33.re = (double)(*u).c33.re;
   (*v).c33.im = (double)(*u).c33.im;
}


extern void import_cnfg_ildg(char *in)
{
   int my_rank,np[4],n,ie,tag,status,prec=0;
   int i,l,ix,iy,iz,it,ixx;
   n_uint64_t read_bytes,nbytes=0;
   double plaq1,eps;
   MPI_Status stat;
   LimeReader *reader=NULL;
   FILE *fin;
   char lime_type_target[100]="ildg-binary-data";
   char *lime_type;
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   /* allocate IO-buffers */
   if (ubuf==NULL)
      alloc_ubuf_ildg(my_rank);

   error(pud[VOLUME/2][0]==NULL,1,"import_cnfg [archive_ildg.c]",
         "Double-precision gauge field is not allocated");
   ie=check_machine();   

   if (my_rank==0)
   {
     fin=fopen(in,"rb");
     error_root(fin==NULL,1,"import_cnfg  [archive_ildg.c]",
               	"Unable to open input file");
     reader=limeCreateReader(fin);
     error( reader == (LimeReader *)NULL ,1,"import_cnfg [archive_ildg.c]",
    		"Unable to open LimeReader");
     nbytes = (n_uint64_t)0;

     /* set the file-pointer to the beginning of the "ildg-binary-data" tag
      * in the LIME file */
     while( (status = limeReaderNextRecord(reader)) != LIME_EOF ){
       nbytes    = limeReaderBytes(reader);
       lime_type = limeReaderType(reader);
       error((status!=LIME_SUCCESS),1,"import_cnfg [archive_ildg.c]",
       	     "limeReaderNextRecord returned status %d and LIME_SUCCESS %d\n",
       status,LIME_SUCCESS);
       if (strcmp(lime_type,lime_type_target) != 0) continue;
        break;
     }
     /* Decide whether the gauge file is stored in single or
      * double precision format */
     if ((int)nbytes==8*VOLUME*NPROC*3*3*4*2){
       prec=8; 
       message("ILDG-file contains double precision gauge field\n");
     }
     else if ((int)nbytes==4*VOLUME*NPROC*3*3*4*2){
       prec=4; 
       message("ILDG-file contains single precision gauge field\n");
     }
     else
     error(1!=0,1,"import_cnfg [archive_ildg]",
		"Lattice geometry doesn't match size of gauge field %d %d\n");
   }
   /* Now read the gauge file in steps of 4 SU(3) matrices */
   for (it=0;it<N0;it++){
    for (iz=0;iz<N3;iz++){
     for (iy=0;iy<N2;iy++){
      for (ix=0;ix<N1;ix++){
	if(my_rank==0){
         if(prec==8){
          read_bytes = (n_uint64_t)(4*18*sizeof(double));
          status = limeReaderReadData((double *)vbuf,&read_bytes, reader);
	  bswap_double(4*18,(double *)(vbuf));/* ILDG always big endian */
         }
	 else if(prec==4){
          read_bytes = (n_uint64_t)(4*18*sizeof(float));
          status = limeReaderReadData((void *)vbuf_s,&read_bytes, reader);
	  bswap_float(4*18,(float *)(vbuf_s)); /* ILDG always big endian */
         for (i=0;i<4;i++)
	  su3_copy_singletodouble((vbuf_s+i),(vbuf+i));
         }
	}
        np[0]=it;
        np[1]=ix;
        np[2]=iy;
        np[3]=iz;
        n=ipr_global(np);
        MPI_Barrier(MPI_COMM_WORLD);
            if (n>0)
            {
               tag=mpi_tag();

               if (my_rank==0)
                  MPI_Send((double*)(vbuf),4*18,MPI_DOUBLE,n,tag,
                           MPI_COMM_WORLD);

               if (my_rank==n)
                  MPI_Recv((double*)(ubuf),4*18,MPI_DOUBLE,0,tag,
                           MPI_COMM_WORLD,&stat);
            }
            else if (my_rank==0)
               for (l=0;l<(4);l++)
                  ubuf[l]=vbuf[l];

            ixx = ((iy%L2)*L3 + (ix%L1)*L3*L2 + (it%L0)*L1*L2*L3);
            if (my_rank==n)
               set_links_ildg(ixx,iz);
      }
     }
    }
   }

   plaq1=plaq_sum_dble()/(double)(6*NPROC*VOLUME);
   message("Plaquette   %f\n",plaq1);
   message("Plaquette/3 %f\n",plaq1/3);
   message("Please check consistency of average plaquette by hand\n");
   eps=sqrt((double)(6*NPROC*VOLUME))*DBL_EPSILON;
}
