/*******************************************************************************
*
* File IOarchive.c
*
* Copyright (C) 2008-11 Andreas Juettner
*                       Stefano Capitani
*                       Bastian Knippschild
*                 
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Programs to read and save propagators and correlation functions 
* in LIME file format 
*
* The externally accessible functions are
*
* void export_init(char *out_file, int argc,char *argv[])
*  	set up files for correlation function export in LIME format
*
* void baryon_IOloop(char *tagname,char *out_file, complex_dble *corr, 
*	 int *mom, int *mu, int num_c)
*	export the baryon 3pt function to LIME format
*
* void meson_IOloop(char *tagname,char *out_file, complex_dble *corr, 
*	 int *mom, int *mu, int num_c)
*	export the hadrond 2- and 3pt-functions to LIME format
* 
* The ILDG format uses the LIME file format. The LIME library is used
* in order to access the gauge configuration. Details on the ILDG format
* and the LIME library can be found here:
* http://www-zeuthen.desy.de/~pleiter/ildg/ildg-file-format-1.1.pdf
* http://usqcd.jlab.org/usqcd-docs/c-lime/lime_1p2.pdf
*
*******************************************************************************/

#define IOARCHIVE_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "misc.h"
#include "start.h"
#include "global.h"
#include "mpi.h"
#include "measure.h"
#include "linalg.h"
#include <lime_config.h>
#include <lime.h>


#define MAXBUF 1048576
n_uint64_t mino(n_uint64_t i, n_uint64_t j){
  return i < j ? i : j;
}
/* Discover how many bytes there are in the file */
/* taken over 
*  from http://usqcd.jlab.org/usqcd-software/c-lime/c-lime/examples/lime_pack.c
*/
long file_size(FILE *fp)
{
  long oldpos = ftell(fp);
  long length;
  
  if (fseek(fp, 0L,SEEK_END) == -1)
    return -1;
  
  length = ftell(fp);
  
  return ( fseek(fp,oldpos,SEEK_SET) == -1 ) ? -1 : length;
}
void export_init(char *out_file, int argc,char *argv[])
{
    int MB_flag,ME_flag,ifile,ir;
    FILE *fin=NULL;
    char input[MAX_INFILE_SIZE];
    LimeWriter *w;
    LimeRecordHeader *h;
    n_uint64_t bytes;


    /* first read the input-file of this run and determine its size */
     ifile = find_opt(argc,argv,"-i");
     fin=freopen(argv[ifile+1],"r",stdin);
     error_root(fin==NULL,1,"read_infile [measure1.c]",
                 "Unable to open input file");
     bytes=file_size(fin);
     error_root(bytes>=MAX_INFILE_SIZE,1,"read_infile [measure1.c]",
                "Increase MAX_INFILE_SIZE");
     ir=fread(&input,sizeof(char), bytes, fin);
     fclose(fin);
     message("Writing %d bytes of header to out_file: %s \n",
                (int)bytes,out_file);
    /* now open the file that will contain the correlator */
     fin=freopen(out_file,"w",stdin);
     w = limeCreateWriter(fin);
      /* write the input file */
        MB_flag = 1; ME_flag = 1;
        h = limeCreateHeader(MB_flag,ME_flag,"input",bytes);
        limeWriteRecordHeader(h,w);
        limeWriteRecordData(input,&bytes,w);
    limeDestroyWriter(w);
    fclose(fin);
}
void export_meson(char *tagname, char *out_file,twopt_mom *icorr2pt, 
				unsigned long length)
{
    int MB_flag,ME_flag;
    int status=EXIT_SUCCESS;
    FILE *fin=NULL;
    LimeWriter *w;
    LimeRecordHeader *h;
    n_uint64_t i;
    i=(n_uint64_t)length;

    /* now open the file that will contain the correlator */
     fin=fopen(out_file,"a");
     w = limeCreateWriter(fin);
      /* write the correlators */
        MB_flag = 1; ME_flag = 1;
        h = limeCreateHeader(MB_flag,ME_flag,tagname,i);
        limeWriteRecordHeader(h,w);
        limeWriteRecordData(icorr2pt,&i,w);
        error_root(status != LIME_SUCCESS,1,"read_infile [measure1.c]",
                 "LIME write error %d\n", status);
    limeDestroyWriter(w);
    fclose(fin);
}

void meson_IOloop(char *tagname,char *out_file, complex_dble *corr,int nmom_max, int *mom, int *mu, int num_c)
{   
   int i,ii,imu1,imu2,ip,t,my_rank;
   double wt1,wt2;
   twopt_mom *corr2pt=NULL;
   corr2pt   = amalloc(num_c*num_c*nmom_max*sizeof(twopt_mom),3);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   if (my_rank==0){
     i=0;
     ii=-1;
     wt1=MPI_Wtime();
     for(imu1=0;imu1<num_c;imu1++)
     for(imu2=0;imu2<num_c;imu2++)
     for(ip=0;ip<nmom_max;ip++){
       ii+=1;
       for(t=0;t<NPROC0*L0;t++){
        /* preparation for binary output */
        (*(corr2pt+ii)).musrc = mu[imu1];
        (*(corr2pt+ii)).musnk = mu[imu2];
        (*(corr2pt+ii)).p[0] = (*(mom+ip*3));
        (*(corr2pt+ii)).p[1] = (*(mom+ip*3+1));
        (*(corr2pt+ii)).p[2] = (*(mom+ip*3+2));
        (*(corr2pt+ii)).corr[t].re = (*(corr+i)).re;
        (*(corr2pt+ii)).corr[t].im = (*(corr+i)).im;
        /* ASCII output */
 /*      message("%d %s %d %d %d %d %d %d %.16e %.16e\n",
         	ii,tagname,mu[imu1],mu[imu2],
                 (*(mom+ip*3)),
                 (*(mom+ip*3+1)),
                 (*(mom+ip*3+2)),
                 t,(*(corr+i)).re,(*(corr+i)).im);
  */
       i++;
       }
     }
     export_meson(tagname,out_file,corr2pt,
         	(num_c*num_c*nmom_max)*sizeof(twopt_mom));
     wt2=MPI_Wtime();
   /*  message("time for IO %.2e sec\n",wt2-wt1);*/
   }
   afree(corr2pt);
}

void meson_ASCII_IO(char *filename, int rec_seek, int msg_seek,int nmom_max)
{ 
    int i,j,k=0,l=0;
    char buf[MAXBUF];
    LimeReader *r;
    int status;
    twopt_mom *corr;
    n_uint64_t nbytes, bytes_left, bytes_to_copy, read_bytes;
    int rec, msg;
    char *lime_type;
    size_t bytes_pad;
    int MB_flag, ME_flag;
    FILE *fin=NULL;
  
    fin=fopen(filename,"r");
    error_root(fin == (FILE *)NULL,0,"meson_ASCII_IO [IOarchive.c]",
               "IOarchive: Unable to open file %s for reading\n", filename);
    /* Open LIME reader */
    r = limeCreateReader(fin);
    error_root(r == (LimeReader *)NULL,0,"meson_ASCII_IO [IOarchive.c]",
               "Unable to open LimeReader\n");
    /* Loop over records */
    rec = 0;
    msg = 0;
    while( (status = limeReaderNextRecord(r)) != LIME_EOF ){
    error_root(status != LIME_SUCCESS,0,"meson_ASCII_IO [IOarchive.c]",
             "IOarchive: limeReaderNextRecord returned status = %d\n",status);
    

    nbytes    = limeReaderBytes(r);
    lime_type = limeReaderType(r);
    bytes_pad = limeReaderPadBytes(r);
    MB_flag   = limeReaderMBFlag(r);
    ME_flag   = limeReaderMEFlag(r);

    /* Update message and record numbers */
    if(MB_flag == 1){
      msg++;
      rec = 0;
    }

    rec++;

    /* Skip to next record until target record is reached */
    if (msg != msg_seek || rec != rec_seek) continue;
    
    /* Buffered copy */
    bytes_left = nbytes;
    while(bytes_left > (n_uint64_t)0){
      bytes_to_copy = mino((n_uint64_t)(nmom_max*sizeof(twopt_mom)),bytes_left);
      read_bytes = bytes_to_copy;
      error_root((int)read_bytes>MAXBUF,0,"meson_ASCII_IO [IOarchive.c]",
             "read_bytes exceeds MAXBUF");
      status = limeReaderReadData((void *)buf, &read_bytes, r);
      error_root(status<0 && status !=LIME_EOR,1,"meson_ASCII_IO [IOarchive.c]",
             "LIME read error occurred: status= %d", status);
      error_root(read_bytes != bytes_to_copy,1,"meson_ASCII_IO [IOarchive.c]",
             "Read error %lld bytes wanted,%lld read\n",
                (unsigned long long)nbytes, (unsigned long long)read_bytes);
    
      /* Print to stdout */
      corr=(twopt_mom*)buf;  
      l=0;
      for(i=0;i<nmom_max;i++){
       for(j=0;j<NPROC0*L0;j++)
        {
            printf("%d %s %d %d %d %d %d %d %.16e %.16e\n",
  		k,lime_type,
		(*(corr+l)).musrc,
		(*(corr+l)).musnk,
		(*(corr+l)).p[0],
		(*(corr+l)).p[1],
		(*(corr+l)).p[2],
		j,
		(*(corr+l)).corr[j].re,
		(*(corr+l)).corr[j].im);
        }
       l++;
       k++;
      }
      bytes_left -= bytes_to_copy;
    }
    /* Quit at this record */
    break;
  }

  limeDestroyReader(r);
  fclose(fin);

}
void baryon_IOloop(char *tagname,char *out_file, complex_dble *corr, 
	int nmom_max, int *mom, int *mu, int num_c)
{   
   int i,ii,imu1,ip,t,my_rank;
   twopt_mom *corr3pt=NULL;
   corr3pt = amalloc(num_c*nmom_max*sizeof(twopt_mom),3);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   if (my_rank==0){
    i=0;
    ii=-1;
    for(imu1=0;imu1<num_c;imu1++)
    for(ip=0;ip<nmom_max;ip++){
      ii+=1;
      for(t=0;t<NPROC0*L0;t++){
       /* preparation for binary output */
       (*(corr3pt+ii)).musrc = 0;
       (*(corr3pt+ii)).musnk = mu[imu1];
       (*(corr3pt+ii)).p[0] = (*(mom+ip*3));
       (*(corr3pt+ii)).p[1] = (*(mom+ip*3+1));
       (*(corr3pt+ii)).p[2] = (*(mom+ip*3+2));
       (*(corr3pt+ii)).corr[t].re = (*(corr+i)).re;
       (*(corr3pt+ii)).corr[t].im = (*(corr+i)).im;
       /* ASCII output */
/*    message("%d %s %d %d %d %d %d %d %.16e %.16e\n",
	      ii,tagname,5,mu[imu1],
                (*(mom+ip*3)),
                (*(mom+ip*3+1)),
                (*(mom+ip*3+2)),
                t,(*(corr+i)).re,(*(corr+i)).im);
*/
    i++;
    }}
    export_meson(tagname,out_file,corr3pt,(num_c*nmom_max)*sizeof(twopt_mom));
    }
    afree(corr3pt);
}


void write_fsv(char *out,full_spinor_dble *sv)
{
   int ldat[9],iw,i,j;
   double norm;
   FILE *fout=NULL;
   spinor_dble **wsd;
   wsd=reserve_wsd(1);

   MPI_Comm_rank(MPI_COMM_WORLD,&(ldat[8]));
   sprintf(out,"%s_rank%d",out,ldat[8]);
   fout=fopen(out,"wb");
   error_loc(fout==NULL,1,"write_fsv [IOarchive.c]",
             "Unable to open output file");
   error_chk();

   ldat[0]=NPROC0;
   ldat[1]=NPROC1;
   ldat[2]=NPROC2;
   ldat[3]=NPROC3;

   ldat[4]=L0;
   ldat[5]=L1;
   ldat[6]=L2;
   ldat[7]=L3;

   iw=fwrite(ldat,sizeof(int),9,fout);
   /* compute 12 norm squares for spinor_dble's */
   for (i=1;i<=4;i++){ 
    for (j=1;j<=3;j++){
    copy_fs_sv(VOLUME,sv,wsd[0],i,j);
    norm=norm_square_dble(VOLUME,1,wsd[0]);
    iw+=fwrite(&norm,sizeof(double),1,fout);
    }
   }
   iw+=fwrite(sv,sizeof(full_spinor_dble),VOLUME,fout);

   error_loc((iw!=(21+VOLUME))||(ferror(fout)!=0),1,
             "write_fsv [IOarchive.c]","Incorrect write count or write error");
   error_chk();
   fclose(fout);
   release_wsd();
}


void read_fsv(char *in,full_spinor_dble *sv)
{
   int ldat[9],ir,n,i,j,k;
   double norm0[12],norm1,eps;
   FILE *fin=NULL;
   spinor_dble **wsd;
   wsd=reserve_wsd(1);

   MPI_Comm_rank(MPI_COMM_WORLD,&n);
   sprintf(in,"%s_rank%d",in,n);
   
   fin=fopen(in,"rb");
   error_loc(fin==NULL,1,"read_fsv [IOarchive.c]",
             "Unable to open input file");
   error_chk();

   ir=fread(ldat,sizeof(int),9,fin);

   error((ldat[0]!=NPROC0)||(ldat[1]!=NPROC1)||
         (ldat[2]!=NPROC2)||(ldat[3]!=NPROC3)||
         (ldat[4]!=L0)||(ldat[5]!=L1)||(ldat[6]!=L2)||(ldat[7]!=L3)||
         (ldat[8]!=n),1,"read_fsv [sarchive.c]","Unexpected lattice data");
   for (i=0;i<12;i++){
    ir+=fread(&(norm0[i]),sizeof(double),1,fin);
   }
   ir+=fread(sv,sizeof(full_spinor_dble),VOLUME,fin);

   error_loc((ir!=(21+VOLUME))||(ferror(fin)!=0),1,
             "read_fsv [IOarchive.c]","Incorrect read count or read error");
   error_chk();
   fclose(fin);


   eps=sqrt((double)(NPROC*VOLUME))*DBL_EPSILON;
   k=0;
   for (i=1;i<=4;i++){ 
    for (j=1;j<=3;j++){
    copy_fs_sv(VOLUME,sv,wsd[0],i,j);
    norm1=norm_square_dble(VOLUME,1,wsd[0]);
    error(fabs(norm1-norm0[k])>(eps*norm0[i]),1,"read_fsv [IOarchive.c]",
         "Norm test failed: norm1=%.4e norm0[%d]=%.4e eps=%.4e",
		norm1,norm0[k],eps);
    k++;
    }
   }
   release_wsd();

   
}
void read_fsv2(char *stem, int prop_dump, int *pim, int nprop,int *isv,
	full_spinor_dble **prop_pt, full_spinor_dble **sv, 
	full_spinor_dble **svaux)
{ 
  int i;
  char in_fsv[NAME_SIZE];
  for (i=0;i<nprop;i++){
   if (prop_dump == 0){
    prop_pt[i]=sv[*(isv+i)];
   } 
   else{
    if((*(pim+i))==(*(isv+i))){ /* prop already loaded? */
     prop_pt[i]=svaux[i];
    }
    else{ /* load it */
     sprintf(in_fsv,"%s%d",stem,*(isv+i));
     read_fsv(in_fsv,svaux[i]);
     prop_pt[i]=svaux[i];
     (*(pim+i))=(*(isv+i));
    }
   } 
  }
}
