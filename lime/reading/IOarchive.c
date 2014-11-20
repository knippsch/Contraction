/*******************************************************************************
*
* File IOarchive.c
*
* Copyright (C) 2008 Andreas Juettner
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Programs to read and save propagators and correlation functions 
* in LIME file format 
*
* The externally accessible functions are
*
* void meson_ASCII_IO(char *filename, int rec_seek, int msg_seek)
*	converts entries in lime file with 
*	record number 		rec_sec 	and
*	message number 		msg_seek 	to ASCII.

* void meson_ASCII_IO_one(char *filename, int rec_seek, int msg_seek, 
*	int chan_seek, int mom_seek)
*	converts entries in lime file with 
*	record number 		rec_sec 	and
*	message number 		msg_seek 	and
*	channel number  	chan_seek	and
*	momentum channel	mom_seek	to ASCII.
* 
* The ASCII ouput is formatted as:
*  /-- running integer number
*  |     /-- tag name
*  |     |  /-- source gamma index
*  |     |  | /-- sink gamma index
*  |     |  | |  /-- three momentum 
*  |     |  | |  |    /-- time slice 
*  |     |  | |  |    |        /-- real part
*  |     |  | |  |    |        |                      /-- imag part
*  |     |  | |  |    |        |                      |
*  |     |  | | ----- |        |                      |
* 18464 2pt 5 5 0 0 0 0 8.5233228411805362e-01 0.0000000000000000e+00
* 18465 2pt 5 5 0 0 0 1 5.6879604419396254e-02 0.0000000000000000e+00
* 18466 2pt 5 5 0 0 0 2 1.2484384424550077e-02 0.0000000000000000e+00
* 18467 2pt 5 5 0 0 0 3 6.7376590573064755e-03 0.0000000000000000e+00
* 18468 2pt 5 5 0 0 0 4 5.9425073611925660e-03 0.0000000000000000e+00
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
#include <stdarg.h>
#include "measure.h"
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



void error_root(int test,int no,char *name,char *format,...)
{
   va_list args;
   if ((test!=0))
   {

      printf("\nError in %s:\n",name);
      va_start(args,format);
      vprintf(format,args);
      va_end(args);

      printf("\nProgram aborted\n\n");
      fflush(stdout);
      exit(no);
  }
}


void meson_ASCII_IO_one_bin(double *pc, int reim, char *filename, int rec_seek,
			    int msg_seek, int chan_seek, int mom_seek)
{ 
    int j,l=0;
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
  
    corr=malloc(NMOM_MAX*sizeof(twopt_mom));
  
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
    l=0;
    while(bytes_left > (n_uint64_t)0){
      bytes_to_copy=mino((n_uint64_t)(NMOM_MAX*sizeof(twopt_mom)),bytes_left);
      read_bytes =  bytes_to_copy;
      error_root((int)read_bytes>MAXBUF,0,"meson_ASCII_IO [IOarchive.c]",
             "read_bytes exceeds MAXBUF");
      status = limeReaderReadData((void *)buf, &read_bytes, r);
      error_root(status<0 && status !=LIME_EOR,1,"meson_ASCII_IO [IOarchive.c]",
             "LIME read error occurred: status= %d", status);
      error_root(read_bytes != bytes_to_copy,1,"meson_ASCII_IO [IOarchive.c]",
             "Read error %lld bytes wanted,%lld read\n",
                (unsigned long long)nbytes, (unsigned long long)read_bytes);
    
      /* Print to stdout */
      if (l==chan_seek){
      corr=(twopt_mom*)buf;  
      for(j=0;j<TT;j++)
      {
	/*	printf("%d %s %d %d %d %d %d %d %.16e %.16e\n",
	 *       	k,lime_type,
	 *       	(*(corr+mom_seek)).musrc,
	 *       	(*(corr+mom_seek)).musnk,
	 *       	(*(corr+mom_seek)).p[0],
	 *       	(*(corr+mom_seek)).p[1],
	 *       	(*(corr+mom_seek)).p[2],
	 *       	j,
	 *       	(*(corr+mom_seek)).corr[j].re,
	 *       	(*(corr+mom_seek)).corr[j].im);
	 */
	if (reim==0)
	  pc[j]=(*(corr+mom_seek)).corr[j].re;
	else
	  pc[j]=(*(corr+mom_seek)).corr[j].im;
      }
      }
      bytes_left -= bytes_to_copy;
      l=l+1;
    }
    /* Quit at this record */
    break;
  }
  limeDestroyReader(r);
  fclose(fin);

}


void meson_ASCII_IO_one(char *filename, int rec_seek, int msg_seek, 
	int chan_seek, int mom_seek)
{ 
    int j,k=0,l=0;
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
  
    corr=malloc(NMOM_MAX*sizeof(twopt_mom));
  
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
    l=0;
    while(bytes_left > (n_uint64_t)0){
      bytes_to_copy=mino((n_uint64_t)(NMOM_MAX*sizeof(twopt_mom)),bytes_left);
      read_bytes =  bytes_to_copy;
      error_root((int)read_bytes>MAXBUF,0,"meson_ASCII_IO [IOarchive.c]",
             "read_bytes exceeds MAXBUF");
      status = limeReaderReadData((void *)buf, &read_bytes, r);
      error_root(status<0 && status !=LIME_EOR,1,"meson_ASCII_IO [IOarchive.c]",
             "LIME read error occurred: status= %d", status);
      error_root(read_bytes != bytes_to_copy,1,"meson_ASCII_IO [IOarchive.c]",
             "Read error %lld bytes wanted,%lld read\n",
                (unsigned long long)nbytes, (unsigned long long)read_bytes);
    
      /* Print to stdout */
      if (l==chan_seek){
      corr=(twopt_mom*)buf;  
      for(j=0;j<TT;j++)
       {
           printf("%d %s %d %d %d %d %d %d %.16e %.16e\n",
       	k,lime_type,
       	(*(corr+mom_seek)).musrc,
       	(*(corr+mom_seek)).musnk,
       	(*(corr+mom_seek)).p[0],
       	(*(corr+mom_seek)).p[1],
       	(*(corr+mom_seek)).p[2],
       	j,
       	(*(corr+mom_seek)).corr[j].re,
       	(*(corr+mom_seek)).corr[j].im);
       }
       }
      bytes_left -= bytes_to_copy;
      l=l+1;
    }
    /* Quit at this record */
    break;
  }
  limeDestroyReader(r);
  fclose(fin);

}


void meson_ASCII_IO(char *filename, int rec_seek, int msg_seek)
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
  
    corr=malloc(NMOM_MAX*sizeof(twopt_mom));
  
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
      bytes_to_copy = mino((n_uint64_t)(NMOM_MAX*sizeof(twopt_mom)),bytes_left);
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
      for(i=0;i<NMOM_MAX;i++){
       for(j=0;j<TT;j++)
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
