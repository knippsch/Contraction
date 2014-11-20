
/*******************************************************************************
*
* File twoptASCII.c
*
* Copyright (C) 2008 Andreas Jüttner
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*
* Syntax: twoptASCII  <lime_file> <# Message>
* or    : twoptASCII  <lime_file> <# Message>  <# channel> <# momentum>
*
* Converts 2pt and 3pt correlation functions stored in binary format
* in a lime-file into formatted ASCII ouput, e.g. 
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
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "measure.h"

int main(int argc,char *argv[])
{
  if( (argc < 3) || (argc == 4) || (argc > 5) ) {
    fprintf(stderr,"\n");
    fprintf(stderr,"Usage: \n");
    fprintf(stderr,"   %s <lime_file> <# Message>\n",argv[0]); 
    fprintf(stderr,"or alternatively \n");
    fprintf(stderr,"   %s <lime_file> <# Message> <# channel> <# momentum>\n", 
	argv[0]);
    fprintf(stderr,"\n");
  }
   if(argc==3){
    meson_ASCII_IO(argv[1],1,atoi(argv[2]));
   }
   else{
    meson_ASCII_IO_one(argv[1],1,atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));
   }
  exit(0);
}
