
/*******************************************************************************
*
* File measure8.c
*
* Copyright (C) 2008-11 Andreas Juettner 
*                   and Stefano Capitani
*                   and Bastian Knippschild
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*
* Syntax: measure8 -i <filename> 
*
* Sample input file is measure8.in
*
* Does mesonic 2pt and 3pt correlation functions
*
* And also 2pt and 3pt correlations of baryons
*
*******************************************************************************/

#define MAIN_PROGRAM
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


static int my_rank,endian;
static int id,pid,bc,nc_first,nc_delta,nc_last,level,seed;
static int bs_sap[4],nkv_sap,nmr_sap,ncy_sap,nmx;
static int bs_dfl[4],Ns,rinv,rnmr,rncy,rvnmx;
static int vnkv,vnmx,vdnmx,nmx;
static int nohits,dtsrc,noprops,noextprops,propidfirstleg[NMAX_PARM],
	          extpropid[NMAX_PARM],propid[NMAX_PARM];
static int notwopt,twoptid[NMAX_PARM],twoptprop1[NMAX_PARM],twoptprop2[NMAX_PARM];
static int nothreept,threeptid[NMAX_PARM],threeptprop1[NMAX_PARM],threeptprop2[NMAX_PARM];
static int pos[NMAX_PARM][4],t_snk[NMAX_PARM],prop_dump;

static int mom_insert[NMAX_PARM],gamma_insert[NMAX_PARM];

static double beta,kappa[NMAX_PARM],csw,rkappa,kappasea;
static double rvres,vres,vdres,res;
static double theta[NMAX_PARM][4];

static int source_shift_id[NMAX_PARM], source_sh_no[NMAX_PARM], 
           no_source_sh, source_shift[NMAX_PARM][4];

static int nosmearparm, smearid[NMAX_PARM], smearno[NMAX_PARM];
static double alpha[NMAX_PARM][3];
static char smearing_type[NMAX_PARM][7];

static int no_jac_parm, jac_parm_id[NMAX_PARM], 
           jac_N[NMAX_PARM], jac_link_smear_no[NMAX_PARM];
static double jac_kappa[NMAX_PARM]; 
static char jac_link_smear_type[NMAX_PARM][7];

static int max_mom,nmom_max;
static int mu2pt[NMAX_DIRAC], mu3pt[NMAX_DIRAC], mub2pt[NMAX_DIRAC], mub3pt[NMAX_DIRAC];
static int num2pt, num3pt, numb2pt, numb3pt;

static int twopt_sk_sm[NMAX_PARM], twopt_b_sk_sm[NMAX_PARM], threept_sk_sm[NMAX_PARM],
           threept_b_sk_sm[NMAX_PARM], extprop_sk_sm[NMAX_PARM], extprop_b_sk_sm[NMAX_PARM],
           no_sk_sm[21][NMAX_PARM],sk_sm_nb;
static int noextprops_b,propidleg1_b[NMAX_PARM],extpropid_b[NMAX_PARM],
           propidleg2_b[NMAX_PARM],propidleg3_b[NMAX_PARM],mom_insert_b[NMAX_PARM],pol3pt[NMAX_PARM];
static int notwopt_b,twoptid_b[NMAX_PARM],twoptprop1_b[NMAX_PARM],
           twoptprop2_b[NMAX_PARM], twoptprop3_b[NMAX_PARM],pol2pt[NMAX_PARM];
static int nothreept_b,threeptid_b[NMAX_PARM],threeptprop1_b[NMAX_PARM],
           threeptprop2_b[NMAX_PARM];
static int t_snk_b[NMAX_PARM];

static double kappa_b[NMAX_PARM];
static double theta_b[NMAX_PARM][4]; 

static char log_dir[NAME_SIZE],cnfg_dir[NAME_SIZE],corr_dir[NAME_SIZE],
		          sfld_dir[NAME_SIZE],start[7];
static char log_file[NAME_SIZE],log_save[NAME_SIZE];
static char cnfg_file[NAME_SIZE],end_file[NAME_SIZE];
static char src_type[NMAX_PARM][7];
static char nbase[NAME_SIZE],ncnfg[NAME_SIZE];
static FILE *fin=NULL,*flog=NULL,*fend=NULL;

static lat_parms_t lat;
static gcr_parms_t gcr;
static dfl_parms_t dfl;

static void read_infile(int argc,char *argv[])
{
   int ifile,i;
   
   if (my_rank==0)
   {

      flog=freopen("MPI_Errors","w",stdout);

      ifile=find_opt(argc,argv,"-i");
      endian=endianness();
      error_root((ifile==0)||(ifile==(argc-1)),1,"read_infile [measure8.c]",
                 "Syntax: measure8 -i <filename>");

      error_root(endian==UNKNOWN_ENDIAN,1,"read_infile [measure8.c]",
                 "Machine has unknown endianness");
      
      fin=freopen(argv[ifile+1],"r",stdin);
      error_root(fin==NULL,1,"read_infile [measure8.c]",
                 "Unable to open input file");

      read_line("id","%d\n",&id);
      read_line("pid","%d\n",&pid);
      read_line("start","%s\n",start);
      read_line("log_dir","%s\n",log_dir);
      read_line("cnfg_dir","%s\n",cnfg_dir);
      read_line("corr_dir","%s\n",corr_dir);
      read_line("sfld_dir","%s\n",sfld_dir);
      read_line("prop_dump","%d\n",&prop_dump);
      read_line("beta","%lf\n",&beta);
      read_line("kappasea","%lf\n",&kappasea);
      read_line("csw","%lf\n",&csw);
      read_line("bc","%d\n",&bc);
      read_line("nc_first","%d\n",&nc_first);
      read_line("nc_delta","%d\n",&nc_delta);
      read_line("nc_last","%d\n",&nc_last);
      read_line("level","%d\n",&level);
      read_line("seed","%d\n",&seed);

      read_line("bs_sap","%d %d %d %d\n",
                &bs_sap[0],&bs_sap[1],&bs_sap[2],&bs_sap[3]);
      read_line("nkv_sap","%d\n",&nkv_sap);
      read_line("nmr_sap","%d\n",&nmr_sap);      
      read_line("ncy_sap","%d\n",&ncy_sap);
      read_line("nmx","%d\n",&nmx);
      read_line("res","%lf\n",&res);

      read_line("bs_dfl","%d %d %d %d\n",
                &bs_dfl[0],&bs_dfl[1],&bs_dfl[2],&bs_dfl[3]);
      read_line("Ns","%d\n",&Ns);
      read_line("rkappa","%lf\n",&rkappa);
      read_line("rinv","%d\n",&rinv);
      read_line("rnmr","%d\n",&rnmr);
      read_line("rncy","%d\n",&rncy);
      read_line("rvnmx","%d\n",&rvnmx);      
      read_line("rvres","%lf\n",&rvres);
      read_line("vnkv","%d\n",&vnkv);
      read_line("vnmx","%d\n",&vnmx);
      read_line("vres","%lf\n",&vres);
      read_line("vdnmx","%d\n",&vdnmx);
      read_line("vdres","%lf\n",&vdres);
      read_line("no_hits","%d\n",&nohits);
      read_line("dtsrc","%d\n",&dtsrc);
	   
      /**** sourece shifts ****/
      read_line("no_source_sh", "%d\n", &no_source_sh);
      for(i=0;i<no_source_sh;i++){
	       read_line("source_shift", "%d %d %d %d %d\n", &source_shift_id[i],
                   &source_shift[i][0],&source_shift[i][1],&source_shift[i][2],
                   &source_shift[i][3]);
      }
	  
      /**** smearings ****/
      read_line("no_smear_parm", "%d\n", &nosmearparm);
      for(i=0;i<nosmearparm;i++){
	       read_line("smearing","%d %s",&smearid[i],smearing_type[i]);
	       if(strcmp(smearing_type[i], "APE")==0){
		        read_line("","%lf\n",&alpha[i][0]);
		        alpha[i][1] = 0.0;
		        alpha[i][2] = 0.0;
	       }
	       else if(strcmp(smearing_type[i], "HYP")==0){
	         read_line("","%lf %lf %lf\n",&alpha[i][0],&alpha[i][1],&alpha[i][2]);
        }
      }
		
      read_line("no_jac_parm", "%d\n", &no_jac_parm);
      for(i=0;i<no_jac_parm;i++){
	       read_line("jacobi", "%d %s %d %lf %d", &jac_parm_id[i], 
		                         jac_link_smear_type[i],
		                         &jac_N[i], &jac_kappa[i], &jac_link_smear_no[i]);
      }

      /**** momenta and gamma matrices for contractions ****/
      read_line("num2pt","%d\n",&num2pt);
      for (i=0;i<num2pt;i++) scanf("%d\n",&mu2pt[i]);

      read_line("num3pt","%d\n",&num3pt);
      for (i=0;i<num3pt;i++) scanf("%d\n",&mu3pt[i]);

      read_line("numb2pt","%d\n",&numb2pt);
      for (i=0;i<numb2pt;i++) scanf("%d\n",&mub2pt[i]);

      read_line("numb3pt","%d\n",&numb3pt);
      for (i=0;i<numb3pt;i++) scanf("%d\n",&mub3pt[i]);

      read_line("max_mom","%d\n",&max_mom);
      nmom_max = (2*max_mom+1)*(2*max_mom+1)*(2*max_mom+1);


      /**** propagators and contractions ****/
      read_line("no_props","%d\n",&noprops);
      for(i=0;i<noprops;i++){
	       theta[i][0]=0.0;
	       read_line("prop","%d %lf %s %d %d %d %d %lf %lf %lf %d %d\n",
		                       &propid[i],&kappa[i],src_type[i], 
		                       &pos[i][0],&pos[i][1],&pos[i][2],&pos[i][3],
		                       &theta[i][1],&theta[i][2],&theta[i][3],
		                       &smearno[i],&source_sh_no[i]);
      }

      read_line("no_ext_props","%d",&noextprops);
      for(i=0;i<noextprops;i++){
        theta[i+noprops][0]=0.0;
        read_line("ext_prop","%d %d %d %lf %d %d %d %lf %lf %lf\n",&extpropid[i],
		                           &propidfirstleg[i], &extprop_sk_sm[i], &kappa[i+noprops], 
		                           &t_snk[i],&mom_insert[i],&gamma_insert[i],
		                           &theta[i+noprops][1],&theta[i+noprops][2],&theta[i+noprops][3]);
      }

      read_line("no_2pt","%d\n",&notwopt);
      for(i=0;i<notwopt;i++)
        read_line("2pt","%d %d %d %d\n",&twoptid[i],&twoptprop1[i],&twoptprop2[i],&twopt_sk_sm[i]);
        
      read_line("no_3pt","%d",&nothreept);
      for(i=0;i<nothreept;i++)
        read_line("3pt","%d %d %d %d",&threeptid[i],&threeptprop1[i],
		                                    &threeptprop2[i],&threept_sk_sm[i]);

      read_line("no_ext_props_b","%d\n",&noextprops_b);
      for (i=0;i<noextprops_b;i++){
        theta_b[i][0]=0.0;
        read_line("ext_prop_b","%d %d %d %d %d %lf %d %d %lf %lf %lf %d\n",
		                    &extpropid_b[i],&propidleg1_b[i],&propidleg2_b[i],
		                    &propidleg3_b[i],&extprop_b_sk_sm[i],&kappa_b[i],
		                    &t_snk_b[i],&mom_insert_b[i],&theta_b[i][1],&theta_b[i][2],
                      &theta_b[i][3],&pol3pt[i]);
      }


      read_line("no_2pt_b","%d\n",&notwopt_b);
      for(i=0;i<notwopt_b;i++)
        read_line("2pt_b","%d %d %d %d %d %d\n",&twoptid_b[i],&twoptprop1_b[i],
                 &twoptprop2_b[i],&twoptprop3_b[i],&twopt_b_sk_sm[i],&pol2pt[i]);

      read_line("no_3pt_b","%d",&nothreept_b);
      for(i=0;i<nothreept_b;i++)
        read_line("3pt_b","%d %d %d %d",&threeptid_b[i],&threeptprop1_b[i],
                  &threeptprop2_b[i],&threept_b_sk_sm[i]);

      fclose(fin);
   }

   /**** sending everything to all other processes ****/
   MPI_Bcast(&endian,1,MPI_INT,0,MPI_COMM_WORLD);
   
   MPI_Bcast(&id,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&pid,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(start,7,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(log_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(cnfg_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(corr_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(sfld_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   
   MPI_Bcast(&prop_dump,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&beta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&kappasea,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&csw,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nc_first,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nc_delta,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nc_last,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&level,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&seed,1,MPI_INT,0,MPI_COMM_WORLD);
 
   MPI_Bcast(bs_sap,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nkv_sap,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmr_sap,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ncy_sap,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmx,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&res,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   MPI_Bcast(bs_dfl,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&Ns,1,MPI_INT,0,MPI_COMM_WORLD);   
   MPI_Bcast(&rkappa,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&rinv,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&rnmr,1,MPI_INT,0,MPI_COMM_WORLD);   
   MPI_Bcast(&rncy,1,MPI_INT,0,MPI_COMM_WORLD);   
   MPI_Bcast(&rvnmx,1,MPI_INT,0,MPI_COMM_WORLD);   
   MPI_Bcast(&rvres,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&vnkv,1,MPI_INT,0,MPI_COMM_WORLD);   
   MPI_Bcast(&vnmx,1,MPI_INT,0,MPI_COMM_WORLD);   
   MPI_Bcast(&vres,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&vdnmx,1,MPI_INT,0,MPI_COMM_WORLD);   
   MPI_Bcast(&vdres,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&nohits,1,MPI_INT,0,MPI_COMM_WORLD);   
   MPI_Bcast(&dtsrc,1,MPI_INT,0,MPI_COMM_WORLD); 
      
   MPI_Bcast(&no_source_sh, 1, MPI_INT, 0, MPI_COMM_WORLD); 
   MPI_Bcast(&source_sh_no, no_source_sh, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&source_shift, 4*no_source_sh, MPI_INT, 0, MPI_COMM_WORLD);
   
   MPI_Bcast(&nosmearparm, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&smearid, nosmearparm, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&alpha, 3*nosmearparm, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   for(i=0;i<NMAX_PARM;i++)
	    MPI_Bcast(smearing_type[i],7,MPI_CHAR,0,MPI_COMM_WORLD);
	
   MPI_Bcast(&no_jac_parm, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&jac_parm_id, no_jac_parm, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&jac_N, no_jac_parm, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&jac_link_smear_no, no_jac_parm, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&jac_kappa, no_jac_parm, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   for(i=0; i<NMAX_TWIST; i++)
	    MPI_Bcast(jac_link_smear_type[i], 7, MPI_CHAR, 0, MPI_COMM_WORLD);
	
   MPI_Bcast(&num2pt,1,MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&mu2pt,num2pt,MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&num3pt,1,MPI_INT, 0, MPI_COMM_WORLD); 
   MPI_Bcast(&mu3pt,num3pt,MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&numb2pt,1,MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&mub2pt,numb2pt,MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&numb3pt,1,MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&mub3pt,numb3pt,MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&max_mom,1,MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&nmom_max,1,MPI_INT, 0, MPI_COMM_WORLD);	
	
   MPI_Bcast(&noprops,1,MPI_INT,0,MPI_COMM_WORLD); 
   MPI_Bcast(&smearno, noprops, MPI_INT, 0, MPI_COMM_WORLD);
   for (i=0;i<NMAX_PARM;i++)
     MPI_Bcast(src_type[i],7,MPI_CHAR,0,MPI_COMM_WORLD);

   MPI_Bcast(&source_sh_no,noprops,MPI_INT,0, MPI_COMM_WORLD);
   MPI_Bcast(&noextprops,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&kappa,noprops+noextprops,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&theta,4*(noprops+noextprops),MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&pos,4*noprops,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&t_snk,noextprops,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&propid,noprops,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&propidfirstleg,noextprops,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&extprop_sk_sm,noextprops,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&mom_insert,noextprops,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&gamma_insert,noextprops,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&notwopt,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nothreept,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&twoptprop1,notwopt,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&twoptprop2,notwopt,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&twopt_sk_sm,notwopt,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&threeptprop1,nothreept,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&threeptprop2,nothreept,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&threept_sk_sm,nothreept,MPI_INT,0,MPI_COMM_WORLD);

   /**** baryon block ****/
   MPI_Bcast(&noextprops_b,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&extpropid_b,noextprops_b,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&propidleg1_b,noextprops_b,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&propidleg2_b,noextprops_b,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&propidleg3_b,noextprops_b,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&extprop_b_sk_sm,noextprops_b,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&kappa_b,noextprops_b,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&theta_b,4*noextprops_b,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&t_snk_b,noextprops_b,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&mom_insert_b,noextprops_b,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&pol3pt,noextprops_b,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&notwopt_b,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nothreept_b,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&twoptid_b,notwopt_b,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&twoptprop1_b,notwopt_b,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&twoptprop2_b,notwopt_b,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&twoptprop3_b,notwopt_b,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&twopt_b_sk_sm,notwopt_b,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&pol2pt,notwopt_b,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&threeptprop1_b,nothreept_b,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&threeptprop2_b,nothreept_b,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&threept_b_sk_sm,nothreept_b,MPI_INT,0,MPI_COMM_WORLD);
}


static void setup_files(void)
{
   check_dir_root(log_dir);
   /*check_dir_root(cnfg_dir);*/
 /*  check_dir_root(corr_dir);
 */  check_dir_root(sfld_dir);
   
   sprintf(nbase,"%dx%dx%dx%db%1.2fk%.5fc%.5fid%d",
           NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3,beta,kappasea,csw,id);

   sprintf(ncnfg,"/.%d.log~",pid);   
   error_root((strlen(log_dir)+strlen(nbase)+strlen(ncnfg))>=NAME_SIZE,1,
              "setup_files [measure8.c]","log_dir name is too long");
      
   sprintf(log_file,"%s/%s.%d.log",log_dir,nbase,pid);
   sprintf(end_file,"%s/%s.%d.end",log_dir,nbase,pid);
   sprintf(log_save,"%s~",log_file);

   sprintf(ncnfg,"/n%d",nc_last);
   error_root((strlen(cnfg_dir)+strlen(ncnfg))>=NAME_SIZE,1,
              "setup_files [measure8.c]","cnfg_dir name is too long");   

   sprintf(ncnfg,"/.%dn%dc43",pid,nc_last);
   error_root((strlen(corr_dir)+strlen(ncnfg))>=NAME_SIZE,1,
              "setup_files [measure8.c]","corr_dir name is too long");   

   sprintf(ncnfg,"/.%dn%dc43",pid,nc_last);
   error_root((strlen(sfld_dir)+strlen(ncnfg))>=NAME_SIZE,1,
              "setup_files [measure8.c]","sfld_dir name is too long");   

   fend=fopen(log_file,"r");
   error_root(fend!=NULL,1,"check_files [measure8.c]",
              "Attempt to overwrite old *.log");
}


static void print_info(void)
{
   int ip,i;
   
   if (my_rank==0)
   {
      ip=ftell(flog);
      fclose(flog);

      if (ip==0)
         remove("MPI_Errors");
      
      flog=freopen(log_file,"w",stdout);
      error_root(flog==NULL,1,"print_info [measure8.c]",
		                            "Unable to open log file");
      
      printf("\n");
      printf("Computation of quark propagators\n");
      printf("--------------------------------\n\n");

      printf("Program version %s\n",DD_HMC_RELEASE);
      printf("with qcd-measure %s\n",QCD_MEASURE_RELEASE);

      if (endian==LITTLE_ENDIAN)
         printf("The machine is little endian\n");
      else
         printf("The machine is big endian\n");
         
      printf("\n");      
      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
      printf("prop_dump = %d\n",prop_dump);

      printf("O(a) improved theory\n");
      printf("beta = %.6f\n",lat.beta);
      printf("kappasea = %.6f\n",kappasea);
      printf("csw = %.6f\n",lat.csw);
      printf("bc = %+d\n\n",bc);
      printf("start = %s\n",start);
      printf("level = %d, seed = %d, effective seed = %d\n\n",
                level,seed,seed^nc_first);

      printf("bs_sap = %d %d %d %d\n",
             bs_sap[0],bs_sap[1],bs_sap[2],bs_sap[3]);
      printf("nkv_sap,nmr_sap,ncy_sap = %d,%d,%d\n",
             gcr.nkv,gcr.nmr,gcr.ncy);
      printf("nmx = %d\n",nmx);
      printf("res = %.1e\n\n",res);   

      printf("bs_dfl = %d %d %d %d\n",
             bs_dfl[0],bs_dfl[1],bs_dfl[2],bs_dfl[3]);
      printf("Ns = %d\n",dfl.Ns);
      printf("rkappa = %.6f\n",dfl.rkappa);
      printf("rinv = %d\n",dfl.rinv);
      printf("rnmr = %d\n",dfl.rnmr);
      printf("rncy = %d\n",dfl.rncy);         
      printf("rvnmx = %d\n",dfl.rvnmx);
      printf("rvres = %.1e\n",dfl.rvres);
      printf("vnkv = %d\n",dfl.vnkv);         
      printf("vnmx = %d\n",dfl.vnmx);
      printf("vres = %.1e\n",dfl.vres);
      printf("vdnmx = %d\n",dfl.vdnmx);
      printf("vdres = %.1e\n\n",dfl.vdres);
      printf("vdres = %.1e\n\n",dfl.vdres);
      printf("no_hits = %d\n",nohits);
      printf("dtsrc = %d\n",dtsrc);
	  
      printf("no_source_sh = %d\n", no_source_sh);
      for(i=0;i<no_source_sh;i++)
       printf("source_shift = %d %d %d %d %d\n", source_shift_id[i],
               source_shift[i][0],source_shift[i][1],source_shift[i][2],
               source_shift[i][3]);
	   
      printf("no_smear_parm = %d\n", nosmearparm);
      for(i=0;i<nosmearparm;i++)
	       printf("smearing = %d %s %.2e %.2e %.2e\n", smearid[i],
		                    smearing_type[i],alpha[i][0],alpha[i][1],alpha[i][2]);
	  
      printf("no_jac_parm = %d\n", no_jac_parm);
      for(i=0; i<no_jac_parm; i++)
	       printf("smearing = %d %s %d %.1e %d\n",jac_parm_id[i],jac_link_smear_type[i],
		                                        jac_N[i],jac_kappa[i],jac_link_smear_no[i]); 

	     printf("num2pt = %d\n",num2pt);
      printf("         ");
      for (i=0;i<num2pt;i++) printf("%d ",mu2pt[i]); 
      printf("\n");
      printf("num3pt = %d\n",num3pt);
      printf("         ");
      for (i=0;i<num3pt;i++) printf("%d ",mu3pt[i]);
      printf("\n");
      printf("numb2pt = %d\n",numb2pt);
      printf("         ");
      for (i=0;i<numb2pt;i++) printf("%d ",mub2pt[i]);
      printf("\n");
      printf("numb3pt = %d\n",numb3pt);
      printf("         ");
      for (i=0;i<numb3pt;i++) printf("%d ",mub3pt[i]);
      printf("\n");

      printf("max_mom = %d\n",max_mom);
      printf("nmom_max = %d\n",nmom_max);
      		
      printf("no_props = %d\n",noprops);
      for (i=0;i<noprops;i++)
        printf("prop = %d %.6f %s %d %d %d %d %.1e %.1e %.1e %d %d\n",propid[i],
		                kappa[i],src_type[i],pos[i][0],pos[i][1],pos[i][2],pos[i][3],
	  	              theta[i][1],theta[i][2],theta[i][3],smearno[i],source_sh_no[i]);
 
      printf("no_ext_props = %d\n",noextprops);
      for (i=0;i<noextprops;i++)
        printf("ext_prop = %d %d %d %.6f %d %d %d %.1e %.1e %.1e\n",extpropid[i],
		         propidfirstleg[i],extprop_sk_sm[i],kappa[i+noprops],t_snk[i],mom_insert[i],gamma_insert[i],
		                                       theta[i+noprops][1],theta[i+noprops][2],theta[i+noprops][3]);		

      printf("no_2pt = %d\n",notwopt);
      for(i=0;i<notwopt;i++)
        printf("2pt = %d %d %d %d\n",
		               twoptid[i],twoptprop1[i],twoptprop2[i], twopt_sk_sm[i]);

      printf("no_3pt = %d\n",nothreept);
      for(i=0;i<nothreept;i++)
        printf("3pt = %d %d %d %d\n",
		               threeptid[i],threeptprop1[i],threeptprop2[i],threept_sk_sm[i]);

      printf("no_ext_props_b = %d\n",noextprops_b);
      for(i=0;i<noextprops_b;i++)
        printf("ext_prop_b = %d %d %d %d %d %.6f %d %d %.1e %.1e %.1e %d\n",
		              extpropid_b[i],propidleg1_b[i],propidleg2_b[i],propidleg3_b[i],
		              extprop_b_sk_sm[i],kappa_b[i],t_snk_b[i], 
		              mom_insert_b[i],theta_b[i][1],theta_b[i][2],theta_b[i][3],pol3pt[i]);

      printf("no_2pt_b = %d\n",notwopt_b);
      for(i=0;i<notwopt_b;i++)
        printf("2pt_b = %d %d %d %d %d %d\n",twoptid_b[i],twoptprop1_b[i],
                    twoptprop2_b[i], twoptprop3_b[i], twopt_b_sk_sm[i],pol2pt[i]);

      printf("no_3pt_b = %d\n",nothreept_b);
      for (i=0;i<nothreept_b;i++)
       printf("3pt_b = %d %d %d %d\n",threeptid_b[i],threeptprop1_b[i],
                threeptprop2_b[i],threept_b_sk_sm[i]);

       printf("total number of sink smearings:%d\n", sk_sm_nb);

       printf("Configurations no %d -(delta=%d)-> %d\n",nc_first,nc_delta,nc_last); 

      fflush(flog);
   }
}


static void check_endflag(int *iend)
{
   if (my_rank==0)
   {
      fend=fopen(end_file,"r");

      if (fend!=NULL)
      {
         fclose(fend);
         remove(end_file);
         (*iend)=1;
         printf("\nEnd flag set, run stopped\n");
      }
   }

   MPI_Bcast(iend,1,MPI_INT,0,MPI_COMM_WORLD);
}

void static flush_global_cmplx_double(int mom, int num_c,complex_dble *k)
{
 int i;
 for(i=0;i<num_c*num_c*mom*NPROC0*L0;i++){
   (*(k+i)).re=0.0; 
   (*(k+i)).im=0.0; 
 }
}

void static sink_smearing_counter(int counter, int *array_sk_sm, int *a1, int *a2, int *a3){

   int i;

   for(i=0;i<counter;i++){
    switch(array_sk_sm[i]){
     case 0:
      no_sk_sm[0][a1[i]]++;
      no_sk_sm[0][a2[i]]++;
      no_sk_sm[0][a3[i]]++;
      break;
     case 1:
      no_sk_sm[1][a1[i]]++;
      no_sk_sm[1][a2[i]]++;
      no_sk_sm[1][a3[i]]++;
      break;
     case 2:
      no_sk_sm[2][a1[i]]++;
      no_sk_sm[2][a2[i]]++;
      no_sk_sm[2][a3[i]]++;
      break;
     case 3:
      no_sk_sm[3][a1[i]]++;
      no_sk_sm[3][a2[i]]++;
      no_sk_sm[3][a3[i]]++;
      break;
     case 4:
      no_sk_sm[4][a1[i]]++;
      no_sk_sm[4][a2[i]]++;
      no_sk_sm[4][a3[i]]++;
     case 5:
      no_sk_sm[5][a1[i]]++;
      no_sk_sm[5][a2[i]]++;
      no_sk_sm[5][a3[i]]++;
     case 6:
      no_sk_sm[6][a1[i]]++;
      no_sk_sm[6][a2[i]]++;
      no_sk_sm[6][a3[i]]++;
     case 7:
      no_sk_sm[7][a1[i]]++;
      no_sk_sm[7][a2[i]]++;
      no_sk_sm[7][a3[i]]++;
     case 8:
      no_sk_sm[8][a1[i]]++;
      no_sk_sm[8][a2[i]]++;
      no_sk_sm[8][a3[i]]++;
     case 9:
      no_sk_sm[9][a1[i]]++;
      no_sk_sm[9][a2[i]]++;
      no_sk_sm[9][a3[i]]++;
     case 10:
      no_sk_sm[10][a1[i]]++;
      no_sk_sm[10][a2[i]]++;
      no_sk_sm[10][a3[i]]++;
     case 11:
      no_sk_sm[11][a1[i]]++;
      no_sk_sm[11][a2[i]]++;
      no_sk_sm[11][a3[i]]++;
     case 12:
      no_sk_sm[12][a1[i]]++;
      no_sk_sm[12][a2[i]]++;
      no_sk_sm[12][a3[i]]++;
     case 13:
      no_sk_sm[13][a1[i]]++;
      no_sk_sm[13][a2[i]]++;
      no_sk_sm[13][a3[i]]++;
     case 14:
      no_sk_sm[14][a1[i]]++;
      no_sk_sm[14][a2[i]]++;
      no_sk_sm[14][a3[i]]++;
     case 15:
      no_sk_sm[15][a1[i]]++;
      no_sk_sm[15][a2[i]]++;
      no_sk_sm[15][a3[i]]++;
     case 16:
      no_sk_sm[16][a1[i]]++;
      no_sk_sm[16][a2[i]]++;
      no_sk_sm[16][a3[i]]++;
     case 17:
      no_sk_sm[17][a1[i]]++;
      no_sk_sm[17][a2[i]]++;
      no_sk_sm[17][a3[i]]++;
     case 18:
      no_sk_sm[18][a1[i]]++;
      no_sk_sm[18][a2[i]]++;
      no_sk_sm[18][a3[i]]++;
     case 19:
      no_sk_sm[19][a1[i]]++;
      no_sk_sm[19][a2[i]]++;
      no_sk_sm[19][a3[i]]++;
     case 20:
      no_sk_sm[20][a1[i]]++;
      no_sk_sm[20][a2[i]]++;
      no_sk_sm[20][a3[i]]++;
      break;
    }
   }
}


int main(int argc,char *argv[])
{
  int icnfg,stochnorm=1, idirac, icolour, L[4], shift[4];
  int iend,iw0,iw2,iw1,iw3;
  int pim[]={-1,-1,-1,-1}; /* "props in memory" look up table */
  int i,j,k,n,ihit,isv[4],ifail,status[3];
  int *momenta=NULL;
  double wt1,wt2,wt3,wt4,wtavg,avpl,wt5,wt6;
  float rnd;
  int rnd_pos[4];
  char out_file[NAME_SIZE];
  char out_fsv[NAME_SIZE];
  char in_prop[NAME_SIZE];
  complex_dble *sst2, *sst3; 
  complex_dble *dum2pt,*dum3pt; /* *dum,*/
  full_spinor_dble **sv,**fsvaux,**prop_pt;
  spinor_dble **wsd;
  srcdef src;
  smearparm parm;
  propinfo prop;
  complex_dble **s_pointer = NULL;
  complex_dble *s_pointer_3pt = NULL;
  /**************************************************************************/
  /* some preparatory stuff *************************************************/
  /**************************************************************************/
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  read_infile(argc,argv);
  setup_files();
  geometry();
  momenta   = amalloc(nmom_max*3*sizeof(int),3);
  nmom_max=gen_momenta(max_mom,momenta);

  /* sink smearing memory stuff */
  sk_sm_nb = 0;
  for(j=0;j<21;j++) for(i=0;i<NMAX_TWIST;i++) no_sk_sm[j][i] = 0;
   sink_smearing_counter(notwopt,twopt_sk_sm,twoptprop1,twoptprop2,twoptprop1);
   sink_smearing_counter(nothreept,threept_sk_sm,threeptprop1,
       	threeptprop1,threeptprop1);
   sink_smearing_counter(noextprops,extprop_sk_sm,propidfirstleg,
       	propidfirstleg,propidfirstleg);
   sink_smearing_counter(notwopt_b,twopt_b_sk_sm,twoptprop1_b,
       	twoptprop2_b,twoptprop3_b);
   sink_smearing_counter(nothreept_b,threept_b_sk_sm,threeptprop1_b,
       	threeptprop1_b,threeptprop1_b);
   sink_smearing_counter(noextprops_b,extprop_b_sk_sm,propidleg1_b,
       	propidleg2_b,propidleg3_b);
  for(j=0;j<21;j++) for(i=0;i<NMAX_TWIST;i++) if(no_sk_sm[j][i] != 0)sk_sm_nb++;

  /**************************************************************************/
  /* dynamical-memory allocation ********************************************/
  /**************************************************************************/
   dum2pt = amalloc(num2pt*num2pt*nmom_max*NPROC0*L0*sizeof(complex_dble),3);
   dum3pt = amalloc(num3pt*num3pt*nmom_max*NPROC0*L0*sizeof(complex_dble),3);

   sst2   = amalloc(num2pt*num2pt*nmom_max*NPROC0*L0*sizeof(complex_dble),3);
   sst3   = amalloc(num3pt*num3pt*nmom_max*NPROC0*L0*sizeof(complex_dble),3);
   print_info();
   flush_global_cmplx_double(nmom_max,num2pt,sst2);
   flush_global_cmplx_double(nmom_max,num3pt,sst3);

   if(notwopt_b > 0){
    s_pointer = amalloc(6*sizeof(complex_dble *),3);
    for(i=0; i<6; i++) *(s_pointer+i) = 
        amalloc(numb2pt*nmom_max*NPROC0*L0*sizeof(complex_dble),3);
   }
   if(nothreept_b > 0){
    s_pointer_3pt = amalloc(16*numb3pt*nmom_max*NPROC0*L0*sizeof(complex_dble),3);
   }


   /****************************************************************/
   /* allocate memory space for all   props if prop_dump = 0       */
   /* allocate memory space for three props if prop_dump = 1       */
   /* 			  keep props afterwards                    */
   /* allocate memory space for three props if prop_dump = 2       */
   /* 			delete props afterwards            	   */
   error_root((prop_dump!=0)&&(prop_dump!=1)&&(prop_dump!=2),
        1,"main [measure8.c]",
       "prop_dump parameter is out of range, should be 0, 1 or 2 but is %d",
        prop_dump);
   prop_pt = amalloc(4*sizeof(full_spinor_dble *),3); 
   for(i=0;i<4;i++) prop_pt[i] = amalloc(sizeof(full_spinor_dble),3);
   sv = amalloc((noprops+noextprops+noextprops_b+sk_sm_nb)*
        		sizeof(full_spinor_dble *),3);

   for(i=0;i<noprops+noextprops+noextprops_b+sk_sm_nb;i++){
    if (prop_dump==0)    sv[i]=amalloc(VOLUME*(sizeof(full_spinor_dble)),3);
    else if(prop_dump>0) sv[i]=amalloc((sizeof(full_spinor_dble)),3);
   }
   if(prop_dump>0){ /* allocate auxiliary fields */
    fsvaux   = amalloc(4*sizeof(full_spinor_dble *),3); 
    for(i=0;i<4;i++) fsvaux[i] = amalloc(VOLUME*sizeof(full_spinor_dble),3);
   } 

   alloc_ud();
   alloc_swd();
   alloc_u();
   alloc_sw();
   if (Ns<(2*nkv_sap))
      alloc_s(2*nkv_sap+1);
   else
   alloc_s(Ns+2);
   alloc_sd(3+4+1+2); /* dfl_sap_gcr + correlator code + residue */
   wsd=reserve_wsd(4);
   iw0=(wsd[0]-psd[0][0])/NSPIN;
   iw1=(wsd[1]-psd[0][0])/NSPIN;
   iw2=(wsd[2]-psd[0][0])/NSPIN;
   iw3=(wsd[3]-psd[0][0])/NSPIN; 

   /**************************************************************************/
   /* other preparatory stuff ************************************************/ 
   /**************************************************************************/
   L[0]=L0*NPROC0;
   L[1]=L1*NPROC1;
   L[2]=L2*NPROC2;
   L[3]=L3*NPROC3;
   init_coords();
   gcr=set_gcr_parms(bs_sap,nmr_sap,ncy_sap,nkv_sap);
   lat=set_lat_parms(beta,kappa[0],csw);
   dfl=set_dfl_parms(bs_dfl,Ns,rinv,rnmr,rncy,rvnmx,vnkv,vnmx,vdnmx,
                     rkappa,rvres,vres,vdres,0,0);
   print_info();
   iend=0;
   wtavg=0.0;
   start_ranlux(level,seed^nc_first);

   message("prop_dump=%d\n", prop_dump); fflush(stdout);


  /**************************************************************************/
  /**************************************************************************/
  /* start actual computation ***********************************************/
  /**************************************************************************/
  /**************************************************************************/
  message("number of momentum modes: %d\n",nmom_max);
  fflush(stdout);

  for (icnfg=nc_first;(iend==0)&&(icnfg<=nc_last);icnfg+=nc_delta)
  {

   /********************* sink smearing lookuptable *************************/
   sk_sm_nb = 0;
   for(j=0;j<21;j++) for(i=0;i<NMAX_TWIST;i++) no_sk_sm[j][i] = 0;
   sink_smearing_counter(notwopt,twopt_sk_sm,twoptprop1,twoptprop2,twoptprop1);
   sink_smearing_counter(nothreept,threept_sk_sm,threeptprop1,
       			 threeptprop1,threeptprop1);
   sink_smearing_counter(noextprops,extprop_sk_sm,propidfirstleg,
       			 propidfirstleg,propidfirstleg);
   sink_smearing_counter(notwopt_b,twopt_b_sk_sm,twoptprop1_b,
       			 twoptprop2_b,twoptprop3_b);
   sink_smearing_counter(nothreept_b,threept_b_sk_sm,threeptprop1_b,
       			 threeptprop1_b,threeptprop1_b);
   sink_smearing_counter(noextprops_b,extprop_b_sk_sm,propidleg1_b,
       			 propidleg2_b,propidleg3_b);
   for(j=0;j<21;j++) for(i=0;i<NMAX_TWIST;i++) if(no_sk_sm[j][i] != 0)sk_sm_nb++;
   /*************************************************************************/  

   sprintf(in_prop,"%s/prop_%s.%dn%dprop",sfld_dir,nbase,pid,icnfg);
   MPI_Barrier(MPI_COMM_WORLD);
   wt1=MPI_Wtime();      
   message("\nConfiguration no %d:\n",icnfg);
   sprintf(cnfg_file,"%s/%sn%d",cnfg_dir,nbase,icnfg);
   message("Config name: %s\n",cnfg_file);
   fflush(stdout);

   /* source shift for all propagators *************************************/
    for(i=0;i<no_source_sh;i++){
      message("source_sh_no: %d; source_shif: t=%d, x=%d, y=%d, z=%d\n",
               i, source_shift[i][0], source_shift[i][1], source_shift[i][2], source_shift[i][3]);
    }

    for(j=0;j<4;j++){
    if(my_rank == 0) ranlxs(&rnd, 1);
     MPI_Bcast(&rnd, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
     rnd_pos[j] = (int) (rnd*L[j]);
    }
    for(i=0;i<noprops;i++){
      for(j=0;j<4;j++){
        if((source_shift[source_sh_no[i]][j] != 999) && ((icnfg-nc_first) == 0) )
              pos[i][j] = (pos[i][j] + source_shift[source_sh_no[i]][j])/*%L[j]*/;
        else if(source_shift[source_sh_no[i]][j] == 999)
                pos[i][j] = rnd_pos[j];
      }
    }

    for(i=0;i<noprops;i++){
      message("prop: %d; source position: t=%d, x=%d, y=%d, z=%d\n", 
               i, pos[i][0], pos[i][1], pos[i][2], pos[i][3]);
    }

   start_ranlux(level,seed^nc_first);
   /* load and check and assign gauge field *******************************/
    if((strcmp(start,"config"))==0){
     message("loading configuration\n");
     fflush(stdout);
     import_cnfg(cnfg_file);
    }
    else{
     message("running free field\n");
     fflush(stdout);
    }

    /*******************************/
    /* random gauge transformation */
/*    transform_ud(pud);
*/
    avpl=plaq_sum_dble()/(double)(6*VOLUME*NPROC);
    message("Average plaquette %e\n", avpl);  
    fflush(stdout);    

    if (bc==-1){
      flipbcd();
      message("flipped bcs\n");
      fflush(stdout);
    }
    /* initialize the correlator export to LIME file */
    sprintf(out_file,"%s/correlator_%s.%dn%d",corr_dir,nbase,pid,icnfg);
    if(my_rank==0) export_init(out_file,argc,argv);
    assign_ud2u();
    assign_swd2sw();
    assign_ud2ubgr(GCR_BLOCKS);
   /* assign some variables and flush correlator arrays ********************/
    prop.csw   	= csw;
    prop.beta  	= beta; 
   /************************************************************************/
   /* propagator computation ***********************************************/
   /************************************************************************/
   for (ihit=0;ihit<nohits;ihit++){
    message("\nThis is hit number %d\n\n",ihit);
    fflush(stdout);
    set_sd2zero(VOLUME,wsd[0]);
    random_Z4(VOLUME,wsd[0]);
    /***********************************************************************/
    /* generate standard propagators **************************************/ 
    for (i=0;i<noprops;i++){
      if((strcmp(src_type[i],"Z4spin")==0)||(strcmp(src_type[i],"Z4")==0)||
      	(ihit==0)){
       /* define source and prop specifications */
       strcpy((src.type),(src_type[i]));
       k = jac_link_smear_no[smearno[i]];
       strcpy((parm.type),(smearing_type[k]));
       parm.smearid  	=  smearid[k];
       parm.alpha[0] 	=  alpha[k][0];
       parm.alpha[1] 	=  alpha[k][1];
       parm.alpha[2] 	=  alpha[k][2];
       src.pos[0]	=  (pos[i][0]/*+ihit*dtsrc*/)/*%L[0]*/;
       src.pos[1]	=  pos[i][1];
       src.pos[2]	=  pos[i][2];
       src.pos[3]	=  pos[i][3];
       src.n		=  jac_N[smearno[i]];
       src.kappa	=  jac_kappa[smearno[i]];
       src.smear	=  parm;
       prop.kappa 	=  kappa[i];
       prop.theta[0]  	=  theta[i][0]; 
       prop.theta[1]  	=  theta[i][1]; 
       prop.theta[2]  	=  theta[i][2]; 
       prop.theta[3]  	=  theta[i][3]; 
       prop.src   	=  src;
       prop.nmx		=  nmx;
       prop.res		=  res;

       /* shifting */ 
       shift[0] = -prop.src.pos[0];
       shift[1] = -prop.src.pos[1];
       shift[2] = -prop.src.pos[2];
       shift[3] = -prop.src.pos[3];
       shift_ud(shift);
       sw_term();
       copy_bnd_ud();
       message("\n\nShifted gauge field by (%d,%d,%d,%d)\n",
          shift[0], shift[1], shift[2], shift[3]);

       /* creating the deflation modes if necessary */
       if( (i==0) || (kappa[i]!=kappa[i-1]) ||
           (pos[i][0]!=pos[i-1][0]) || (pos[i][1]!=pos[i-1][1]) ||
           (pos[i][2]!=pos[i-1][2]) || (pos[i][3]!=pos[i-1][3])){
         lat=set_lat_parms(prop.beta,prop.kappa,prop.csw);
         wt5=MPI_Wtime();
         dfl_modes(0,status,&ifail);
         wt6=MPI_Wtime();
         message("Deflation subspace generation: status = %d, ifail = %d, t = %.1lf s\n",
           status[0],ifail,wt6-wt5);
       }
       /* invert the thing */
       if (prop_dump==0) prop_pt[0]=sv[i];
       else              prop_pt[0]=fsvaux[0];
       propagator(prop_pt[0],wsd[0],prop);
       message("\n\n"); fflush(stdout);
       sprintf(out_fsv,"%s/prop_%s.%dn%dprop%d",sfld_dir,nbase,pid,icnfg,i);
       if (prop_dump>0) write_fsv(out_fsv,prop_pt[0]);
       /* shifting back */
       shift[0] = -shift[0];
       shift[1] = -shift[1];
       shift[2] = -shift[2];
       shift[3] = -shift[3];
       shift_ud(shift);
      }
    }
    /* end standard propagators ********************************************/
    /***********************************************************************/

    /***********************************************************************/
    /* start sink-smearing currently limited to just 5 smearing parm sets **/
    n = 0;
    for(i=0;i<noprops;i++){ 
     for(j=0;j<21;j++){ 
      if(no_sk_sm[j][i] != 0){
       isv[0]=i;
       read_fsv2(in_prop,prop_dump,&pim[0],1,&isv[0],prop_pt,sv,fsvaux);
       if(prop_dump>0) prop_pt[1]=fsvaux[1];
       else 	       prop_pt[1]=sv[noprops+noextprops+noextprops_b+n];
       message("\nSink smearing no %d of propagator no %d: \n",j,i);
       wt3 = MPI_Wtime();
       k = jac_link_smear_no[j];
       strcpy((src.type),(jac_link_smear_type[j]));
       strcpy((parm.type),(smearing_type[k]));
       parm.smearid  	=  smearid[k];
       parm.alpha[0] 	=  alpha[k][0];
       parm.alpha[1] 	=  alpha[k][1];
       parm.alpha[2] 	=  alpha[k][2];
       src.pos[0]	= (pos[i][0]+ihit*dtsrc)%L[0];
       src.pos[1]	=  pos[i][1];
       src.pos[2]	=  pos[i][2];
       src.pos[3]	=  pos[i][3];
       src.n		=  jac_N[j];
       src.kappa	=  jac_kappa[j];
       src.smear	=  parm;

       shift[0] = -src.pos[0];
       shift[1] = -src.pos[1];
       shift[2] = -src.pos[2];
       shift[3] = -src.pos[3];
       shift_ud(shift);
       copy_bnd_ud();
       message("gauge field shifted by (%d, %d, %d, %d)\n",
           -src.pos[0], -src.pos[1], -src.pos[2], -src.pos[3] );
        
       alloc_pud_sm1();
        
       /*defines the way of linksmearing*/
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

       if(my_rank == 0) {
         message("\nJacobi-smearing with n=%d and kappa=%.2f\n", src.n, src.kappa);
         if(strcmp(src.type, "J_APE") == 0)
           message("and APE-smearing with alpha = %.2f\n", src.smear.alpha[0]);
         else if (strcmp(src.type, "J_HYP") == 0)
           message("and HYP-smearing with alpha1 = %.2f, alpha2 = %.2f and alpha3 = %.2f \n",
                                  src.smear.alpha[0], src.smear.alpha[1], src.smear.alpha[2]);
       }

       for(idirac=1; idirac<=4; idirac++){
        for(icolour=1; icolour<=3; icolour++){
         copy_fs_sv(VOLUME, prop_pt[0], wsd[0], idirac, icolour);
         jacobi_sink(wsd[0], src, idirac, icolour);
         copy_sv_fs(VOLUME, wsd[0], prop_pt[1], idirac, icolour);
        }
       }
       if (prop_dump>0){/* Write smeared prop to disk */
        sprintf(out_fsv,"%s/prop_%s.%dn%dprop%d",sfld_dir,nbase,pid,icnfg,
		      	noprops+noextprops+noextprops_b+n);
        write_fsv(out_fsv,prop_pt[1]);
       }
       wt4 = MPI_Wtime();
       message("Time for sink smearing %d: %.2e\n",n,wt4-wt3);
       no_sk_sm[j][i] = n;
       n++;
       shift[0] = src.pos[0];
       shift[1] = src.pos[1];
       shift[2] = src.pos[2];
       shift[3] = src.pos[3];
       shift_ud(shift);
       free_pud_sm1();
      }
     }
    }
    /* end sink-smearing ***************************************************/
    /***********************************************************************/

    /***********************************************************************/
    /* generate sequential source propagators ******************************/ 
    for(i=noprops;i<noextprops+noprops;i++){
     if((strcmp(src_type[propidfirstleg[i-noprops]],"Z4spin")==0)||
      	(strcmp(src_type[propidfirstleg[i-noprops]],"Z4")==0)||(ihit==0))
     {
      /* define source and prop specifications */
      src.pos[0]	= (pos[propidfirstleg[i-noprops]][0]+ihit*dtsrc)%L[0];
      src.pos[1]	=  pos[propidfirstleg[i-noprops]][1];
      src.pos[2]	=  pos[propidfirstleg[i-noprops]][2];
      src.pos[3]	=  pos[propidfirstleg[i-noprops]][3];
      prop.t_snk	=  t_snk[i-noprops];
      prop.kappa 	=  kappa[i];
      prop.theta[0]  	=  theta[i][0]; 
      prop.theta[1]  	=  theta[i][1]; 
      prop.theta[2]  	=  theta[i][2]; 
      prop.theta[3]  	=  theta[i][3]; 
      prop.src   	=  src;
      prop.nmx		=  nmx;
      prop.res		=  res;
      prop.mom_ins	=  momenta+mom_insert[i-noprops]*3;
      prop.gam_ins	=  gamma_insert[i-noprops];
      if(extprop_sk_sm[i-noprops] >= 0){
       k = noprops+noextprops+noextprops_b;
       message(strcat("Extended source with sinksmeared propagator:%d",
      		" with parameters: %s %d %lf\n"),
       		propidfirstleg[i-noprops],
      		smearing_type[extprop_sk_sm[i-noprops]],
       		jac_N[extprop_sk_sm[i-noprops]],
      		jac_kappa[extprop_sk_sm[i-noprops]]);
       isv[0] = k+no_sk_sm[extprop_sk_sm[i-noprops]][propidfirstleg[i-noprops]];
      }
      else{
       isv[0] = propidfirstleg[i-noprops];
      }
      read_fsv2(in_prop,prop_dump,&pim[0],1,&isv[0],prop_pt,sv,fsvaux);
      if (prop_dump==0) prop_pt[1]=sv[i];
      else              prop_pt[1]=fsvaux[1];
      seq_propagator(prop_pt[1], prop_pt[0], prop);
      seq_src_test(prop_pt[1],5,&pos[propidfirstleg[i-noprops]][0]);
      sprintf(out_fsv,"%s/prop_%s.%dn%dprop%d",
             sfld_dir,nbase,pid,icnfg,i);
      if (prop_dump>0) write_fsv(out_fsv,prop_pt[1]);
     }
    }
    /* end sequential source propagators ***********************************/ 
    /***********************************************************************/

    /***********************************************************************/
    /* generate the extended propagators for the baryons *******************/
    for(j=0;j<noextprops_b;j++){
     if(j==0) message("\ncalculation of extended propagator for the nucleon:\n\n");
     i = propidleg1_b[j];
     strcpy((src.type),(src_type[i]));
     k = jac_link_smear_no[smearno[i]];
     strcpy((parm.type),(smearing_type[k]));

     parm.smearid  	 =  smearid[k];
     parm.alpha[0] 	 =  alpha[k][0];
     parm.alpha[1] 	 =  alpha[k][1];
     parm.alpha[2] 	 =  alpha[k][2];
     src.n		 =  jac_N[smearno[i]];
     src.kappa	         =  jac_kappa[smearno[i]];
     src.smear	         =  parm;
     prop.src   	 =  src;
     prop.src.pos[0]     = (pos[propidleg1_b[j]][0]+ihit*dtsrc)%L[0];
     prop.src.pos[1]     =  pos[propidleg1_b[j]][1];
     prop.src.pos[2]     =  pos[propidleg1_b[j]][2];
     prop.src.pos[3]     =  pos[propidleg1_b[j]][3];
     prop.theta[0]  =  theta_b[j][0];
     prop.theta[1]  =  theta_b[j][1];
     prop.theta[2]  =  theta_b[j][2];
     prop.theta[3]  =  theta_b[j][3];
     prop.kappa     =  kappa_b[j];
     prop.nmx       =  nmx;
     prop.res       =  res;
     prop.t_snk     =  t_snk_b[j];
     prop.mom_ins   =  3*mom_insert_b[j]+momenta;

     /* shifting the gauge fields */
     shift[0] = -prop.src.pos[0];
     shift[1] = -prop.src.pos[1];
     shift[2] = -prop.src.pos[2];
     shift[3] = -prop.src.pos[3];
     shift_ud(shift);
     sw_term();
     copy_bnd_ud();   
     message("\n\ngauge field shifted by (%d, %d, %d, %d)\n",
           -prop.src.pos[0], -prop.src.pos[1],
           -prop.src.pos[2], -prop.src.pos[3] );

     /* creating the deflation modes of necessary */
     if(  (j==0) || (kappa[j]!=kappa[j-1]) ||
          (pos[propidleg1_b[j]][0]!=pos[propidleg1_b[j]-1][0]) || 
          (pos[propidleg1_b[j]][1]!=pos[propidleg1_b[j]-1][1]) ||
          (pos[propidleg1_b[j]][2]!=pos[propidleg1_b[j]-1][2]) || 
          (pos[propidleg1_b[j]][3]!=pos[propidleg1_b[j]-1][3])){
       lat=set_lat_parms(prop.beta,prop.kappa,prop.csw);
       wt5=MPI_Wtime();
       dfl_modes(0,status,&ifail);
       wt6=MPI_Wtime();
       message("Deflation subspace generation: status = %d, ifail = %d, t = %.1lf s\n",
           status[0],ifail,wt6-wt5);
     }

     if(extprop_b_sk_sm[j] >= 0){
       k = noprops+noextprops+noextprops_b;
       isv[0]=k+no_sk_sm[extprop_b_sk_sm[j]][propidleg1_b[j]];
       isv[1]=k+no_sk_sm[extprop_b_sk_sm[j]][propidleg2_b[j]];
       isv[2]=k+no_sk_sm[extprop_b_sk_sm[j]][propidleg3_b[j]];
       message("\nUsing sink-smeared propagators");
       message(" no %d and %d with parameters: %s %d %.2lf\n",
             isv[0], isv[1],src.type,src.n,src.kappa);
     }
     else{
       isv[0]=propidleg1_b[j];
       isv[1]=propidleg2_b[j];
       isv[2]=propidleg3_b[j];
     }
     isv[3]=noextprops+noprops+j;
     message("seq isv %d %d %d %d\n",isv[0],isv[1],isv[2],isv[3]);fflush(stdout);
     read_fsv2(in_prop,prop_dump,&pim[0],3,&isv[0],prop_pt,sv,fsvaux);
     if (prop_dump==0) prop_pt[3]=sv[isv[3]];
     else              prop_pt[3]=fsvaux[3];
     seq_propagator_b(prop_pt[3],prop_pt[0],prop_pt[1],prop_pt[2],prop,extprop_b_sk_sm[j],pol3pt[j]);

     sprintf(out_fsv,"%s/prop_%s.%dn%dprop%d",
             sfld_dir,nbase,pid,icnfg,j+noprops+noextprops);
     if (prop_dump>0) write_fsv(out_fsv,prop_pt[3]);

     /* shifting the gauge fields back */
     shift[0] = +prop.src.pos[0];
     shift[1] = +prop.src.pos[1];
     shift[2] = +prop.src.pos[2];
     shift[3] = +prop.src.pos[3];
     shift_ud(shift);
    }
    /* end the extended propagators for the baryons ************************/
    /***********************************************************************/

    /***********************************************************************/
    /* do 2pt and 3pt contractions for mesons ******************************/ 
    for(i=0;i<notwopt;i++){
     fflush(stdout);
     if(i==0){
       message("\n\nDoing mesonic two-point contractions\n");
       message("-----------------------------------------------------\n");
     }
     if((strcmp(src_type[twoptprop1[i]],"Z4spin")==0)||
        (strcmp(src_type[twoptprop1[i]],"Z4")==0)||(ihit==0)){
      flush_global_cmplx_double(nmom_max,num2pt,sst2);
      if(strcmp(src_type[twoptprop1[i]],"Z4spin")==0) stochnorm=nohits;
      else 				  stochnorm=1;
      flush_global_cmplx_double(nmom_max,num2pt,dum2pt);
      if(twopt_sk_sm[i] >= 0){
       k = noprops+noextprops+noextprops_b;
       isv[0]=k+no_sk_sm[twopt_sk_sm[i]][twoptprop1[i]];
       isv[1]=k+no_sk_sm[twopt_sk_sm[i]][twoptprop2[i]];
      }
      else{
       isv[0]=twoptprop1[i];
       isv[1]=twoptprop2[i];
      }
      read_fsv2(in_prop,prop_dump,&pim[0],2,&isv[0],prop_pt,sv,fsvaux);
      fflush(stdout);
      C2(prop_pt[0],prop_pt[1],num2pt,mu2pt,nmom_max,momenta,dum2pt);
      for(j=0;j<num2pt*num2pt*nmom_max*NPROC0*L0;j++){
       (*(sst2+j)).re += (*(dum2pt+j)).re/stochnorm;
       (*(sst2+j)).im += (*(dum2pt+j)).im/stochnorm;
      }
      meson_IOloop("2pt",out_file,sst2,nmom_max,momenta,&mu2pt[0],num2pt);
     }
    }

    for(i=0;i<nothreept;i++){
     fflush(stdout);
     if(i==0){
       message("\n\nDoing mesonic three-point contractions\n");
       message("-----------------------------------------------------\n");
     }
     if((strcmp(src_type[propidfirstleg[i]],"Z4spin")==0)||
       	(strcmp(src_type[propidfirstleg[i]],"Z4")==0)||(ihit==0)){
      flush_global_cmplx_double(nmom_max,num3pt,sst3);
      if(strcmp(src_type[propidfirstleg[i]],"Z4spin")==0) stochnorm=nohits;
      else				  stochnorm=1; 
      flush_global_cmplx_double(nmom_max,num3pt,dum3pt);
      if(threept_sk_sm[i] >= 0){
       k = noprops+noextprops+noextprops_b;
       isv[0]=k+no_sk_sm[twopt_sk_sm[i]][threeptprop1[i]];
       isv[1]=threeptprop2[i]+noprops;
      }
      else{
       isv[0]=threeptprop1[i];
       isv[1]=threeptprop2[i]+noprops;
      }
      read_fsv2(in_prop,prop_dump,&pim[0],2,&isv[0],prop_pt,sv,fsvaux);
      C2(prop_pt[0],prop_pt[1],num3pt,mu3pt,nmom_max,momenta,dum3pt);
      for(j=0;j<num3pt*num3pt*nmom_max*NPROC0*L0;j++){
       (*(sst3+j)).re += (*(dum3pt+j)).re/stochnorm;
       (*(sst3+j)).im += (*(dum3pt+j)).im/stochnorm;
      }
      meson_IOloop("3pt",out_file,sst3,nmom_max,momenta,&mu3pt[0],num3pt);
     }
    } 
    /* end 2pt and 3pt contractions for mesons *****************************/ 
    /***********************************************************************/
    }/* end loop over ihit */

    /***********************************************************************/
    /* begin 2pt and 3pt contractions for baryons***************************/
    for(j=0;j<notwopt_b;j++){
      fflush(stdout);
      if(j==0){
        message("\n\nDoing baryonic two-point contractions:\n");
        message("-----------------------------------------------------\n");
      }
      if(twopt_b_sk_sm[j] >= 0){
        k = noprops+noextprops+noextprops_b;	
        isv[0]=k+no_sk_sm[twopt_b_sk_sm[j]][twoptprop1_b[j]];
        isv[1]=k+no_sk_sm[twopt_b_sk_sm[j]][twoptprop2_b[j]];
        isv[2]=k+no_sk_sm[twopt_b_sk_sm[j]][twoptprop3_b[j]];
      }
      else{
        isv[0]=twoptprop1_b[j];
        isv[1]=twoptprop2_b[j];
        isv[2]=twoptprop3_b[j];
      }
     read_fsv2(in_prop,prop_dump,&pim[0],3,&isv[0],prop_pt,sv,fsvaux);

     message("seq isv %d %d %d\n",isv[0],isv[1],isv[2]);fflush(stdout);

     C2_b(prop_pt[0],prop_pt[1],prop_pt[2],
          numb2pt,mub2pt,nmom_max,momenta,s_pointer,pol2pt[j]);

     baryon_IOloop("2pt-uds",out_file,*(s_pointer+3),
                    nmom_max,momenta,&mub2pt[0],numb2pt);
     baryon_IOloop("2pt-uud",out_file,*(s_pointer+4),
                    nmom_max,momenta,&mub2pt[0],numb2pt);
     baryon_IOloop("2pt-uuu",out_file,*(s_pointer+5),
                    nmom_max,momenta,&mub2pt[0],numb2pt);
     
    }

    for(j=0;j<nothreept_b;j++){
      fflush(stdout);
      if(j==0){
        message("\n\nDoing baryonic three-point contractions:\n");
        message("-----------------------------------------------------\n");
      }
      prop.mom_ins   =  3*mom_insert_b[threeptprop2_b[j]]+momenta;
      isv[0] = threeptprop1_b[j];
      isv[1] = noextprops+noprops+threeptprop2_b[j];
      /* shifting the gauge fields for the diplaced operators */
      shift[0] = -pos[threeptprop1_b[j]][0];
      shift[1] = -pos[threeptprop1_b[j]][1];
      shift[2] = -pos[threeptprop1_b[j]][2];
      shift[3] = -pos[threeptprop1_b[j]][3];
      shift_ud(shift);
      copy_bnd_ud();

      message("seq isv %d %d\tshifted gauge field by (%d,%d,%d,%d)",
               isv[0],isv[1],shift[0],shift[1],shift[2],shift[3]);
      message("\tfinal momentum = (%d, %d, %d)\n ",*(prop.mom_ins),*(prop.mom_ins+1),*(prop.mom_ins+2));
      fflush(stdout);
      read_fsv2(in_prop,prop_dump,&pim[0],2,&isv[0],prop_pt,sv,fsvaux);

      C3_b(prop_pt[0],prop_pt[1],numb3pt,mub3pt,nmom_max,momenta,s_pointer_3pt,prop);

      baryon_IOloop("3pt-uud",out_file, s_pointer_3pt,
                     nmom_max,momenta,&mub3pt[0],16*numb3pt);
      /* shifting back */
      shift[0] = pos[threeptprop1_b[j]][0];
      shift[1] = pos[threeptprop1_b[j]][1];
      shift[2] = pos[threeptprop1_b[j]][2];
      shift[3] = pos[threeptprop1_b[j]][3];
      shift_ud(shift);             
      copy_bnd_ud();
    }
    /* end 2pt and 3pt contractions for baryons****************************/
    /**********************************************************************/


   /**********************************************************************/
   /* finalize ***********************************************************/ 
   /**********************************************************************/
   MPI_Barrier(MPI_COMM_WORLD);
   wt2=MPI_Wtime();
   message("\n\nConfiguration fully processed in %.2e sec\n",wt2-wt1);
   message("----------------------------------------------------------------\n\n");
   wtavg+=(wt2-wt1);
   check_endflag(&iend);
   error_chk();
   if (my_rank==0)
   {
      if ((iend==1)||(icnfg==nc_last))
      {
         printf("\n\n");
         printf("Processed configurations no %d -> %d\n",nc_first,icnfg);
         printf("Average time per configuration = %.2e sec\n\n",
                wtavg/(double)(icnfg-nc_first+1));
      }
      fflush(flog);
      copy_file(log_file,log_save);
   }
   /* remove all node-dumped propagators from disk */
   if(prop_dump==2){
    for (i=0;i<noprops+noextprops+noextprops_b;i++){
       sprintf(out_fsv,"%s/prop_%s.%dn%dprop%d_rank%d",
        	sfld_dir,nbase,pid,icnfg,i,my_rank);
       ifail=remove(out_fsv);
       error_root(ifail!=0,1,"main [measure8.c]",
              "Could not delete prop %s",out_fsv);
    }
   }
  }
  if (my_rank==0) fclose(flog);
  MPI_Finalize();
  exit(0);
}

