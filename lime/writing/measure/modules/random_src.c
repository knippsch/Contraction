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
#include "version.h"
#include "global.h"
#include "mygeometry.h"

#if 0
void  src_field_random(int vol, spinor_dble *pk, int t)
{
    spinor_dble *rpk;
    spinor_dble r;
    double sigma=1;
    int ix=0;

    set_sd2zero(vol,pk);
    for (rpk=pk;rpk<(pk+vol);rpk++)
    {
	if (coords[ix].t==t)
	{
	    gauss_dble((double*)(&r),24);
	    *rpk=r;
	}
	ix++;
    }
}


void  src_field_random(int vol, spinor_dble *pk, int t)
{
    spinor_dble *rpk;
    double r[24];
    double sigma=1;
    int ix=0, i;

    set_sd2zero(vol,pk);
    for (rpk=pk;rpk<(pk+vol);rpk++)
    {
	if (coords[ix].t==t)
	{
	    ranlxd(r,24);
	    for (i=0;i<24;i++) r[i]=(r[i]<0.5?-1:1);
	    *rpk=*((spinor_dble*)r);
	}
	ix++;
    }
}


#endif

/* randome U(1) noise on timeslice t */
void  src_field_random(int vol, spinor_dble *pk, int t)
{
    spinor_dble *rpk;
    complex_dble *cpt;
    double r[12];
    int ix=0, i;
    double twopi=8*atan(1.);

    set_sd2zero(vol,pk);
    for (rpk=pk;rpk<(pk+vol);rpk++)
    {
	if (coords[ix].t==t)
	{
	    ranlxd(r,12);
	    cpt=(complex_dble*)rpk;
	    for (i=0;i<12;i++) 
	    {
		cpt->re=cos(twopi*r[i]);
		cpt->im=sin(twopi*r[i]);
		cpt++;
	    }
	}
	ix++;
    }
}

/* randome U(1) noise on timeslice t, only one Dirac component */
void  src_field_dt_random(int vol, spinor_dble *pk, int d,  int t)
{
    spinor_dble *rpk;
    complex_dble *cpt;
    double r[12];
    int ix=0, i;
    double twopi=8*atan(1.);

    set_sd2zero(vol,pk);
    for (rpk=pk;rpk<(pk+vol);rpk++)
    {
	if (coords[ix].t==t)
	{
	    ranlxd(r,3);
	    cpt=(complex_dble*)(&rpk->c1+d);
	    for (i=0;i<3;i++) 
	    {
		cpt->re=cos(twopi*r[i]);
		cpt->im=sin(twopi*r[i]);
		cpt++;
	    }
	}
	ix++;
    }
}




