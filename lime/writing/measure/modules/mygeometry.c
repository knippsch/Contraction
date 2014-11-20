
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "misc.h"
#include "mpi.h"
#include "start.h"
#include "version.h"
#include "global.h"
#include "measure.h"

/* initialize the field which gives for each lexicographical index the 4 coordinates 
 * on the lattice
 */
void init_coords()
{
    int x[4];
    int ip, ix;
    int my_rank;

    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    for (x[0]=0;x[0]<L0*NPROC0;x[0]++)
    for (x[1]=0;x[1]<L1*NPROC1;x[1]++)
    for (x[2]=0;x[2]<L2*NPROC2;x[2]++)
    for (x[3]=0;x[3]<L3*NPROC3;x[3]++)
    {
	 ipt_global(x,&ip,&ix);
	 if (ip==my_rank) {
	     coords[ix].t=x[0];
	     coords[ix].x=x[1];
	     coords[ix].y=x[2];
	     coords[ix].z=x[3];
	 }
    }
}
