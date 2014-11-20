#include <math.h>
/*#include "nr.h"*/
#include "nrutil.h"

#define SWAP(a,b) {double tempre=(a.re);\
	           double tempim=(a.im);\
			  (a.re)=(b.re);\
			  (a.im)=(b.im);\
			  (b.re)=tempre;\
			  (b.im)=tempim}

static double cabs(complex_dble c){
 return sqrt(c.re*c.re+c.im*c.im);
}
static complex_dble cinv(complex_dble c){
 complex_dble d;
 d.re=c.re;
 d.im=c.im;
 c.re =  d.re/(d.re*d.re+d.im*d.im);
 c.im = -d.im/(d.re*d.re+d.im*d.im);
 return c;
}
static complex_dble cmul(complex_dble c1,complex_dble c2){
 c.re = c1.re*c2.re+c1.im*c2.im;
 c.im = c1.re*c2.im-c1.im*c2.re;
 return c;
}
void gaussj(complex_dble **a, int n,complex_dble **b,int m)
{
        int *indxc,*indxr,*ipiv;
        int i,icol,irow,j,k,l,ll;
        double big;
        complex_dble dum,pivinv,c;

        indxc=ivector(1,n);
        indxr=ivector(1,n);
        ipiv=ivector(1,n);
        for (j=1;j<=n;j++) ipiv[j]=0;
        for (i=1;i<=n;i++) {
                big=0.0;
                for (j=1;j<=n;j++)
                        if (ipiv[j] != 1)
                                for (k=1;k<=n;k++) {
                                        if (ipiv[k] == 0) {
                                                if (cabs(a[j][k]) >= big) {
                                                        big=cabs(a[j][k]);
                                                        irow=j;
                                                        icol=k;
                                                }
                 } else if (ipiv[k] > 1) nrerror("GAUSSJ: Singular Matrix-1");
                                }
                ++(ipiv[icol]);
                if (irow != icol) {
                        for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
                        for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
                }
                indxr[i]=irow;
                indxc[i]=icol;
                if ((a[icol][icol].re == 0.0)&&(a[icol][icol].im == 0.0)) 
			nrerror("GAUSSJ: Singular Matrix-2");
                pivinv.re=cinv(a[icol][icol]);
                a[icol][icol].re=1.0;
                a[icol][icol].im=0.0;
                for (l=1;l<=n;l++) a[icol][l] = cmul(a[icol][l],pivinv);
                for (l=1;l<=m;l++) b[icol][l] = cmul(b[icol][l],pivinv);
                for (ll=1;ll<=n;ll++)
                        if (ll != icol) {
                         dum.re=a[ll][icol].re;
                         dum.im=a[ll][icol].im;
                         a[ll][icol].re=0.0;
                         a[ll][icol].im=0.0;
                         for (l=1;l<=n;l++) {
                         c.re=a[ll][icol].re;
                         c.im=a[ll][icol].im;
			 a[ll][l].re = c.re-cmul(c,dum);
			 a[ll][l].im = c.im-cmul(c,dum);
			 }
                        for (l=1;l<=m;l++)
                         c.re=b[ll][icol].re;
                         c.im=b[ll][icol].im;
			 b[ll][l].re = c.re-cmul(c,dum);
			 b[ll][l].im = c.im-cmul(c,dum);
			}
                        }
        }
        for (l=n;l>=1;l--) {
                if (indxr[l] != indxc[l])
                        for (k=1;k<=n;k++)
                                SWAP(a[k][indxr[l]],a[k][indxc[l]]);
        }
        free_ivector(ipiv,1,n);
        free_ivector(indxr,1,n);
        free_ivector(indxc,1,n);
}

#undef SWAP

