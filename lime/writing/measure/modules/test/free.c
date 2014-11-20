#include "measure.h"
static const full_spinor_dble fsd0={{{0.0}}};
static const fsd_to_mat(full_spinor_dble *A, complex_dble **a){
 int d1,d2,c1,c2; /* Dirac and Color indices */ 
 for (d1=0;d1<4;d1++){
  for (d2=0;d2<4;d2++){
   for (c1=0;c1<3;c1++){
    for (c2=0;c2<3;c2++){
     i = d1*3 + c1;
     j = d2*3 + c2; 
     (*(a[i][j])).re = (*(A+n)).re;
     (*(a[i][j])).im = (*(A+n)).im;
     n++;
    }
   }
  }
 }
}
static void unit_fsd(full_spinor_dble fsd){  
 /* generate unit full_spinor_dble matrix */
 fsd=fsd0;
 fsd.c11.c11.re=1.0;
 fsd.c11.c22.re=1.0;
 fsd.c11.c33.re=1.0;
 fsd.c22.c11.re=1.0;
 fsd.c22.c22.re=1.0;
 fsd.c22.c33.re=1.0;
 fsd.c33.c11.re=1.0;
 fsd.c33.c22.re=1.0;
 fsd.c33.c33.re=1.0;
 fsd.c44.c11.re=1.0;
 fsd.c44.c22.re=1.0;
 fsd.c44.c33.re=1.0;
 return fsd;
}
static void add_fsd(full_spinor_dble fsd1, fulfull_spinor_dble fsd2){
 int i;
 double *dum1,*dum2;
 full_spinor_dble *fsd3;
 dum1=(double*)&fsd1;
 dum2=(double*)&fsd2;
 for (i=0;i<144;i++){
  *(fsd3+i)=*(dum1+i) + *(dum2+i);
 }
 return (*fsd3);
}
static full_spinor_dble pspacedenom(double kappa,int *L, int *n){
 int i,mu;
 double dp[4]; 
 complex_dble kappa_c;
 full_spinor_dble unit,dum,dum1,dum2,dum3;
 dum  = fsd0;

 dp[0] = 2*pi/L[0]*n[0];
 dp[1] = 2*pi/L[1]*n[1];
 dp[2] = 2*pi/L[2]*n[2];
 dp[3] = 2*pi/L[3]*n[3];

 unit = unit_fsd(unit);
 kappa_c.re=2.*kappa;
 kappa_c.im=0.0;
 f.re = -1.*cos(dp[mu]);
 f.im = 0.0;
 mul_cmplx_fsv_dble(f,unit,dum1);
 for(mu=0;mu<4;mu++){
  f.re = 0.0;
  f.im = +1.*sin(dp[mu]);
  mul_cmplx_fsv_dble(f,unit,dum2);
  mul_gamma_l(mu,dum2,dum3);
  dum2 = add_fsd(dum1,dum3);
  dum  = add_fsd(dum,dum2); 
 }
 mul_cmplx_fsv_dble(kappa_c,unit,dum2)
 return add_fsd(unit,dum); 
}


void main(){
 double kappa;
 complex_dble **a,**b,**x; 
 int i,L[4],n[4];
 complex_dble kappa_c,f;
 full_spinor_dble A,B,X;

 X  = fsd0;
 a = amalloc(12*sizeof(*complex_dble));
 b = amalloc(12*sizeof(*complex_dble));
 c = amalloc(12*sizeof(*complex_dble));
 for (i=0;i<12;i++){
  a[i] = amalloc(12*sizeof(complex_dble));
  b[i] = amalloc(12*sizeof(complex_dble));
  x[i] = amalloc(12*sizeof(complex_dble));
 }

 n[0]=0;
 n[1]=0;
 n[2]=0;
 n[3]=0;
 B = unit_fsd(B);
 A = pspacedenom(kappa,L,n);
 fsd_to_mat(A,a);
 fsd_to_mat(B,b);
 fsd_to_mat(X,x);

 gaussj(a,12,b,12);  


}
