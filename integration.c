#include <stdlib.h>
//#include <malloc.h>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <string>

#define QSIMP_EPS     1.0e-7       // accuracy for numerical integration
#define QMIDINF_EPS   1.0e-5       // accuracy for numerical integration
#define QMIDINF2_EPS  1.0e-3       // accuracy for numerical integration
#define QSIMPMID_EPS  1.0e-3       // accuracy for numerical integration

#define QSIMP_MAX     40           // max loops for numerical integration
#define QMIDINF_MAX   40           // max loops for numerical integration 
#define QMIDINF2_MAX  20           // max loops for numerical integration
#define QSIMPMID_MAX  20           // max loops for numerical integration      

#define FUNC_SIMP(x) ((*func)(x))
#define FUNC_MIDINF(a,b) ((*func)(1.0/(a),(b))/((a)*(a)))
#define FUNC_MIDINF2(a,b,c) ((*func)(1.0/(a),(b),(c))/((a)*(a)))
#define FUNC_SIMPMID(v,w,x) ((*func)(v,w,x))
#define FUNC_SIMPMID2(v,w) ((*func)(v,w))

void err_handler(const char*);

double qsimp(double (*func)(double), double a, double b)
{
  int j;
  double s,st,ost,os;
  double trapzd(double (*func)(double), double, double, int);
  
  ost = os = -1.0e30;
  for (j=1;j<=QSIMP_MAX;j++) {
    st=trapzd(func,a,b,j);
    s=(4.0*st-ost)/3.0;
    if (fabs(s-os) < QSIMP_EPS*fabs(os)) return s; 
    if (s == 0.0 && os == 0.0 && j > 6) return s;
    os=s;
    ost=st;
  }
  err_handler("Too many steps in routine qsimp");
  return 0.0;
}

double trapzd(double (*func)(double), double a, double b, int n)
{
  double x,tnm,sum,del;
  static double s;
  int it,j;
  
  if (n == 1) {
    return (s=0.5*(b-a)*(FUNC_SIMP(a)+FUNC_SIMP(b)));
  } else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC_SIMP(x);
    s=0.5*(s+(b-a)*sum/tnm);
    return s;
  }
}

double qsimp2(double (*func)(double), double a, double b)
{
  int j;
  double s,st,ost,os;
  double trapzd2(double (*func)(double), double, double, int);
  
  ost = os = -1.0e30;
  for (j=1;j<=QSIMP_MAX;j++) {
    st=trapzd2(func,a,b,j);
    s=(4.0*st-ost)/3.0;
    if (fabs(s-os) < QSIMP_EPS*fabs(os)) return s; 
    if (s == 0.0 && os == 0.0 && j > 6) return s;
    os=s;
    ost=st;
  }
  err_handler("Too many steps in routine qsimp");
  return 0.0;
}

double trapzd2(double (*func)(double), double a, double b, int n)
{
  double x,tnm,sum,del;
  static double s;
  int it,j;
  
  if (n == 1) {
    return (s=0.5*(b-a)*(FUNC_SIMP(a)+FUNC_SIMP(b)));
  } else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC_SIMP(x);
    s=0.5*(s+(b-a)*sum/tnm);
    return s;
  }
}

double qmidinf(double (*func)(double,double),double a,double b,double x)
{
  int j;
  double s,olds;
  double midinf(double (*func)(double,double),double,double,int,double);
  
  olds = -1.0e30;
  for (j=1;j<=QMIDINF_MAX;j++) {
    s=midinf(func,a,b,j,x);
    if (fabs(s-olds) < QMIDINF_EPS*fabs(olds)) return s;
    if (s == 0.0 && olds == 0.0 && j > 6) return s;
    olds=s;
    // printf("%d %g\n",j,s);
  }
  err_handler("Too many steps in routine qtrap");
  return 0.0;
}

double midinf(double (*func)(double,double),double aa,double bb,int n,double u)
{
  double x,tnm,sum,del,ddel,b,a;
  static double s;
  int it,j;

  b=1.0/aa;
  a=1.0/bb;
  if (n == 1) {
    return (s=(b-a)*FUNC_MIDINF(0.5*(a+b),u));
  } else {
    for(it=1,j=1;j<n-1;j++) it *= 3;
    tnm=it;
    del=(b-a)/(3.0*tnm);
    ddel=del+del;
    x=a+0.5*del;
    sum=0.0;
    for (j=1;j<=it;j++) {
      sum += FUNC_MIDINF(x,u);
      x += ddel;
      sum += FUNC_MIDINF(x,u);
      x += del;
    }
    return (s=(s+(b-a)*sum/tnm)/3.0);
  }
}

double qmidinf2(double (*func)(double,double,double),double a,double b,double x,double y)
{
  int j;
  double s,olds;
  double midinf2(double (*func)(double,double,double),double,double,int,double,double);
  
  olds = -1.0e30;
  for (j=1;j<=QMIDINF2_MAX;j++) {
    s=midinf2(func,a,b,j,x,y);
    if (fabs(s-olds) < QMIDINF2_EPS*fabs(olds)) return s;
    if (s == 0.0 && olds == 0.0 && j > 6) return s;
    printf("step number %d, s=%g, limit %g>%g\n",j,s,
           fabs(s-olds),QMIDINF2_EPS*fabs(olds));
    olds=s;
  }
  printf("Reached max number of steps in routine qtrap\n");

  return s;
}

double midinf2(double (*func)(double,double,double),double aa,double bb,int n,double u,double v)
{
  double x,tnm,sum,del,ddel,b,a;
  static double s;
  int it,j;

  b=1.0/aa;
  a=1.0/bb;
  if (n == 1) {
    return (s=(b-a)*FUNC_MIDINF2(0.5*(a+b),u,v));
  } else {
    for(it=1,j=1;j<n-1;j++) it *= 3;
    tnm=it;
    del=(b-a)/(3.0*tnm);
    ddel=del+del;
    x=a+0.5*del;
    sum=0.0;
    for (j=1;j<=it;j++) {
      sum += FUNC_MIDINF2(x,u,v);
      x += ddel;
      sum += FUNC_MIDINF2(x,u,v);
      x += del;
    }
    return (s=(s+(b-a)*sum/tnm)/3.0);
  }
}

int fflag=0;

double qsimpmid(double (*func)(double,double,double),double a,double b,double c,double d)
{
  int j;
  double s,st,ost,os;
  double midpnt(double (*func)(double,double,double),double,double,int,double,double);

  ost = os = -1.0e30;
  for (j=1;j<=QSIMPMID_MAX;j++) {

    if(j>=16) fflag=1; else fflag=0;

    st=midpnt(func,a,b,j,c,d);
    s=(9.0*st-ost)/8.0;
    // printf("step number %d, st=%g\n",j,s);
    if (fabs(s-os) < QSIMPMID_EPS*fabs(os)) return s;
    if (s == 0.0 && os == 0.0 && j > 6) return s;
    os=s;
    ost=st;
  }
  err_handler("Too many steps in routine qsimpmid");
  return 0.0;
}

double midpnt(double (*func)(double,double,double),double a,double b,int n,double c,double d)
{
  double v,tnm,sum,del,ddel;
  static double s;
  int it,j;
  
  if (n == 1) {
    return (s=(b-a)*FUNC_SIMPMID(0.5*(a+b),c,d));
  } else {
    for(it=1,j=1;j<n-1;j++) it *= 3;
    tnm=it;
    del=(b-a)/(3.0*tnm);
    ddel=del+del;
    v=a+0.5*del;
    sum=0.0;
    for (j=1;j<=it;j++) {
      sum += FUNC_SIMPMID(v,c,d);

      // if(fflag) printf("%g %g %g : %g\n",v,c,d,FUNC_SIMPMID(v,c,d));

      v += ddel;
      sum += FUNC_SIMPMID(v,c,d);

      // if(fflag) printf("%g %g %g : %g\n",v,c,d,FUNC_SIMPMID(v,c,d));

      v += del;
    }
    s=(s+(b-a)*sum/tnm)/3.0;
    return s;
  }
}

double qsimpmid2(double (*func)(double,double),double a,double b,double c)
{
  int j;
  double s,st,ost,os;
  double midpnt2(double (*func)(double,double),double,double,int,double);
 
  ost = os = -1.0e30;
  for (j=1;j<=QSIMPMID_MAX;j++) {
    st=midpnt2(func,a,b,j,c);
    s=(9.0*st-ost)/8.0;
    // printf("step number %d, st=%g, ost=%g, s=%g\n",j,st,ost,s);
    if(fabs(s-os) < QSIMPMID_EPS*fabs(os) && j>6) return s;
    if(s == 0.0 && os == 0.0 && j > 6) return s;
    os=s;
    ost=st;
  }
  err_handler("Too many steps in routine qsimp");
  return 0.0;
}

double midpnt2(double (*func)(double,double),double a,double b,int n,double c)
{
  double v,tnm,sum,del,ddel;
  static double s;
  int it,j;
  
  if (n == 1) {
    return (s=(b-a)*FUNC_SIMPMID2(0.5*(a+b),c));
  } else {
    for(it=1,j=1;j<n-1;j++) it *= 3;
    tnm=it;
    del=(b-a)/(3.0*tnm);
    ddel=del+del;
    v=a+0.5*del;
    sum=0.0;
    for (j=1;j<=it;j++) {
      sum += FUNC_SIMPMID2(v,c);
      v += ddel;
      sum += FUNC_SIMPMID2(v,c);
      v += del;
    }
    s=(s+(b-a)*sum/tnm)/3.0;
    return s;
  }
}
