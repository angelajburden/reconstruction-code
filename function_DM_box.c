//module to link to Recon_RM_box.c nERSC
#include "header_DM.h"
#include<stdio.h>
#include <srfftw.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include <fstream>
#include <iomanip>
#include<complex>
#include<time.h>
#include<unistd.h>

#define ASEED 12346653

//function to find distance given redshift.
double calc_dp(double red) {
  //printf("red = %lf\n",red);
  double qsimp(double (*func)(double), double, double);
  double dpbit(double); 
  return qsimp(dpbit,0.,red);
}

double dpbit(double z) { 
  return 2997.92458/sqrt((1.+z)*(1.+z)*(1.+z)*Om_0 + (1-Om_0)); 
} 
double Hubble(double z) {
  return H0*sqrt((1.+z)*(1.+z)*(1.+z)*Om_0 + (1-Om_0)); 
}

double Growth(double red) {
  double Ok = 1. - Om_0 - OL_0;
  double H2 = Om_0*(1.+red)*(1.+red)*(1.+red) + Ok*(1.+red)*(1.+red) + OL_0;
  double Om = Om_0*(1.+red)*(1.+red)*(1.+red)/H2;
  printf("Om = %lf\n",Om);
  double OL = OL_0/H2;
  printf("OL = %lf\n",OL);
  double pow_val = exp((4./7)*log(Om));
  double D_z = (2.5*Om)/(((pow_val)-OL + (1. + 0.5*Om)*(1. + (OL/70.)))*(1.+red));
  return D_z;
}

void bin_part_NGP(struct basic_particle *DMpart, int  NG, double dr[], int NPtot, double L, int flag) {

  int NGK = NG*NG*2*(NG/2+1);

  for (int i=1; i<NPtot; i++) {
    int xp = (int)(((DMpart[i].cp[0])/L) * NG);
    int yp = (int)(((DMpart[i].cp[1])/L) * NG);
    int zp = (int)(((DMpart[i].cp[2])/L) * NG);

    if (xp <  0)  xp = NG + xp;
    if (xp >= NG) xp = xp - NG;

    if (yp <  0)  yp = NG + yp;
    if (yp >= NG) yp = yp - NG;

    if (zp <  0)  zp = NG + zp;
    if (zp >= NG) zp = zp - NG;
   
    int ind = zp + 2*(NG/2 +1)*(yp + NG*xp); 

    if (ind >NGK || ind<0) {
      printf("problem with binning gals/rands, xp=%d, yp=%d, zp=%d\n",xp,yp,zp);
      break;
    }
    dr[ind] +=1.0;    
  }  
}

void calc_displacement_field(double RG, double dr[], double dkx[], double dky[], double dkz[], int NG, double L, int flag) {

  for (int i=0;i<NG;i++) 
    for (int l=0;l<NG;l++)
      for (int m=0;m<=(NG/2);m++) {
        double fx, fy, fz;
        if (i<=NG/2) fx =(double)i*2*pi/L;
        if (i> NG/2) fx = (double)(i-NG)*2*pi/L;
       
        if (l<=NG/2) fy =(double)l*2*pi/L;
        if (l> NG/2) fy = (double)(l-NG)*2*pi/L;

        if (m<=NG/2) fz =(double)m*2*pi/L;

        double ks = (fx*fx + fy*fy + fz*fz);  

        int RE = (2*m) +2*(NG/2+1)*(l+NG*i);
        int IM = (2*m+1)+2*(NG/2+1)*(l+NG*i);
        double GaussFT = exp(-0.5*ks*RG*RG);

        dr[RE] *=GaussFT;
        dr[IM] *=GaussFT;

        double dkr = dr[RE];
	double dki = dr[IM];
        // if (i<20 &&  l<20 && m<20) printf("i=%d, l=%d,m=%d,fx=%lf, fy=%lf,fz=%lf\n",i,l,m,fkx,fy,fz);
        //and calculate the displacement field  
        if (ks>0) {
	    dkx[RE] = (dki*fx)/(ks); 
	    dkx[IM] = -(dkr*fx)/(ks); 
            dky[RE] = (dki*fy)/(ks); 
            dky[IM] = -(dkr*fy)/(ks); 
            dkz[RE] = (dki*fz)/(ks); 
            dkz[IM] = -(dkr*fz)/(ks); 

        } 
        //if (i<20 &&  l<20 && m<20) printf("dkx =%lf, dky=%lf, dkz=%lf\n", dkx[RE], dky[RE],dkz[RE]);
        if (ks <=0.0) {
	    dkx[RE] = 0.0;
	    dkx[IM] = 0.0;
            dky[RE] = 0.0;
            dky[IM] = 0.0;
            dkz[RE] = 0.0;
            dkz[IM] = 0.0;
        }        
 
  }
}
void adjust_psi_interpolate(int NG, int NPtot, double dkx[], double dky[], double dkz[], double f_growth, double D_factor, int flag, struct basic_particle *DMpart, double L) {
   double max_shift =0.0;

   for (int p=1; p<NPtot; p++) {

      int i = (int)(((DMpart[p].cp[0])/L) * NG);
      int j = (int)(((DMpart[p].cp[1])/L) * NG);
      int k = (int)(((DMpart[p].cp[2])/L) * NG);
      double dx = ((DMpart[p].cp[0])/L)*(double)NG - (double)i;
      double dy = ((DMpart[p].cp[1])/L)*(double)NG - (double)j;
      double dz = ((DMpart[p].cp[2])/L)*(double)NG - (double)k;

      int i1 = i+1;
      int j1 = j+1;
      int k1 = k+1;

      if (i <  0)  i = NG + i;
      if (i >= NG) i = i - NG;

      if (j <  0)  j = NG + j;
      if (j >= NG) j = j - NG;

      if (k <  0)  k = NG + k;
      if (k >= NG) k = k - NG;

      if (i1 <  0)  i1 = NG + i1;
      if (i1 >= NG) i1 = i1 - NG;

      if (j1 <  0)  j1 = NG + j1;
      if (j1 >= NG) j1 = j1 - NG;

      if (k1 <  0)  k1 = NG + k1;
      if (k1 >= NG) k1 = k1 - NG;

      if(i1<NG && i1>=0 && j1<NG && j1>=0 && k1<NG && k1>=0) {

      int rA = k + 2*(NG/2 +1)*(j +i*NG); 
      int rB = k1 + 2*(NG/2 +1)*(j + i*NG);
      int rC = k + 2*(NG/2 +1)*(j1 + i*NG);
      int rD = k1 + 2*(NG/2 +1)*(j1 + i*NG);
      int rE = k + 2*(NG/2 +1)*(j + i1*NG);
      int rF = k1 + 2*(NG/2 +1)*(j + i1*NG);
      int rG = k + 2*(NG/2 +1)*(j1 + i1*NG);
      int rH = k1 + 2*(NG/2 +1)*(j1 + i1*NG);


      double shiftx = ((1.0-dx)*(1.0-dy)*(1.0-dz)*dkx[rA] +
                      (1.0-dx)*(1.0-dy)*(dz)*dkx[rB] + 
                      (1.0-dx)*(dy)*(1.0-dz)*dkx[rC] + 
                      (1.0-dx)*dy*dz*dkx[rD]+ 
                      dx*(1.0-dy)*(1.0-dz)*dkx[rE] + 
                      dx*(1.0-dy)*dz*dkx[rF] +
                      dx*dy*(1.0-dz)*dkx[rG] +
                      dx*dy*dz*dkx[rH]);

      double shifty = ((1.0-dx)*(1.0-dy)*(1.0-dz)*dky[rA] +
                      (1.0-dx)*(1.0-dy)*(dz)*dky[rB] + 
                      (1.0-dx)*(dy)*(1.0-dz)*dky[rC] + 
                      (1.0-dx)*dy*(dz)*dky[rD]+ 
                      dx*(1.0-dy)*(1.0-dz)*dky[rE] + 
                      dx*(1.0-dy)*dz*dky[rF] +
                      dx*dy*(1.0-dz)*dky[rG] +
                      dx*dy*dz*dky[rH]);

      double shiftz = ((1.0-dx)*(1.0-dy)*(1.0-dz)*dkz[rA] +
                      (1.0-dx)*(1.0-dy)*(dz)*dkz[rB] + 
                      (1.0-dx)*(dy)*(1.0-dz)*dkz[rC] + 
                      (1.0-dx)*dy*(dz)*dkz[rD]+ 
                      dx*(1.0-dy)*(1.0-dz)*dkz[rE] + 
                      dx*(1.0-dy)*dz*dkz[rF] +
                      dx*dy*(1.0-dz)*dkz[rG] +
                      dx*dy*dz*dkz[rH]);
     
   // double dot_psi =  shiftx*cos(gal[p].dec)*sin(gal[p].ra) + shifty*cos(gal[p].dec)*cos(gal[p].ra) + shiftz*sin(gal[p].dec);

    if (fabs(shiftx)>100) printf("galx %d dx=%lf, dy=%lf, dz=%lf,i=%d, j-%d, k=%d \n", p, dkx[rA], dky[rA], dkz[rA],i,j,k );
    if (fabs(shifty)>100) printf("galy %d dx=%lf, dy=%lf, dz=%lf,i=%d, j-%d, k=%d \n", p, dkx[rA], dky[rA], dkz[rA],i,j,k );
    if (fabs(shiftz)>100) printf("galz %d dx=%lf, dy=%lf, dz=%lf,i=%d, j-%d, k=%d \n", p, dkx[rA], dky[rA], dkz[rA],i,j,k );
    /* double x_r = cos(gal[p].dec)*sin(gal[p].ra);
     double y_r = cos(gal[p].dec)*cos(gal[p].ra);
     double z_r = sin(gal[p].dec);

     double shift_realx= shiftx -(f_growth/(bias + f_growth))*(dot_psi*x_r);
     double shift_realy= shifty -(f_growth/(bias + f_growth))*(dot_psi*y_r);
     double shift_realz= shiftz -(f_growth/(bias + f_growth))*(dot_psi*z_r);*/

     //double RSD_x =0;
     //double RSD_y =0;
     double RSD_z =0;


     //RSD_x =f_growth*(shift_realx*x_r + shift_realy*y_r + shift_realz*z_r)*x_r;
     //RSD_y =f_growth*(shift_realx*x_r + shift_realy*y_r + shift_realz*z_r)*y_r;
     RSD_z =f_growth*shiftz;//(shift_realx*x_r + shift_realy*y_r + shift_realz*z_r)*z_r;
    // if ( p > 100 && p <200) printf("RSDx = %lf\n",RSD_x);

     DMpart[p].cp[0] += -D_factor*(shiftx);
     DMpart[p].cp[1] += -D_factor*(shifty);
     DMpart[p].cp[2] += -D_factor*(shiftz + RSD_z);
     
     if (-D_factor*(shiftx)> max_shift) max_shift =-D_factor*(shiftx);
    }
  }
  printf("max shiftx = %lf\n",max_shift); 
}

void psi_rand(int NG, int NPtot, double dkx[], double dky[], double dkz[], double f_growth, double D_factor, int flag, struct basic_rand *RMpart, double L) {
   double max_shift =0.0;

   for (int p=1; p<NPtot; p++) {

      int i = (int)(((RMpart[p].cp[0])/L) * NG);
      int j = (int)(((RMpart[p].cp[1])/L) * NG);
      int k = (int)(((RMpart[p].cp[2])/L) * NG);
      double dx = ((RMpart[p].cp[0])/L)*(double)NG - (double)i;
      double dy = ((RMpart[p].cp[1])/L)*(double)NG - (double)j;
      double dz = ((RMpart[p].cp[2])/L)*(double)NG - (double)k;

      int i1 = i+1;
      int j1 = j+1;
      int k1 = k+1;

      if (i <  0)  i = NG + i;
      if (i >= NG) i = i - NG;

      if (j <  0)  j = NG + j;
      if (j >= NG) j = j - NG;

      if (k <  0)  k = NG + k;
      if (k >= NG) k = k - NG;

      if (i1 <  0)  i1 = NG + i1;
      if (i1 >= NG) i1 = i1 - NG;

      if (j1 <  0)  j1 = NG + j1;
      if (j1 >= NG) j1 = j1 - NG;

      if (k1 <  0)  k1 = NG + k1;
      if (k1 >= NG) k1 = k1 - NG;

      if(i1<NG && i1>=0 && j1<NG && j1>=0 && k1<NG && k1>=0) {

      int rA = k + 2*(NG/2 +1)*(j +i*NG); 
      int rB = k1 + 2*(NG/2 +1)*(j + i*NG);
      int rC = k + 2*(NG/2 +1)*(j1 + i*NG);
      int rD = k1 + 2*(NG/2 +1)*(j1 + i*NG);
      int rE = k + 2*(NG/2 +1)*(j + i1*NG);
      int rF = k1 + 2*(NG/2 +1)*(j + i1*NG);
      int rG = k + 2*(NG/2 +1)*(j1 + i1*NG);
      int rH = k1 + 2*(NG/2 +1)*(j1 + i1*NG);


      double shiftx = ((1.0-dx)*(1.0-dy)*(1.0-dz)*dkx[rA] +
                      (1.0-dx)*(1.0-dy)*(dz)*dkx[rB] + 
                      (1.0-dx)*(dy)*(1.0-dz)*dkx[rC] + 
                      (1.0-dx)*dy*dz*dkx[rD]+ 
                      dx*(1.0-dy)*(1.0-dz)*dkx[rE] + 
                      dx*(1.0-dy)*dz*dkx[rF] +
                      dx*dy*(1.0-dz)*dkx[rG] +
                      dx*dy*dz*dkx[rH]);

      double shifty = ((1.0-dx)*(1.0-dy)*(1.0-dz)*dky[rA] +
                      (1.0-dx)*(1.0-dy)*(dz)*dky[rB] + 
                      (1.0-dx)*(dy)*(1.0-dz)*dky[rC] + 
                      (1.0-dx)*dy*(dz)*dky[rD]+ 
                      dx*(1.0-dy)*(1.0-dz)*dky[rE] + 
                      dx*(1.0-dy)*dz*dky[rF] +
                      dx*dy*(1.0-dz)*dky[rG] +
                      dx*dy*dz*dky[rH]);

      double shiftz = ((1.0-dx)*(1.0-dy)*(1.0-dz)*dkz[rA] +
                      (1.0-dx)*(1.0-dy)*(dz)*dkz[rB] + 
                      (1.0-dx)*(dy)*(1.0-dz)*dkz[rC] + 
                      (1.0-dx)*dy*(dz)*dkz[rD]+ 
                      dx*(1.0-dy)*(1.0-dz)*dkz[rE] + 
                      dx*(1.0-dy)*dz*dkz[rF] +
                      dx*dy*(1.0-dz)*dkz[rG] +
                      dx*dy*dz*dkz[rH]);
     
   // double dot_psi =  shiftx*cos(gal[p].dec)*sin(gal[p].ra) + shifty*cos(gal[p].dec)*cos(gal[p].ra) + shiftz*sin(gal[p].dec);

    if (fabs(shiftx)>100) printf("galx %d dx=%lf, dy=%lf, dz=%lf,i=%d, j-%d, k=%d \n", p, dkx[rA], dky[rA], dkz[rA],i,j,k );
    if (fabs(shifty)>100) printf("galy %d dx=%lf, dy=%lf, dz=%lf,i=%d, j-%d, k=%d \n", p, dkx[rA], dky[rA], dkz[rA],i,j,k );
    if (fabs(shiftz)>100) printf("galz %d dx=%lf, dy=%lf, dz=%lf,i=%d, j-%d, k=%d \n", p, dkx[rA], dky[rA], dkz[rA],i,j,k );
    /* double x_r = cos(gal[p].dec)*sin(gal[p].ra);
     double y_r = cos(gal[p].dec)*cos(gal[p].ra);
     double z_r = sin(gal[p].dec);

     double shift_realx= shiftx -(f_growth/(bias + f_growth))*(dot_psi*x_r);
     double shift_realy= shifty -(f_growth/(bias + f_growth))*(dot_psi*y_r);
     double shift_realz= shiftz -(f_growth/(bias + f_growth))*(dot_psi*z_r);*/

     //double RSD_x =0;
     //double RSD_y =0;
     double RSD_z =0;


     //RSD_x =f_growth*(shift_realx*x_r + shift_realy*y_r + shift_realz*z_r)*x_r;
     //RSD_y =f_growth*(shift_realx*x_r + shift_realy*y_r + shift_realz*z_r)*y_r;
     RSD_z =f_growth*shiftz;//(shift_realx*x_r + shift_realy*y_r + shift_realz*z_r)*z_r;
    // if ( p > 100 && p <200) printf("RSDx = %lf\n",RSD_x);

     RMpart[p].cp[0] += -D_factor*(shiftx);
     RMpart[p].cp[1] += -D_factor*(shifty);
     RMpart[p].cp[2] += -D_factor*(shiftz + RSD_z);
     
     if (-D_factor*(shiftx)> max_shift) max_shift =-D_factor*(shiftx);
    }
  }
  printf("max shiftx = %lf\n",max_shift); 
}







