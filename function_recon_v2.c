//module to link to Recon_DR12_weight.c nERSC
#include "header_recon_NERSC.h"
#include<stdio.h>
#include <drfftw.h>
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

//void read_in_rands(const char *data, int NGAL_MAX, int *NGAL, struct weight_gal *gal, int flag, double minx, double maxx, double miny, double maxy, double minz, double maxz) {
double read_in_data(const char *data, int NGAL_MAX, int *NGAL, struct weight_gal *gal, int flag, int flag2){
  double tgal_weights=0.0;
  double minx=1000.0;
  double miny=1000.0;
  double minz=1000.0;
  FILE *fp_rand;
  if((fp_rand=fopen(data,"r"))==NULL) printf("data filewith weightflag %d not opened\n", flag);
  const int bsz=300; char buf[bsz];  
  if (flag2 ==2) {
    while((fgets(buf, bsz, fp_rand))!=NULL) { 
     if (buf[0] =='#'|| buf[0] ==0) continue; 
     double ra, dec, cz;
     double w1,w2,w3,w4,w5,w6,w7;
     sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&ra,&dec,&cz,&w1,&w2,&w3,&w4,&w5,&w6,&w7);
     if(++(*NGAL) > NGAL_MAX) { *NGAL--; break; }
          if(*NGAL==1) printf("ra=%lf,dec=%lf,z=%lf,w1=%lf,w2=%lf,w3=%lf,w4=%lf, w5=%lf, w6=%lf, w7=%lf\n",ra,dec,cz,w1,w2,w3,w4,w5,w6,w7);
	  gal[*NGAL].ra = ra*pi/180.;
          gal[*NGAL].dec = dec*pi/180.;
          gal[*NGAL].z = cz;
          gal[*NGAL].dist = calc_dp(gal[*NGAL].z);   
          if(*NGAL==1) printf("dist=%lf\n",gal[*NGAL].dist);       
          gal[*NGAL].cp[0] = (gal[*NGAL].dist*cos(gal[*NGAL].dec)*cos(gal[*NGAL].ra)); 
          gal[*NGAL].cp[1] = (gal[*NGAL].dist*cos(gal[*NGAL].dec)*sin(gal[*NGAL].ra));
          gal[*NGAL].cp[2] = (gal[*NGAL].dist*sin(gal[*NGAL].dec));
          gal[*NGAL].wn[0]=  w1;
          gal[*NGAL].wn[1]=  w2;
          gal[*NGAL].wn[2]=  w3;
          gal[*NGAL].wn[3]=  w4;
          gal[*NGAL].wn[4]=  w5;
          gal[*NGAL].wn[5]=  w6;
          gal[*NGAL].wn[6]=  w7;
          if (flag>0){
            int wI = flag-1;
            gal[*NGAL].weight=gal[*NGAL].wn[wI];
          }
          if (flag==0) gal[*NGAL].weight=1.0;
          tgal_weights += gal[*NGAL].weight;
        /*  if (gal[*NGAL].cp[0]>maxx) maxx=gal[*NGAL].cp[0];
          if (gal[*NGAL].cp[1]>maxy) maxy=gal[*NGAL].cp[1];
          if (gal[*NGAL].cp[2]>maxz) maxz=gal[*NGAL].cp[2];*/

          if (gal[*NGAL].cp[0]<minx) minx=gal[*NGAL].cp[0];
          if (gal[*NGAL].cp[1]<miny) miny=gal[*NGAL].cp[1];
          if (gal[*NGAL].cp[2]<minz) minz=gal[*NGAL].cp[2];   //just to get dimensions initially*/
     }
   }
   if (flag2 ==0) {
    while((fgets(buf, bsz, fp_rand))!=NULL) { 
     if (buf[0] =='#'|| buf[0] ==0) continue; 
     double ra, dec, cz;
     double fkp, nz, cp, noz, star, see, sytot;
     sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %*lf\n",&ra, &dec, &cz, &fkp, &cp, &noz, &star, &see, &sytot,&nz);
     if (cz>0.43 && cz<0.70){
     if(++(*NGAL) > NGAL_MAX) { *NGAL--; break; }
          if(*NGAL==1) printf("ra=%lf,dec=%lf,z=%lf,fkp=%lf,cp=%lf,noz=%lf,star=%lf, see=%lf, nz=%lf\n",ra,dec,cz,fkp, cp, noz, star, see, sytot, nz);
	  gal[*NGAL].ra = ra*pi/180.;
          gal[*NGAL].dec = dec*pi/180.;
          gal[*NGAL].z = cz;
          gal[*NGAL].dist = calc_dp(gal[*NGAL].z);            
          gal[*NGAL].cp[0] = (gal[*NGAL].dist*cos(gal[*NGAL].dec)*cos(gal[*NGAL].ra)); 
          gal[*NGAL].cp[1] = (gal[*NGAL].dist*cos(gal[*NGAL].dec)*sin(gal[*NGAL].ra));
          gal[*NGAL].cp[2] = (gal[*NGAL].dist*sin(gal[*NGAL].dec));
          gal[*NGAL].fkp=  fkp;
          gal[*NGAL].nz = nz;
          gal[*NGAL].weight=fkp*sytot;
          tgal_weights += gal[*NGAL].weight;
          
          if (gal[*NGAL].cp[0]<minx) minx=gal[*NGAL].cp[0];
          if (gal[*NGAL].cp[1]<miny) miny=gal[*NGAL].cp[1];
          if (gal[*NGAL].cp[2]<minz) minz=gal[*NGAL].cp[2];   //just to get dimensions initially
    }
  }
  }

  if (flag2 ==3) {
    while((fgets(buf, bsz, fp_rand))!=NULL) { 
     if (buf[0] =='#'|| buf[0] ==0) continue; 
     double ra, dec, cz;
     double fkp, nz;
     sscanf(buf,"%lf %lf %lf %lf %lf\n",&ra, &dec, &cz, &fkp,&nz);
     if (cz>0.43 && cz<0.70){
     if(++(*NGAL) > NGAL_MAX) { *NGAL--; break; }
          if(*NGAL==1) printf("ra=%lf,dec=%lf,z=%lf,fkp=%lf,nz=%lf\n",ra,dec,cz,fkp,nz);
	  gal[*NGAL].ra = ra*pi/180.;
          gal[*NGAL].dec = dec*pi/180.;
          gal[*NGAL].z = cz;
          gal[*NGAL].dist = calc_dp(gal[*NGAL].z);  
          if(*NGAL==1) printf("dist=%lf\n",gal[*NGAL].dist);          
          gal[*NGAL].cp[0] = (gal[*NGAL].dist*cos(gal[*NGAL].dec)*cos(gal[*NGAL].ra)); 
          gal[*NGAL].cp[1] = (gal[*NGAL].dist*cos(gal[*NGAL].dec)*sin(gal[*NGAL].ra));
          gal[*NGAL].cp[2] = (gal[*NGAL].dist*sin(gal[*NGAL].dec));
          gal[*NGAL].fkp=  fkp;
          gal[*NGAL].nz = nz;
          gal[*NGAL].weight=fkp;
          tgal_weights += gal[*NGAL].weight;

          if (gal[*NGAL].cp[0]<minx) minx=gal[*NGAL].cp[0];
          if (gal[*NGAL].cp[1]<miny) miny=gal[*NGAL].cp[1];
          if (gal[*NGAL].cp[2]<minz) minz=gal[*NGAL].cp[2];   //just to get dimensions initially
    }
  }
  }
  fclose(fp_rand);
  printf("minx =%lf, miny=%lf, minz=%lf\n",minx,miny,minz);
  return tgal_weights;
}


void make_mask_from_rand(int NG, int NGAL, int mask[],double min_x, double min_y, double min_z, struct basic_gal *gal, double L) {

   for (int p=1; p<NGAL; p++) {
      int i = (int)(((gal[p].cp[0]-min_x)/L) * NG);
      int j = (int)(((gal[p].cp[1]-min_y)/L) * NG);
      int k = (int)(((gal[p].cp[2]-min_z)/L) * NG);

      int i2 = i+1;
      int j2 = j+1;
      int k2 = k+1;

      int i1 = i-1;
      int j1 = j-1;
      int k1 = k-1;
      int kint, iint, jint;
      for (int iv=0;iv<3;iv++)
        for (int jv=0;jv<3;jv++)
          for (int kv=0;kv<3;kv++){
            
            if (kv==0) kint = k1;
            if (kv==1) kint = k;
            if (kv==2) kint =k2;

            if (jv==0) jint = j1;
            if (jv==1) jint = j;
            if (jv==2) jint =j2;

            if (iv==0) iint = i1;
            if (iv==1) iint = i;
            if (iv==2) iint =i2;

            int Nijk = kint + 2*(NG/2 +1)*(jint +iint*NG); 
            mask[Nijk] ++;
      }
  }  
}      
void bin_gals_NGP(struct basic_gal *gal, double  min_x, double min_y, double min_z, int  NG, float dr[], int NGAL, double L, int flag) {

  int NGK = NG*NG*2*(NG/2+1);

  for (int i=1; i<NGAL; i++) {
    int xp = (int)(((gal[i].cp[0]-min_x)/L) * NG);
    int yp = (int)(((gal[i].cp[1]-min_y)/L) * NG);
    int zp = (int)(((gal[i].cp[2]-min_z)/L) * NG);
   
    int ind = zp + 2*(NG/2 +1)*(yp + NG*xp); 

    if (ind >NGK || ind<0) {
      printf("problem with binning gals/rands, xp=%d, yp=%d, zp=%d\n",xp,yp,zp);
      break;
      //printf("xp= %d, yp =%d, zp=%d\n",xp,yp,zp);
      //printf("i =%d, galx = %lf, dist =%lf, z= %lf\n",i,gal[i].cp[0],gal[i].dist,gal[i].z);
    }
    dr[ind] +=gal[i].fkp;
    
  }  
}
void smooth_field(double RG, double dr[],double dr_uni[], int NG, double L, int flag) {

  double GaussFT=0.0;
  for (int i=0;i<NG;i++) 
    for (int l=0;l<NG;l++)
      for (int m=0;m<=(NG/2);m++) {
        double fx, fy, fz;
        if (i<=NG/2) fx =(double)i*2*pi/L;
        if (i> NG/2) fx = (double)(i-NG)*2*pi/L;
       
        if (l<=NG/2) fy =(double)l*2*pi/L;
        if (l> NG/2) fy = (double)(l-NG)*2*pi/L;

        if (m<=NG/2 ) fz =(double)m*2*pi/L;

        double ks = (fx*fx + fy*fy + fz*fz);  
        GaussFT = exp(-0.5*ks*RG*RG);

        int RE = (2*m) +2*(NG/2+1)*(l+NG*i);
        int IM = (2*m+1)+2*(NG/2+1)*(l+NG*i);

        if (flag ==0) {      
          dr[RE] *=GaussFT;
          dr[IM] *=GaussFT;
	  dr_uni[RE] *=GaussFT;
          dr_uni[IM] *=GaussFT;
        }

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
void adjust_psi_interpolate(int NG, int NGAL, double dkx[], double dky[], double dkz[], double min_x, double min_y, double min_z, double f_growth, double bias, double D_factor, int flag, struct basic_gal *gal, double L) {
   double max_shift =0.0;
   //double max_RSD =-100.0;
   //int NGK = NG*NG*2*(NG/2+1);
   for (int p=1; p<NGAL; p++) {

      int i = (int)(((gal[p].cp[0]-min_x)/L) * NG);
      int j = (int)(((gal[p].cp[1]-min_y)/L) * NG);
      int k = (int)(((gal[p].cp[2]-min_z)/L) * NG);
      double dx = ((gal[p].cp[0]-min_x)/L)*(double)NG - (double)i;
      double dy = ((gal[p].cp[1]-min_y)/L)*(double)NG - (double)j;
      double dz = ((gal[p].cp[2]-min_z)/L)*(double)NG - (double)k;

      int i1 = i+1;
      int j1 = j+1;
      int k1 = k+1;

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
     
     double dot_psi =  shiftx*cos(gal[p].dec)*sin(gal[p].ra) + shifty*cos(gal[p].dec)*cos(gal[p].ra) + shiftz*sin(gal[p].dec);
     //double fb = f_growth/bias;
    if (fabs(shiftx)>100) printf("galx %d dx=%lf, dy=%lf, dz=%lf,i=%d, j-%d, k=%d \n", p, dkx[rA], dky[rA], dkz[rA],i,j,k );
    if (fabs(shifty)>100) printf("galy %d dx=%lf, dy=%lf, dz=%lf,i=%d, j-%d, k=%d \n", p, dkx[rA], dky[rA], dkz[rA],i,j,k );
    if (fabs(shiftz)>100) printf("galz %d dx=%lf, dy=%lf, dz=%lf,i=%d, j-%d, k=%d \n", p, dkx[rA], dky[rA], dkz[rA],i,j,k );
     double x_r = cos(gal[p].dec)*cos(gal[p].ra);
     double y_r = cos(gal[p].dec)*sin(gal[p].ra);
     double z_r = sin(gal[p].dec);

     double shift_realx= shiftx -(f_growth/(bias + f_growth))*(dot_psi*x_r);
     double shift_realy= shifty -(f_growth/(bias + f_growth))*(dot_psi*y_r);
     double shift_realz= shiftz -(f_growth/(bias + f_growth))*(dot_psi*z_r);

     double RSD_x =0;
     double RSD_y =0;
     double RSD_z =0;


     RSD_x =f_growth*(shift_realx*x_r + shift_realy*y_r + shift_realz*z_r)*x_r;
     RSD_y =f_growth*(shift_realy*y_r + shift_realy*y_r + shift_realz*z_r)*y_r;
     RSD_z =f_growth*(shift_realz*z_r + shift_realy*y_r + shift_realz*z_r)*z_r;
    // if ( p > 100 && p <200) printf("RSDx = %lf\n",RSD_x);
     if (flag ==1) {
         gal[p].cp[0] += -D_factor*(shift_realx);
         gal[p].cp[1] += -D_factor*(shift_realy);
         gal[p].cp[2] += -D_factor*(shift_realz);
         gal[p].new_dist = sqrt(gal[p].cp[0]*gal[p].cp[0] + gal[p].cp[1]*gal[p].cp[1] +gal[p].cp[2]*gal[p].cp[2]);
         if (-D_factor*(shift_realx)> max_shift) max_shift =-D_factor*(shift_realx);
     }

     if (flag ==3) { 
       gal[p].cp[0] += -D_factor*(shift_realx+ RSD_x); 
       gal[p].cp[1] += -D_factor*(shift_realy+ RSD_y); 
       gal[p].cp[2] += -D_factor*(shift_realz+ RSD_z); 
       if (-D_factor*(shift_realx)> max_shift) max_shift =-D_factor*(shift_realx);
     }   

/*     if (flag ==4) { 
       gal[p].cp[0] += -D_factor*(shift_realx+ RSD_x); 
       gal[p].cp[1] += -D_factor*(shift_realy+ RSD_y); 
       gal[p].cp[2] += -D_factor*(shift_realz+ RSD_z); 
       if (-D_factor*(shift_realx)> max_shift) max_shift =-D_factor*(shift_realx);
       if (shift_realx>50.0) printf("p=%d, shift_realx=%lf, x =%lf, y=%lf, z=%lf, i=%d, j=%d, k=%d\n", p, shift_realx, gal[p].cp[0], gal[p].cp[1],gal[p].cp[2],i,j,k);
       gal[p].wn[0] = -D_factor*(shift_realx+ RSD_x); 
       gal[p].wn[1] = -D_factor*(shift_realy+ RSD_y); 
       gal[p].wn[2] = -D_factor*(shift_realz+ RSD_z); 

     }     
    }
  }*/
 
  }
  }
  printf("max shiftx = %lf\n",max_shift); 
}

void distribute_randoms(struct basic_gal *uni, struct weight_gal *gal, int NGAL, int NRAN) {
 printf("NRAN = %d, NGAL = %d\n",NRAN,NGAL);
 //int replace_gal =6000;
 srand(time(NULL) + getpid());
 for (int i=1; i<NRAN; i++) {

    double dgal = (double)NGAL * rand()/(RAND_MAX + 1.0);
    int igal = (int) dgal;
    if (igal==0) igal = NGAL-1;
    double new_dist = gal[igal].new_dist;
    uni[i].cp[0] = new_dist*cos(uni[i].dec)*cos(uni[i].ra); 
    uni[i].cp[1] = new_dist*cos(uni[i].dec)*sin(uni[i].ra);
    uni[i].cp[2] = new_dist*sin(uni[i].dec);
    if(i ==35546)printf("uni[%d]x = %lf\n",i,uni[i].cp[0]);
    uni[i].weight = gal[igal].weight;
    uni[i].new_dist = new_dist;
  }

}
void output_ra_dec(struct basic_gal *gal, int NGAL) {
  for (int i=1; i<NGAL; i++) {
    gal[i].new_dist = sqrt(gal[i].cp[0]*gal[i].cp[0] + gal[i].cp[1]*gal[i].cp[1] +gal[i].cp[2]*gal[i].cp[2]);
    gal[i].dec = asin(gal[i].cp[2]/gal[i].new_dist);
    gal[i].ra = atan2(gal[i].cp[1],gal[i].cp[0]);
    if (gal[i].ra<0.000000)gal[i].ra +=(2*pi);            
  }
}
void calc_red_spline(int NGAL, struct basic_gal *gal, double spline_dist [], double spline_red [],double y2_dist[], int flag) {

  for (int i=1; i<NGAL; i++) {
    double galred;
    double galdist;
    galdist = gal[i].new_dist;
    splint (spline_dist, spline_red, y2_dist, 7, galdist, &galred);
    gal[i].z = galred;
  }
}
double read_in_xyzw(const char *data, int NGAL_MAX, int *NGAL, struct basic_gal *gal, int flag){
  double tgal_weight =0.0;
  FILE *fp_rand;
  if((fp_rand=fopen(data,"r"))==NULL) printf("data filewith weightflag %d not opened\n", flag);
  const int bsz=300; char buf[bsz];  
  while((fgets(buf, bsz, fp_rand))!=NULL) { 
    // if (buf[0] =='#'|| buf[0] ==0) continue; 
     double xp,yp,zp,wp;
     sscanf(buf,"%lf %lf %lf %lf\n",&xp,&yp,&zp,&wp);
     if(++(*NGAL) > NGAL_MAX) { *NGAL--; break; }
          if(*NGAL==1) printf("x=%lf,y=%lf,z=%lf,w1=%lf\n",xp,yp,zp,wp);
          gal[*NGAL].dist = sqrt(xp*xp + yp*yp + zp*zp);   
          if(*NGAL==1) printf("dist=%lf\n",gal[*NGAL].dist);       
          gal[*NGAL].cp[0] = xp; 
          gal[*NGAL].cp[1] = yp;
          gal[*NGAL].cp[2] = zp;
          gal[*NGAL].fkp=  wp;
          tgal_weight+= wp;
     }
   return tgal_weight;
}






