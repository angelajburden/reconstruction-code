//Basic reconstruction algorithm for mocks.NERSC to show comparisons of mocks 
//module load fftw/2.1.5.7 
//CC Recon_DR12_mocks_v2.c integration.c function_recon_v2.c util.c $FFTW_POST_LINK_OPTS -ldrfftw -ldfftw -lm -o reconDR12
//qsub -I -q debug -l walltime=00:30:00 -l mppwidth=12
//aprun -n 1 ./reconDR12 9 1
//g++ Recon_DR12_mock.c function_DR12.c util.c integration.c -lrfftw -lfftw -lm -o recon_DR12_NERSC
// ./reconDR12_NERSC 9 1
#include "header_recon_NERSC.h"
#include<stdio.h>
#include <drfftw.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include <fstream>
#include <iomanip>
#include<complex>
#include<omp.h>
#include<time.h>
#define pi2 (pi*pi)
#define eps 1.0e-3// numerical accuracy for integrations
#define NEVAL 10000

const int NG     = 512;

struct basic_gal *gal;  
struct basic_gal *uni; 
 
#define CHUNKSIZE 1; 

//parameters for box
const int NGMAX = NG*NG*NG;	      
const int NGK = NG*NG*2*(NG/2+1);
//const int NRB = 1000;	           //number of redshift bins.
char *data = cvector(0,200);
char *Ran =cvector(0,200);
char *mask_file = cvector(0,200);
char *recon_data = cvector(0,200); 
char *recon_rand = cvector(0,200);

int main(int argc, char *argv[]) {

  if(argc<1) err_handler("need arguments: region=9/N 7/S, weightflag \n");
  int region;
  int fileno;

  sscanf(argv[1],"%d",&region);  //region = 7 = south, region = 9 = north
  sscanf(argv[2],"%d",&fileno); 
  double RG= 15.0; //smoothing length in Mpc/h
  double L =3500.0;
  int NGAL_MAX = 800000;
  int NRAN_MAX = 26005000;

  if (region == 9) sprintf(data,"/global/project/projectdirs/boss/galaxy/mmagana/QPM_MOCKS_2014_DR12/north/nonrec/Mocks_and_randoms/xyz%d_red_test.txt",fileno); 
  printf("data = %s\n",data);

  double min_x = 1000.0;
  double min_y = 1000.0;
  double min_z = 1000.0;


  min_x = -1770.0  - 200.0; //this is just the minimum position of the survey in cartesian padded by 200 for FFTW, I pre-computed this but can be done
  min_y = -1600.0 - 200.0;
  min_z = -110.0 - 200.0;

  double dr_uni_sum =0.;
  double dr_sum =0.;

  //allocate memory for struct
  if(!(gal = (struct basic_gal*)malloc(NGAL_MAX*sizeof(struct basic_gal))-1)) 
    printf("memory allocation problem for galaxies\n");
  if(!(uni = (struct basic_gal*)malloc(NRAN_MAX*sizeof(struct basic_gal))-1)) 
    printf("memory allocation problem for random galaxies\n");

    //allocate memory for grid point vectors
  int *mask = (int*)malloc(sizeof(int)*NGK);
  double *dr = (double*)malloc(sizeof(double)*NGK);
  double *dr_uni = (double*)malloc(sizeof(double)*NGK);

  int NRAN =0;
  int NGAL =0;

  for (int i=0; i<NGK; i++) {
      dr[i]=0.0;
      dr_uni[i]=0.0;
      mask[i]=0;
  }
  

  //weights and densities
  double tgal_weight=0.0;
  double tran_weight=0.0;

  FILE *fp_dat;
  if((fp_dat=fopen(data,"r"))==NULL) printf("data filewith weightflagnot opened\n");
  const int bsz=300; char buf[bsz];  
  while((fgets(buf, bsz, fp_dat))!=NULL) { 
     double xp,yp,zp,wp;
     sscanf(buf,"%lf %lf %lf %lf\n",&xp,&yp,&zp,&wp);
     if(++(NGAL) > NGAL_MAX) { NGAL--; break; }
          if(NGAL==1) printf("x=%lf,y=%lf,z=%lf,w1=%lf\n",xp,yp,zp,wp);
          gal[NGAL].dist = sqrt(xp*xp + yp*yp + zp*zp);   
          if(NGAL==1) printf("dist=%lf\n",gal[NGAL].dist);       
          gal[NGAL].cp[0] = xp; 
          gal[NGAL].cp[1] = yp;
          gal[NGAL].cp[2] = zp;
          gal[NGAL].fkp=  wp;
          tgal_weight+= wp;
	  
  }
  
  printf("NGAL = %d, tgal_weights=%lf \n",NGAL,tgal_weight); 

  sprintf(Ran,"/global/project/projectdirs/boss/galaxy/mmagana/QPM_MOCKS_2014_DR12/north/nonrec/rand20x.txt"); 
  FILE *fp_ran;
  if((fp_ran=fopen(Ran,"r"))==NULL) printf("ran filewith weightflagnot opened\n");
  const int bsz2=300; char buf2[bsz2];  
  while((fgets(buf2, bsz2, fp_ran))!=NULL) { 
     double xp,yp,zp,wp;
     sscanf(buf2,"%lf %lf %lf %lf\n",&xp,&yp,&zp,&wp);
     if(++(NRAN) > NRAN_MAX) { NRAN--; break; }
          if(NRAN==1) printf("x=%lf,y=%lf,z=%lf,w1=%lf\n",xp,yp,zp,wp);
          uni[NRAN].dist = sqrt(xp*xp + yp*yp + zp*zp);   
          if(NGAL==1) printf("dist=%lf\n",uni[NRAN].dist);       
          uni[NRAN].cp[0] = xp; 
          uni[NRAN].cp[1] = yp;
          uni[NRAN].cp[2] = zp;
          uni[NRAN].fkp=  wp;
          tran_weight+= wp;
  }
  
  printf("NRAN = %d, tran_weights=%lf \n",NRAN,tran_weight); 

  //Read in mask.//this shouls really be used but bug in magle mask?
 /*sprintf(mask_file,"/users/angela/lustre/angela/masks/window_vector_cmass_DR12_N.txt");
  FILE *fp_mask;
  if((fp_mask=fopen(mask_file,"r"))==NULL) printf("Mask file not opened\n");
  const int bsz3=100; char buf3[bsz3];
  int inc =0;
  while((fgets(buf3, bsz3, fp_mask))!=NULL) {  
    if(++inc>NGMAX) { inc--; break; }                                     
    int mask_flag; 
    sscanf(buf3,"%d\n",&mask_flag);
    mask[inc] = mask_flag;
    if (mask[inc]>1) printf("mask=%d\n",mask[inc]);
  }*/
  sprintf(mask_file,"/global/project/projectdirs/boss/galaxy/angela/DR12_mask_from_rands_N.txt");
  FILE *fp_mask;
  if((fp_mask=fopen(mask_file,"r"))==NULL) printf("Mask file not opened\n");
  const int bsz3=100; char buf3[bsz3];
  int inc =0;
  while((fgets(buf3, bsz3, fp_mask))!=NULL) {  
    if(++inc>NGK) { inc--; break; }                                     
    int mask_flag; 
    sscanf(buf3,"%d\n",&mask_flag);
    mask[inc] = mask_flag;
    if (mask[inc]>1) printf("mask=%d\n",mask[inc]);
  }


  //find the ratios of nbar and weights for gals and rands
  double alpha = tgal_weight/tran_weight;
  printf("alpha=%lf\n",alpha);

  //bin galaxies and randoms. 

  for (int i=1; i<(NGAL+1); i++) {
    int xp = (int)(((gal[i].cp[0]-min_x)/L) * NG);
    int yp = (int)(((gal[i].cp[1]-min_y)/L) * NG);
    int zp = (int)(((gal[i].cp[2]-min_z)/L) * NG);

    int ind = zp + 2*(NG/2 +1)*(yp + NG*xp); 

   if (isnan(dr[ind])==1){
     printf("NAN at i=%d, dr=%lf\n", ind, dr[ind]);
     printf("xp=%d, yp=%d, zp=%d, ind=%d\n", xp,yp,zp, ind);
    }
    dr[ind] +=gal[i].fkp;  
   if (isnan(dr[ind])==1){
     printf("NAN after fkp added i=%d, dr=%lf\n", ind, dr[ind]);
    }   
  } 
  for (int i=1; i<NRAN; i++) {
    int xp = (int)(((uni[i].cp[0]-min_x)/L) * NG);
    int yp = (int)(((uni[i].cp[1]-min_y)/L) * NG);
    int zp = (int)(((uni[i].cp[2]-min_z)/L) * NG);
   
    int ind = zp + 2*(NG/2 +1)*(yp + NG*xp); 

    if (ind >NGK || ind<0) {
      printf("problem with binning gals/rands, xp=%d, yp=%d, zp=%d\n",xp,yp,zp);
      break;
    }
    dr_uni[ind] +=uni[i].fkp;    
  }  
  printf("gals binned\n");

  //FFT both fields seperately, smooth and RFFT
  rfftwnd_plan dp1;
  dp1 = rfftw3d_create_plan(NG,NG,NG,(fftw_direction)-1,
			       FFTW_MEASURE | FFTW_IN_PLACE);
  rfftwnd_one_real_to_complex(dp1,(fftw_real*)dr,NULL);
  rfftwnd_one_real_to_complex(dp1,(fftw_real*)dr_uni,NULL);

  //free memory
  rfftwnd_destroy_plan(dp1);

  //smooth density field
  smooth_field(RG, dr, dr_uni, NG, L, 0);

  //now reverse Fourier transform
  rfftwnd_plan dp_c2r1;
  dp_c2r1 = rfftw3d_create_plan(NG,NG,NG,(fftw_direction)+1,
			       FFTW_MEASURE | FFTW_IN_PLACE);

  rfftwnd_one_complex_to_real(dp_c2r1,(fftw_complex*)dr,NULL);
  rfftwnd_one_complex_to_real(dp_c2r1,(fftw_complex*)dr_uni,NULL);

  rfftwnd_destroy_plan(dp_c2r1); 
  printf("Reverse FFT done\n");

  //normalise (convension for FFTW)
 for (int i=0; i<NGK; i++) { 
   if (isnan(dr[i])==1){
     printf("NAN from rFFTW at i=%d, dr=%lf\n", i, dr[i]);
     break;
    }
  
   dr[i] /=(double)NGMAX;
   if (isnan(dr[i])==1){
     printf("NAN at i=%d, dr=%lf\n", i, dr[i]);
     break;
    }
    dr_uni[i] /=(double)NGMAX;  
    //if (i >10000 && i <10100) printf( "dr[%d]=%lf\n", i, dr[i]);
  }
  printf("data normalised\n");
  //Apply mask to both fields and calculate the new normalisation.
  for (int i=0; i<NG; i++)
      for (int j=0; j<NG; j++)
        for (int k=0; k<NG; k++) {
          int ind = k + 2*(NG/2 + 1)*(j +NG*i);
          //int ind2 = k + NG*(j +NG*i);
 
         dr[ind] *= (double)mask[ind];
         dr_uni[ind] *= (double)mask[ind];

          dr_sum += dr[ind];
         dr_uni_sum += dr_uni[ind];
          if (isnan(dr_sum)==1 ) {
            printf("nan detected at i=%d, j=%d, k=%d\n",i,j,k);//mask=%d\n,mask[ind]);
            break;
          }
          if (isnan(dr_uni_sum)==1 ) {
          printf("nan detected at i=%d, j=%d, k=%d",i,j,k);//,mask[ind]);
//mask=%d\n"
            break;
           }   
  }
  free(mask);
  printf("dr_sum = %lf\n",dr_sum);
  printf("dr_uni_sum = %lf\n",dr_uni_sum);



  double alpha1 = (dr_sum)/(dr_uni_sum);
  printf("alpha 1 = %lf\n", alpha1);
  for (int i=0; i<NGK; i++) {
    dr_uni[i] *= (alpha1);
  }
  
  double *del= (double*)malloc(sizeof(double)*NGK);

  for (int i=0; i<NGK; i++) del[i] =0.0;
  int D=0;
  for (int i=0; i<NGK; i++) {
    if (dr_uni[i]>0.0) del[i] =(dr[i] - dr_uni[i])/(dr_uni[i]*bias);
    else del[i] =0;
    if (del[i] >5000)  D++;
  }
  free(dr);
  free(dr_uni);
  printf("D =%d\n",D);

  printf("Fourier transform again to calc displacements\n");
  rfftwnd_plan dp;
  dp = rfftw3d_create_plan(NG,NG,NG,(fftw_direction)-1, 
			       FFTW_MEASURE | FFTW_IN_PLACE);
  rfftwnd_one_real_to_complex(dp,(fftw_real*)del,NULL);
  //free memory
  rfftwnd_destroy_plan(dp);

  double *dkx = (double*)malloc(NGK*sizeof(double));
  double *dky = (double*)malloc(NGK*sizeof(double));
  double *dkz = (double*)malloc(NGK*sizeof(double));
  for (int i=0; i<NGK; i++) dkx[i] =dky[i]= dkz[i] =0.0;

  //calculate displacement field
  calc_displacement_field(RG, del, dkx, dky, dkz, NG, L, 1);
  free(del);
  //for (int ii=67371008; ii<67371108; ii++) printf("dkx=%lf, dky=%lf, dkz=%lf\n", dkx[ii], dky[ii], dkz[ii]);
  //now reverse Fourier transform
  rfftwnd_plan dp_c2r;
  dp_c2r = rfftw3d_create_plan(NG,NG,NG,(fftw_direction)+1,
			       FFTW_MEASURE | FFTW_IN_PLACE);

  rfftwnd_one_complex_to_real(dp_c2r,(fftw_complex*)dky,NULL);
  rfftwnd_one_complex_to_real(dp_c2r,(fftw_complex*)dkz,NULL);
  rfftwnd_one_complex_to_real(dp_c2r,(fftw_complex*)dkx,NULL);

  rfftwnd_destroy_plan(dp_c2r); 
 
 for (int i=0; i<NG; i++) 
    for (int j =0; j<NG; j++)
      for (int k=0; k<NG; k++){
        int ind = k + 2*(NG/2 +1)*(j + i*NG);
        dkx[ind] /=(double)NGMAX;
        dky[ind] /=(double)NGMAX;  
        dkz[ind] /=(double)NGMAX;
        if (dkx[ind]>100) {printf("ind=%d, dkx=%lf\n", ind, dkx[ind]);break;}
  }
 
   double D_factor =-1.00;

   //read in second half of randoms and overwrite the first (save memory)
  //Remove RSDs from gals and record the new positions.
  //adjust_psi_interpolate(NG, NGAL, dkx, dky, dkz, min_x, min_y, min_z, f_growth, bias, D_factor, 2, gal, L);

  //Sample the new positions to get a new n(z) for the randoms.
  //distribute_randoms(uni2,gal,NGAL,NRAN2);//this function replicates the n(r) of the data.

  //Displace the gals and remove RSDs
  adjust_psi_interpolate(NG, NGAL, dkx, dky, dkz, min_x, min_y, min_z, f_growth, bias, D_factor, 3, gal, L);

  //Displace the random field
  adjust_psi_interpolate(NG, NRAN, dkx, dky, dkz, min_x, min_y, min_z, f_growth, bias, D_factor, 1, uni, L);

  free(dkx);
  free(dky);
  free(dkz);

  //calculate the new ra, dec, z of the gals and randoms to output.
  //output_ra_dec(gal,NGAL);
  //output_ra_dec(uni2,NRAN2);

  //calculate the new redshift
  /*int no_spline=0;
  no_spline =100;
  double *spline_red = (double*)malloc(sizeof(double)*no_spline);  
  double *spline_dist = (double*)malloc(sizeof(double)*no_spline);  
  double *y2_dist = (double*)malloc(sizeof(double)*no_spline); 

  double red_inc;
  red_inc =0.00300000;
  for (int i =0; i<no_spline;i++) {
    spline_red[i]=4.300000e-01 + (double)i*red_inc;
    spline_dist[i]=calc_dp(spline_red[i]);
  }
  spline (spline_dist, spline_red, no_spline, 1.0e30, 1.0e30, y2_dist);

  calc_red_spline(NGAL, gal, spline_dist, spline_red, y2_dist,1);
  calc_red_spline(NRAN2,uni2, spline_dist, spline_red, y2_dist,2);*/

  //output all the new values;*/

  FILE *fp_recon_outD; //NB the outputs here are in xyz coords, the commented code above should give ra, dec red

  if (region==9)sprintf(recon_data,"/global/project/projectdirs/boss/galaxy/angela/Recon_files/QPM_DR12_tests/CMASS/comparisons/Rgal_DR12v1_N_test_mock%d.txt",fileno);
  if (region==7)sprintf(recon_data,"/global/project/projectdirs/boss/galaxy/angela/Recon_files/QPM_DR12_tests/CMASS/comparisons/Rgal_DR12v1_S_test_mock%d.txt",fileno);
  printf("filename: %s\n",recon_data);
  if((fp_recon_outD = fopen(recon_data,"w"))==0) 
    printf("cannot open pow output file\n"); 
  for (int i=1; i<NGAL;i++) {
    fprintf(fp_recon_outD,"%e %e %e %e\n",gal[i].cp[0], gal[i].cp[1],gal[i].cp[2],gal[i].fkp);
  }

  FILE *fp_recon_outR;
  if (region==9)sprintf(recon_rand,"/global/project/projectdirs/boss/galaxy/angela/Recon_files/QPM_DR12_tests/CMASS/comparisons/Rrand_DR12v1_N_test_mock%d.txt",fileno);
  if (region==7)sprintf(recon_rand,"/global/project/projectdirs/boss/galaxy/angela/Recon_files/QPM_DR12_tests/CMASS/comparisons/Rrand_DR12v1_S_test_mock%d.txt",fileno);
  if((fp_recon_outR = fopen(recon_rand,"w"))==0) 
    printf("cannot open pow output file\n"); 
  for (int i=1; i<NRAN;i++) {
    fprintf(fp_recon_outR,"%e %e %e %e\n",uni[i].cp[0],uni[i].cp[1],uni[i].cp[2],uni[i].fkp);
  }
}


/**********************/
 

      
 



