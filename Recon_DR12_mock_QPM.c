//Basic reconstruction algorithm for mocks.NERSC
//module load fftw/2.1.5.7 
//CC Recon_DR12_mock_QPM.c integration.c function_recon_QPM.c util.c $FFTW_POST_LINK_OPTS -ldrfftw -ldfftw -lm -o recon_xyz_S10
//qsub -I -q debug -l walltime=00:30:00 -l mppwidth=12
//aprun -n 1 ./reconDR12_QPM 7 1 1 2
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
char *recon_mask = cvector(0,200);

int main(int argc, char *argv[]) {

  if(argc<4) err_handler("need arguments: region=9/N 7/S, filestart, fileend, sample \n");
  int weightflag =0;//for diaplcement field weighting, weightflag =1;
  int region;
  //int fileno;
  int filestart;
  int fileend;
  //int w_flag;
  int sample;
  //double RG; //smoothing
  sscanf(argv[1],"%d",&region);  //region = 7 = south, region = 9 = north
  sscanf(argv[2],"%d",&filestart);
  sscanf(argv[3],"%d",&fileend);
  sscanf(argv[4],"%d",&sample);
  //sscanf(argv[3],"%d",&w_flag); 
  //sscanf(argv[4],"%d",&sample);
  //sscanf(argv[4],"%lf",&RG);
  double bias = 1.0;
  double f_growth = 1.0;
  if (sample ==1) bias =2.1; 
  if (sample ==2) bias =2.1;
  if (sample ==3) bias =2.1;
  if (sample ==1) f_growth = 0.7574;
  if (sample ==2) f_growth = 0.7574;
  if (sample ==3) f_growth = 0.7041;
  double RG= 20.0; //smoothing length in Mpc/h
  double L=0.0;
  if (region ==9 && sample ==1) L =3700.0;//extended CMASS
  if (region ==7 && sample ==1) L =3000.0;
  if (region ==9 && sample ==2) L =3500.0;//CMASS
  if (region ==7 && sample ==2) L =3000.0;
  if (sample==3) L=2500.0;
  int NGAL_MAX = 2000000;
  int NRAN_MAX = 58005000;
  double ZMIN =0.43;
  double ZMAX =0.7;
  double area=0.0;
  if (region ==9) area =6851.0;
  if (region ==7) area =2525.0;
  //data files
  //if (region == 9 && sample==1) sprintf(data,"/global/project/projectdirs/boss/galaxy/mmagana/QPM_MOCKS_2014_DR12/north/nonrec///Mocks_and_randoms/xyz%d_red_test.txt",fileno); 
  for (int fileno=filestart; fileno<=fileend; fileno++) {
 /* if (region == 9) sprintf(data,"/global/project/projectdirs/boss/galaxy/QPM/dr12d_cmass_extended/ngc/a0.6452_%04d.dr12d_cmass_extended_ngc.rdz", fileno);
  if (region == 7) sprintf(data,"/global/project/projectdirs/boss/galaxy/QPM/dr12d_cmass_extended/sgc/a0.6452_%04d.dr12d_cmass_extended_sgc.rdz", fileno);*/

/*  if (region == 9 && sample==1) sprintf(data,"/global/project/projectdirs/boss/galaxy/QPM/dr12d_cmass_extended/ngc/a0.6452_%04d.dr12d_cmass_extended_ngc.rdz", fileno);
  if (region == 7 && sample==1) sprintf(data,"/global/project/projectdirs/boss/galaxy/QPM/dr12d_cmass_extended/sgc/a0.6452_%04d.dr12d_cmass_extended_sgc.rdz", fileno);*/

 if  (region == 9 && sample==2) sprintf(data,"/global/project/projectdirs/boss/galaxy/QPM/dr12d_cmass/ngc/a0.6452_%04d.dr12d_cmass_ngc.rdzw", fileno);
 // if (region == 7 && sample==2) sprintf(data,"/global/project/projectdirs/boss/galaxy/QPM/dr12d_cmass/sgc/a0.6452_%04d.dr12d_cmass_sgc.rdzw", fileno);
  //if (region == 9 && sample==2) sprintf(data,"/scratch2/scratchdirs/npadmana/QPM/dr12d_cmass/qpmcosmo/ngc/a0.6452_%04d.dr12d_cmass_ngc.xyzwi", fileno);
/*  if (region == 9 && sample==3) sprintf(data,"/global/project/projectdirs/boss/galaxy/QPM/dr12d_lowz/ngc/a0.7143_%04d.dr12d_lowz_ngc.rdzw",fileno);
  if (region == 7 && sample==3) sprintf(data,"/global/project/projectdirs/boss/galaxy/QPM/dr12d_lowz/sgc/a0.7143_%04d.dr12d_lowz_sgc.rdzw", fileno);

  if (region == 9 && sample==4) sprintf(data,"/global/project/projectdirs/boss/galaxy/QPM/dr12d_combined_beta/bias_mock_test3b_%04d.rdzw",fileno);
 // if (region == 7 && sample==4) sprintf(data,"/global/project/projectdirs/boss/galaxy/QPM/dr12d_combined_beta/ngc/a0.7143_%04d.dr12d_lowz_ngc.rdzw", fileno);
*/
  printf("data = %s\n",data);
  //extended CMASS gc dr12 QPM

  double min_x = 1000.0;
  double min_y = 1000.0;
  double min_z = 1000.0;
 
  if (region ==7 && sample==1) {
    min_x = 400.0-200.0;
    min_y = -1260.0-200.0;
    min_z = -350.0-200.0;
  }
  if (region ==9 && sample==1) {
    min_x = -1850.0-200.0;
    min_y = -1700.0-200.0;
    min_z = -150.0-200.0;
  }

  if (region ==9 && sample==2) {
    min_x = -1770.0  - 200.0;
    min_y = -1600.0 - 200.0;
    min_z = -110.0 - 200.0;

  }
  if (region ==7 && sample==2) {
    min_x = 814.0 - 200.0;
    min_y = -1202.0 - 200.0;
    min_z = -334.5 - 200.0;
  }
  if (region ==7 && sample ==3) {
    min_x = 300.0 - 200.0;
    min_y = -800.0 - 200.0;
    min_z = -200.0 - 200.0;
  }
  if (region ==9 && sample ==3) {
    min_x = -1160.0-200.0;
    min_y = -1100.0-200.0;
    min_z = -100.0 - 200.0;
  }
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
 // int *ipNRAN =&NRAN;
  int NGAL =0;
 // int *ipNGAL =&NGAL;

  for (int i=0; i<NGK; i++) {
      dr[i]=0.0;
      dr_uni[i]=0.0;
      mask[i]=0;
  }
  
  //weights and densities
  double tgal_weight=0.0;
  double tran_weight=0.0;

  FILE *fp_dat;
  if((fp_dat=fopen(data,"r"))==NULL) printf("data not opened\n");
  const int bsz=300; char buf[bsz];  
  while((fgets(buf, bsz, fp_dat))!=NULL) { 
     double ra,dec,cz,fkp, w2,wp;
     if(sample==2 || sample==3 )sscanf(buf,"%lf %lf %lf %lf %lf\n",&ra,&dec,&cz,&fkp,&w2);
     if(sample==1)sscanf(buf,"%lf %lf %lf %lf\n",&ra,&dec,&cz,&wp);
     if(++(NGAL) > NGAL_MAX) { NGAL--; break; }
           if (sample==1) fkp=1.0;
          gal[NGAL].ra = ra*pi/180.;
          gal[NGAL].dec = dec*pi/180.;
          gal[NGAL].z = cz;
          gal[NGAL].dist = calc_dp(gal[NGAL].z);  
          if(NGAL==1) printf("dist=%lf\n",gal[NGAL].dist);          
          gal[NGAL].cp[0] = (gal[NGAL].dist*cos(gal[NGAL].dec)*cos(gal[NGAL].ra)); 
          gal[NGAL].cp[1] = (gal[NGAL].dist*cos(gal[NGAL].dec)*sin(gal[NGAL].ra));
          gal[NGAL].cp[2] = (gal[NGAL].dist*sin(gal[NGAL].dec));
          gal[NGAL].fkp=  fkp;
          gal[NGAL].weight=w2;
          if (sample==2 || sample ==3 || sample ==4 )gal[NGAL].weight=fkp*w2;
          //gal[NGAL].nz = nz;
          //gal[NGAL].weight=fkp;
          tgal_weight += gal[NGAL].weight; 

  }
  
  printf("NGAL = %d, tgal_weights=%lf \n",NGAL,tgal_weight); 

/*  if (region ==9 && sample==1) sprintf(Ran,"/global/project/projectdirs/boss/galaxy/QPM/dr12d_cmass_extended/randoms/a0.6452_rand20x.dr12d_cmass_extended_ngc.rdz"); 
  if (region ==7 && sample==1) sprintf(Ran,"/global/project/projectdirs/boss/galaxy/QPM/dr12d_cmass_extended/randoms/a0.6452_rand20x.dr12d_cmass_extended_sgc.rdz");   */
  if (region ==9 && sample==2) sprintf(Ran,"/scratch1/scratchdirs/angela/a0.6452_rand20x.dr12d_cmass_ngc_vito.txt");
  if (region ==7 && sample==2) sprintf(Ran,"/scratch1/scratchdirs/angela/a0.6452_rand20x.dr12d_cmass_sgc_vito.txt");

/* if (region ==9 && sample==3) sprintf(Ran,"/global/project/projectdirs/boss/galaxy/QPM/dr12d_lowz/randoms/a0.7143_rand20x.dr12d_lowz_ngc.rdz");
  if (region ==7 && sample==3) sprintf(Ran,"/global/project/projectdirs/boss/galaxy/QPM/dr12d_lowz/randoms/a0.7143_rand20x.dr12d_lowz_sgc.rdz");*/
  printf("ran = %s\n",Ran);
  FILE *fp_ran;
  if((fp_ran=fopen(Ran,"r"))==NULL) printf("ran filewith weight flagnot opened\n");
  const int bsz2=300; char buf2[bsz2];  
  while((fgets(buf2, bsz2, fp_ran))!=NULL) { 
     double ra,dec,cz,wp;
     int vtmp;
    // if (sample!=2)sscanf(buf2,"%lf %lf %lf %lf\n",&ra,&dec,&cz,&wp);
     if (sample==2)sscanf(buf2,"%lf %lf %lf %lf %d\n",&ra,&dec,&cz,&wp,&vtmp);
     if(++(NRAN) > NRAN_MAX) { NRAN--; break; }
          uni[NRAN].ra = ra*pi/180.;
          uni[NRAN].dec = dec*pi/180.;
          uni[NRAN].z = cz;
          uni[NRAN].dist = calc_dp(uni[NRAN].z);  
          if(NGAL==1) printf("dist=%lf\n",uni[NRAN].dist);          
          uni[NRAN].cp[0] = (uni[NRAN].dist*cos(uni[NRAN].dec)*cos(uni[NRAN].ra)); 
          uni[NRAN].cp[1] = (uni[NRAN].dist*cos(uni[NRAN].dec)*sin(uni[NRAN].ra));
          uni[NRAN].cp[2] = (uni[NRAN].dist*sin(uni[NRAN].dec));
          uni[NRAN].fkp=  wp;
          uni[NRAN].weight =wp*(1.0-(double)vtmp);
         // if (sample==2)uni[NRAN].fkp=uni[NRAN].weight;
          //gal[NGAL].nz = nz;
          //gal[NGAL].weight=fkp;
          tran_weight += uni[NRAN].weight; 
          
          /*if (uni[NRAN].cp[0]<min_x) min_x=uni[NRAN].cp[0];
          if (uni[NRAN].cp[1]<min_y) min_y=uni[NRAN].cp[1];
          if (uni[NRAN].cp[2]<min_z) min_z=uni[NRAN].cp[2];   //just to get dimensions initially

          if (uni[NRAN].cp[0]>max_x) max_x=uni[NRAN].cp[0];
          if (uni[NRAN].cp[1]>max_y) max_y=uni[NRAN].cp[1];
          if (uni[NRAN].cp[2]>max_z) max_z=uni[NRAN].cp[2];   //just to get dimensio*/

  }
  
  printf("NRAN = %d, tran_weights=%lf\n", NRAN,tran_weight); 
  //create mask.
 /* make_mask_from_rand(NG,NRAN,mask,min_x,min_y,min_z,uni,L);
  for (int i=0; i<NGK; i++) {
    if (mask[i]>0) mask[i]=1;
  } 
  FILE *fp_mask;
  sprintf(recon_mask,"/global/project/projectdirs/boss/galaxy/angela/Recon_files/QPM_DR12_tests/CMASS/mask_DR12_CMASSex_N_from_rands.txt");

  if((fp_mask = fopen(recon_mask,"w"))==0) 
    printf("cannot open pow output file\n"); 
  for (int i=1; i<NGK;i++) {
    fprintf(fp_mask,"%d\n",mask[i]);
  }*/

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
// if (region == 9 && sample ==1)sprintf(mask_file,"/global/project/projectdirs/boss/galaxy/angela/Recon_files/QPM_DR12_tests/CMASS/mask_DR12_CMASSex_N_from_rands.txt");
  //if (region == 7 && sample ==1)sprintf(mask_file,"/global/project/projectdirs/boss/galaxy/angela/Recon_files/QPM_DR12_tests/CMASS/mask_DR12_CMASSex_S_from_rands.txt");


  if (region == 9 && sample ==2)sprintf(mask_file,"/global/project/projectdirs/boss/galaxy/angela/Recon_files/QPM_DR12_tests/CMASS/mask_DR12_CMASS_N_from_rands.txt");
 if (region == 7 && sample ==2)sprintf(mask_file,"/global/project/projectdirs/boss/galaxy/angela/Recon_files/QPM_DR12_tests/CMASS/mask_DR12_CMASS_S_from_rands.txt");
  //if (region == 9 && sample ==3)sprintf(mask_file,"/scratch/scratchdirs/angela/masks/mask_DR12_LOWZ_N_from_rands.txt");
  //if (region == 7 && sample ==3)sprintf(mask_file,"/scratch/scratchdirs/angela/masks/mask_DR12_LOWZ_S_from_rands.txt");
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

  //calculate nbar for the gals and rands.
   // set up estimate of galaxy density & spline
    
    const int NBD=100;
    double *d_vec  = (double*)calloc(NBD,sizeof(double));
    double *n_vec  = (double*)calloc(NBD,sizeof(double));
    double *n2_vec = (double*)calloc(NBD,sizeof(double));
    double D_MIN = calc_dp(ZMIN)- 20.0;
    double D_MAX = calc_dp(ZMAX)+ 20.0;
    double D_BINWIDTH = (D_MAX-D_MIN)/(double)NBD;

    for(long ig=1;ig<NGAL;ig++) {
      double dp =  sqrt(gal[ig].cp[0]*gal[ig].cp[0]+gal[ig].cp[1]*gal[ig].cp[1]+gal[ig].cp[2]*gal[ig].cp[2]);
      int bin = (int)( (double)(dp-D_MIN)/D_BINWIDTH );
      if(bin>=0 && bin<NBD) n_vec[bin]++;
    }

    for(int i=0;i<NBD;i++) {
      double dpmin     = D_MIN + (double)(i)*D_BINWIDTH;
      double dpmax     = D_MIN + (double)(i+1)*D_BINWIDTH;
      double vol_shell = 4./3.*pi*(dpmax*dpmax*dpmax-dpmin*dpmin*dpmin);
      double vol_bin   = vol_shell * (area*pi*pi/(180.*180.)) / (4.*pi);

      d_vec[i] = D_MIN + ((double)(i)+0.5)*D_BINWIDTH;
      n_vec[i] = n_vec[i] / vol_bin;
    }

    spline(d_vec-1,n_vec-1,NBD,1.0e30,1.0e30,n2_vec-1);

    // *********************************************************
    // apply nbar and weights to galaxies & randoms
    // assumes read in weights are systematic only (no FKP)

    for(long ig=1;ig<NGAL;ig++) {	  
      double dp = sqrt(gal[ig].cp[0]*gal[ig].cp[0]+gal[ig].cp[1]*gal[ig].cp[1]+gal[ig].cp[2]*gal[ig].cp[2]);
      splint(d_vec-1,n_vec-1,n2_vec-1,NBD,dp,&gal[ig].nbar);
      //double wfkp = 1.0/(1.0+Pfkp*gal[ig].nbar);
      //gal[ig].wght *= wfkp;
    }

    for(long ir=1;ir<NRAN;ir++) {	  
      double dp = sqrt(uni[ir].cp[0]*uni[ir].cp[0]+uni[ir].cp[1]*uni[ir].cp[1]+uni[ir].cp[2]*uni[ir].cp[2]);
      splint(d_vec-1,n_vec-1,n2_vec-1,NBD,dp,&uni[ir].nbar);
      //double wfkp = 1.0/(1.0+Pfkp*ran[ir].nbar);
      //ran[ir].wght *= wfkp;
    }
  //find the ratios of nbar and weights for gals and rands
  double alpha = tgal_weight/tran_weight;
  printf("alpha=%lf\n",alpha);

  //bin galaxies and randoms. 

  for (int i=1; i<NGAL; i++) {
    int xp = (int)(((gal[i].cp[0]-min_x)/L) * NG);
    int yp = (int)(((gal[i].cp[1]-min_y)/L) * NG);
    int zp = (int)(((gal[i].cp[2]-min_z)/L) * NG);
   
    int ind = zp + 2*(NG/2 +1)*(yp + NG*xp); 

    if (ind >NGK || ind<0) {
      printf("problem with binning gals/rands, xp=%d, yp=%d, zp=%d\n",xp,yp,zp);
      break;
    }
    dr[ind] +=gal[i].weight;     
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
    dr_uni[ind] +=uni[i].weight;    
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
   dr[i] /=(double)NGMAX;
   if (isnan(dr[i])==1){
     printf("NAN at i=%d, dr=%lf\n", i, dr[i]);
    }
    dr_uni[i] /=(double)NGMAX;  
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
            printf("nan detected at i=%d, j=%d, k=%d, mask=%d\n",i,j,k,mask[ind]);
            break;
          }
          if (isnan(dr_uni_sum)==1 ) {
            printf("nan detected at i=%d, j=%d, k=%d, mask=%d\n",i,j,k,mask[ind]);
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

  //float *dr_new = (float*)malloc(sizeof(float)*NGK);
  //for (int i=0; i<NGK; i++) dr_new[i] =0.0;

  //read in second half of randoms and overwrite the first (save memory)
  //Remove RSDs from gals and record the new positions.
  //adjust_psi_interpolate(NG, NGAL, dkx, dky, dkz, min_x, min_y, min_z, f_growth, bias, D_factor, 2, gal, L);

  //Sample the new positions to get a new n(z) for the randoms.
 // distribute_randoms(uni,gal,NGAL,NRAN);//this function replicates the n(r) of the data.

  //Displace the gals and remove RSDs
  adjust_psi_interpolate(NG, NGAL, dkx, dky, dkz, min_x, min_y, min_z, f_growth, bias, D_factor, 3, gal, L,weightflag);//last flag is weightflag

  //Displace the new random field
   adjust_psi_interpolate(NG, NRAN, dkx, dky, dkz, min_x, min_y, min_z, f_growth, bias, D_factor, 1, uni, L,weightflag);

    free(dkx);
    free(dky);
    free(dkz);

  //calculate the new ra, dec, z of the gals and randoms to output.
  output_ra_dec(gal,NGAL);
  output_ra_dec(uni,NRAN);

  //calculate the new redshift
  int no_spline=0;
  no_spline =220;
  double *spline_red = (double*)malloc(sizeof(double)*no_spline);  
  double *spline_dist = (double*)malloc(sizeof(double)*no_spline);  
  double *y2_dist = (double*)malloc(sizeof(double)*no_spline); 

  double red_inc;
  
  red_inc =0.00300000;
  for (int i =0; i<no_spline;i++) {
    spline_red[i]=1.50000e-01 + (double)i*red_inc;
    spline_dist[i]=calc_dp(spline_red[i]);
  }
  spline (spline_dist, spline_red, no_spline, 1.0e30, 1.0e30, y2_dist);

  calc_red_spline(NGAL, gal, spline_dist, spline_red, y2_dist,1);
  calc_red_spline(NRAN,uni,spline_dist, spline_red, y2_dist,2);

  //output all the new values;

  FILE *fp_recon_outD;

/*  if (region==9 && sample ==1)sprintf(recon_data,"/global/project/projectdirs/boss/galaxy/angela/Recon_files/QPM_DR12_tests/eCMASS/north/gals/Rgal_DR12d_N_CMASSex_mock%d_w.txt",fileno);
  if (region==7 && sample ==1)sprintf(recon_data,"/global/project/projectdirs/boss/galaxy/angela/Recon_files/QPM_DR12_tests/eCMASS/south/gals/Rgal_DR12d_S_CMASSex_mock%d_w.txt",fileno);*/

 /* if (region==9 && sample ==2)sprintf(recon_data,"/scratch1/scratchdirs/angela/QPM/dr12d_cmass/ngc/gals/Rgal_DR12d_N_CMASS_mock%d_wvito_weightflag%d.txt",fileno,weightflag);
  if (region==7 && sample ==2)sprintf(recon_data,"/scratch1/scratchdirs/angela/QPM/dr12d_cmass/sgc/gals/Rgal_DR12d_S_CMASS_mock%d_wvito_weightflag%d.txt",fileno,weightflag);*/

  if (region==9 && sample ==2)sprintf(recon_data,"/scratch1/scratchdirs/angela/QPM/dr12d_cmass/ngc/gals/Rgal_DR12d_N_CMASS_mock%d_s10_xyzw.txt",fileno);

 /* if (region==9 && sample ==3)sprintf(recon_data,"/scratch/scratchdirs/angela/Recon_lowz/north/gals/Rgal_DR12d_N_LOWZ_mock%d_w.txt",fileno);
  if (region==7 && sample ==3)sprintf(recon_data,"/scratch/scratchdirs/angela/Recon_lowz/south/gals/Rgal_DR12d_S_LOWZ_mock%d_w.txt",fileno);*/

/*  printf("filename: %s\n",recon_data);
  if((fp_recon_outD = fopen(recon_data,"w"))==0) 
    printf("cannot open pow output file\n"); 
  for (int i=1; i<=NGAL;i++) {
    fprintf(fp_recon_outD,"%e %e %e %e %e\n",gal[i].ra, gal[i].dec, gal[i].z ,gal[i].fkp, gal[i].weight);
  }*/

  printf("filename: %s\n",recon_data);
  if((fp_recon_outD = fopen(recon_data,"w"))==0) 
    printf("cannot open pow output file\n"); 
  for (int i=1; i<=NGAL;i++) {
    fprintf(fp_recon_outD,"%e %e %e %e\n",gal[i].cp[0], gal[i].cp[1], gal[i].cp[2] ,gal[i].weight);
  }

  FILE *fp_recon_outR;
 /* if (region==9 && sample ==1)sprintf(recon_rand,"/global/project/projectdirs/boss/galaxy/angela/Recon_files/QPM_DR12_tests/eCMASS/north/rands/Rrand_DR12v1_N_CMASSex_mock%d_w.txt",fileno);
  if (region==7 && sample ==1)sprintf(recon_rand,"/global/project/projectdirs/boss/galaxy/angela/Recon_files/QPM_DR12_tests/eCMASS/south/rands/Rrand_DR12v1_S_CMASSex_mock%d_w.txt",fileno);*/

/*  if (region==9 && sample ==2)sprintf(recon_rand,"/scratch1/scratchdirs/angela/QPM/dr12d_cmass/ngc/rands/Rrand_DR12vd_N_CMASS_mock%d_wvito_weight%d.txt",fileno,weightflag);
  if (region==7 && sample ==2)sprintf(recon_rand,"/scratch1/scratchdirs/angela/QPM/dr12d_cmass/sgc/rands/Rrand_DR12vd_S_CMASS_mock%d_wvito_weight%d.txt",fileno,weightflag);*/

  if (region==9 && sample ==2)sprintf(recon_rand,"/scratch1/scratchdirs/angela/QPM/dr12d_cmass/ngc/rands/Rrand_DR12vd_N_CMASS_mock%d_s10_xyzw.txt",fileno);

 /* if (region==9 && sample ==3)sprintf(recon_rand,"/scratch/scratchdirs/angela/Recon_lowz/north/rands/Rrand_DR12d_N_LOWZ_mock%d_w.txt",fileno);
  if (region==7 && sample ==3)sprintf(recon_rand,"/scratch/scratchdirs/angela/Recon_lowz/south/rands/Rrand_DR12d_S_LOWZ_mock%d_w.txt",fileno);*/

 /* if((fp_recon_outR = fopen(recon_rand,"w"))==0) 
    printf("cannot open pow output file\n"); 
  for (int i=1; i<=NRAN;i++) {
    fprintf(fp_recon_outR,"%e %e %e %e %e\n",uni[i].ra, uni[i].dec, uni[i].z, uni[i].fkp, uni[i].weight);
  }*/
  if((fp_recon_outR = fopen(recon_rand,"w"))==0) 
    printf("cannot open pow output file\n"); 
  for (int i=1; i<=NRAN;i++) {
    fprintf(fp_recon_outR,"%e %e %e %e\n",uni[i].cp[0], uni[i].cp[1], uni[i].cp[2], uni[i].weight);
  }
  }
}


/**********************/
 

      
 



