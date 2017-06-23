#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <string>
//#include <cpgplot.h>

const int   NR_END     = 1;  

const double pi    = 3.1415926536; // Pi
const double G_cst = 6.67e-11;     // gravitational cst G
const double Mpc   = 3.09e22;      // Mpc in meters
const double M_sun = 1.99e30;      // Solar mass in kg
const double MSTAR = -20.44;

//cosmological parameters //DR12 version
const double Om_0 = 0.31;
const double OL_0 = 0.69;
const double a_0 =1;
const double H0 = 67.6;
//const double bias = 1.87;
//const double f_growth = 0.7441;//;//exp(0.6*log(Om_0)) + (OL_0/70.)*( 1 + 0.5*Om_0);0.85;

// define universe we're working in - read from file
extern double sig8;         // value of sigma_8
extern double rho;          // current density of universe (solar mass/Mpc^3)
extern double G0;           // present day density cst
extern double Gl;           // Cosmological Constant densty
//extern double H0;           // Hubble cst km/s.MPc
extern double pow_n;        // Scalar spectral index
extern double bfrac;        // Baryon fraction of total density
extern double pow_norm;     // power spectrum normalisation z=0 only
extern double a_st;         // Sheth & Tormen constant a
extern double p_st;         // Sheth & Tormen constant p
extern double N_st;         // Sheth & Tormen constant N
extern double gom_m;        // time varying Omega_m
extern double gom_v;        // time varying Omega_v



extern unsigned long ija[];  //for spare matrix alogrithm and 
extern double sa[];          //precon. conjugate grad. method

// set up for complex structures
#ifndef _FCOMPLEX_DECLARE_T_
  typedef struct FCOMPLEX {float r,i;} fcomplex;
#define _FCOMPLEX_DECLARE_T_
#endif



// house keeping functions: util.c
void   err_handler(const char*);
double *dvector(long,long);
float  poidev(float, long*);
void   free_dvector(double*,long,long);
float  *vector(long,long);
void   free_vector(float*,long,long);
int    *ivector(long,long);
void   free_ivector(int*,long,long);
char   *cvector(long,long);
void   free_cvector(char*,long,long);
float  **matrix(long,long,long,long);
void   free_matrix(float**,long,long,long,long);
int    **imatrix(long,long,long,long); 
void   free_imatrix(int**,long,long,long,long);
double **dmatrix(long,long,long,long);
void   free_dmatrix(double**,long,long,long,long);
float  ***f3tensor(long,long,long,long,long,long);
void   free_f3tensor(float***,long,long,long,long,long,long);
float  ****f4tensor(long,long,long,long,long,long,long,long);
void   free_f4tensor(float****,long,long,long,long,long,long,long,long);
int    **vimatrix(long,long,long,int*);
void   free_vimatrix(int**,long,long,long);
float  ran2(long*);
float  gasdev(long*);

// chisq optimisation routines: fit_param.c
double golden(double,double,double,double(*f)(double),double,double*);
void powell(double[],double**,int,double,int*,double*,double (*func)(double []));

// integration functions: integration.c
double qsimp(double (*func)(double),double,double);
double qsimp2(double (*func)(double),double,double);
double qmidinf(double (*func)(double,double),double,double,double);
double qmidinf2(double (*func)(double,double,double),double,double,double,double);
double qsimpmid(double (*func)(double,double,double),double,double,double,double);
double qsimpmid2(double (*func)(double,double),double,double,double);

// root finding functions: root.c
double dfridr(double (*func)(double,double),double,double,double*,double);
double zbrent(double (*func)(double,double),double,double,double,double);
void   zroots(float a[], int, fcomplex roots[], int);
void   laguer(fcomplex a[], int, fcomplex*, int*);

// spline fitting routines: spline.c
void   spline(double[],double[],int,double,double,double[]);
void   splint(double[],double[],double[],int,double,double*);

//  matrix inversion: mat_inv.c
void mat_inv(double**,double**,int,double*);
void dsvdcmp(double**,int,int,double[],double**);

// PS mass functions: ps_mass_func.c
void   nm_power();
double Pk(double); 
double Sigth(double,double);
double dSigthbydR(double,double);

//preconditioned conjugate gradient method for multiplying matrices: linbcg.c
void linbcg(unsigned long,double,double,int,double,int,int,double,double,unsigned long);
void atimes(unsigned long,double,double,int,double,unsigned long,double);
void asolve(unsigned long,double,double,int);
void dsprsax(double,unsigned long,double,double,unsigned long);
void dsprstx(double,unsigned long,double,double,unsigned long);
//void sprsin(float**,int,float,unsigned long,float,unsigned long);
double snrm(unsigned long,double,int);

//power spectrum calculation: power.c
void calc_power(float*,int,int,double, char*,double,double,double);

//binning of data
double calc_dp(double);//calculate distance from redshift
double Hubble(double);//calculate Hubble parameter
double Growth(double);//calculate growth function.

struct basic_gal { 
    double ra, dec, z, dist, fkp, nbar, new_dist, weight, ztrue; 
    double cp[3]; 
    double RSD[3]; 
 };
struct basic_particle { 
    double cp[3]; 
    double vcp[3]; 
 };
struct weight_gal { 
    double ra, dec, z, dist, new_dist,nz,fkp; 
    double cp[3]; 
    double wn[7];
    double RSD[3];  
    double weight;
 }; 

void read_in_rands(const char*, int, int*, struct basic_gal*, int, double,double,double,double,double,double);
double read_in_data(const char*, int, int*, struct weight_gal*, double*, double*, double*,int, int);
void make_mask_from_rand(int,int,int mask[],double,double,double, struct basic_gal*,double); 
void bin_gals_NGP(struct basic_gal*, double, double, double, int, float dr[], int, double, int);
void smooth_field(double, double dr[], double dr_uni[], int, double, int);
void calc_displacement_field(double, double del[], double dkx[], double dky[], double dkz[], int, double, int);
void adjust_psi_interpolate(int, int, double dkx[], double dky[], double dkz[], double, double, double, double, double, double, int, struct basic_gal*, double L, int);
void adjust_psi_interpolateW(int, int, double dkx[], double dky[], double dkz[], double, double, double, double, double, double, int, struct weight_gal*, double L);
void distribute_randoms(struct basic_gal*, struct weight_gal*, int, int);
void output_ra_dec(struct basic_gal*, int);
void output_ra_decW(struct weight_gal*, int);
void calc_red_spline(int,struct basic_gal*,double spline_dist [], double spline_red [],double y2[], int);
void calc_red_splineW(int,struct weight_gal*,double spline_dist [], double spline_red [],double y2[], int);
double read_in_xyzw(const char*, int , int*, struct basic_gal* , int );
void read_in_particle(char*, int, int, struct basic_particle*, int);
double read_in_data_weight(const char*, int , int* , struct weight_gal* , int , int);
void bin_data(struct weight_gal*, int, int, int, double *del, double, double, double, double, int);
