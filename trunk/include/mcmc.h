


#ifndef mcmc_h
#define mcmc_h


#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>   // www.fftw.org                                                   
#include <FrameL.h>  // from LIGOtools package: www.ldas-sw.ligo.caltech.edu/ligotools 
#include <time.h>
#include <remez.h>   // FIR-filter design routine:  www.janovetz.com/jake              
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/time.h>
#include <stdlib.h>

/*
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/GenerateInspiral.h>
*/

#define TRUE (1==1)
#define FALSE (!TRUE)
#define MvdSdebug !TRUE

#define h2r 0.2617993877991494263  //hr -> rad
#define r2h 3.819718634205488207   //rad -> hr
#define d2r 0.01745329251994329509  //deg -> rad
#define r2d 57.29577951308232311    //rad -> deg
#define c3rd  0.3333333333333333333  // 1/3
#define M0 4.926e-6 //Solar mass in seconds

#define max(A,B) ((A)>(B)?(A):(B))
#define min(A,B) ((A)<(B)?(A):(B))













// Global constants etc., assigned in setconstants()

//  *** PLEASE DON'T ADD ANY NEW ONES, BUT USE THE STRUCTS BELOW INSTEAD (e.g. runpar or mcmcvariables) ***
//      (and if you're bored, go ahead and place the variables below in these structs as well)

//The following global arrays have the size of the number of parameters, i.e. 12 or 15, or ~20 to make sure we're always ok:

int fitpar[12],offsetpar[12];
double truepar[12],pdfsigs[12];
int waveformversion;


//Global variables:

char datadir[99];

int npar,iter,skip,screenoutput,adapt;

int offsetmcmc;
double offsetx;

int corrupd,ncorr,prmatrixinfo;

double temp0;
int nburn,nburn0;

int partemp,savehotchains,prpartempinfo;
double tempmax;

int dosnr,domcmc,domatch,intscrout,writesignal;
int printmuch;

double truespin,truetheta,prior_tc_mean;
int downsamplefactor;

int tempi;
  
double Ms,Mpc,G,c,Mpcs,pi,tpi,mtpi;


int useoldmcmcoutputformat;


double tukeywin;

int inject;
double annealfact;

double *chaintemp;                  // vector of temperatures for individual chains (initialised later) 












// Structure with run parameters.  
// This should eventually include all variables in the input files and replace many of the global variables.
// That also means that this struct must be passed throughout much of the code.
struct runpar{
  //int npar;                     // Number of parameters in the MCMC/template
  int ntemps;		          // Size of temperature ladder
  int mcmcseed;                   // Seed for MCMC
  int selectdata;                 // Select which data set to run on
  int networksize;                // Number of IFOs in the detector network
  int maxIFOdbaseSize;            // Maximum number of IFOs for which the properties are read in (e.g. from mcmc.data)
  int selectifos[9];              // Select IFOs to use for the analysis
  
  //int adapt;                    // Use adaptation or not
  //double *fitpar;
  int *setranpar;                 // Set random true (i.e., injection) parameters 
  int ranparseed;                 // Seed to randomise true parameters (i.e., injection)
  
  double blockfrac;               // Fraction of non-correlated updates that is a block update
  double corrfrac;                // Fraction of MCMC updates that used the correlation matrix
  double mataccfr;                // The fraction of diagonal elements that must improve in order to accept a new covariance matrix
  
  double netsnr;                  // Total SNR of the network
  double targetsnr;               // Target total SNR of the network, scale the distance
  double temps[99];               // Temperature ladder for manual parallel tempering
  double startpar[12];            // Starting parameters for the MCMC chains (may be different from true parameters)
  
  double databeforetc;            // Data stretch in the time domain before t_c to use in the analysis
  double dataaftertc;             // Data stretch in the time domain after t_c to use in the analysis
  double lowfrequencycut;         // Lower frequency cutoff to compute the overlap for
  double highfrequencycut;        // Upper frequency cutoff to compute the overlap for
  
  char infilename[99];            // Run input file name
  char outfilename[99];           // Copy of input file name
  char datainfilename[99];        // Run data input file name
};


//Structure for MCMC variables
struct mcmcvariables{
  int iteri;             // State/iteration number
  int npar;              // Number of parameters in the MCMC/template
  int nparfit;           // Number of parameters in the MCMC that is fitted for
  int ntemps;            // Number of chains in the temperature ladder
  int tempi;             // The current temperature index
  int networksize;       // Number of IFOs in the detector network
  
  
  double temp;           // The current temperature
  double mataccfr;       // The fraction of diagonal elements that must improve in order to accept a new covariance matrix
  double basetime;       // Base of time measurement, get rid of long GPS time format
  
  
  double *histmean;      // Mean of hist block of iterations, used to get the covariance matrix
  double *histdev;       // Standard deviation of hist block of iterations, used to get the covariance matrix
  
  int *corrupdate;       // Switch (per chain) to do correlated (1) or uncorrelated (0) updates
  int *acceptelems;      // Count 'improved' elements of diagonal of new corr matrix, to determine whether to accept it

  double *temps;         // Array of temperatures in the temperature ladder
  double *newtemps;      // New temperature ladder, was used in adaptive parallel tempering
  double *tempampl;      // Temperature amplitudes for sinusoid T in parallel tempering
  double *logL;          // Current log(L)
  double *nlogL;         // New log(L)
  double *dlogL;         // log(L)-log(Lo)
  double *maxdlogL;      // Remember the maximum dlog(L)
  double *sumdlogL;      // Sum of the dlogLs, summed over 1 block of ncorr (?), was used in adaptive parallel tempering, still printed?
  double *avgdlogL;      // Average of the dlogLs, over 1 block of ncorr (?), was used in adaptive parallel tempering
  double *expdlogL;      // Expected dlogL for a flat distribution of chains, was used in adaptive parallel tempering                    
  

  double *corrsig;       // Sigma for correlated update proposals
  int *swapTs1;          // Totals for the columns in the chain-swap matrix
  int *swapTs2;          // Totals for the rows in the chain-swap matrix                                               
  int *acceptprior;      // Check boundary conditions and choose to accept (1) or not(0)
  int *ihist;            // Count the iteration number in the current history block to calculate the covar matrix from

  int **accepted;        // Count accepted proposals
  int **swapTss;         // Count swaps between chains
  double **param;        // The current parameters for all chains
  double **nparam;       // The new parameters for all chains
  double **maxparam;     // The best parameters for all chains (max logL)
  double **sig;          // The standard deviation of the gaussian to draw the jump size from
  double **sigout;       // The sigma that gets written to output
  double **scale;        // The rate of adaptation
  
  double ***hist;        // Store a block of iterations, to calculate the covariance matrix
  double ***covar;       // The Cholesky-decomposed covariance matrix
  
  int seed;              // MCMC seed
  gsl_rng *ran;          // GSL random-number seed
  
  FILE *fout;            // Output-file pointer
  FILE **fouts;          // Output-file pointer array
};


// Structure for spin parameter set with 12 parameters
struct parset{
  double m1;         // mass 1                    
  double m2;         // mass 2                    
  double m;          // total mass                
  double mu;         // mass ratio                
  double mc;         // chirp mass                
  double eta;        // sym mass ratio            
  double tc;         // coalescence time          
  double logdl;      // log-distance              
  double spin;       // magnitude of total spin   
  double kappa;      // L^.S^, cos of angle between L^ & S^
  double longi;      // longitude                 
  double sinlati;    // latitude (sin(delta))     
  double phase;      // wave phase   (phi0)       
  double sinthJ0;    // sin Theta_J0              
  double phiJ0;      // Phi_J0                    
  double alpha;      // Alpha_c                   
  
  double NdJ;        // N^.J_o^; inclination of J_o

  // derived quantities (see also `localpar()'):
  double *loctc;     // vector of `local coalescence times' (w.r.t. FT'd data!)         
  double *localti;   // vector of local altitudes                                       
  double *locazi;    // vector of local azimuths                                        
  double *locpolar;  // vector of local polarisations                                   
};


struct interferometer{
        char name[16];
         int index;
      double lati, longi;             // position coordinates                                     
      double rightarm, leftarm;       // arm directions (counter-clockwise (!) from true north)   
      double radius_pole, radius_eqt; // earth radii at poles & equator (in metres)               
      double lowCut, highCut;         // frequency cutoffs (Hz) (`seismic' & `high-frequency' c.) 
      double before_tc;               // seconds of flat (!) window before `prior_tc_mean'        
      double after_tc;                // seconds of flat window after `prior_tc_mean'             

        char ch1name[128];            // specifications for channel 1                             
        char ch1filepath[128];
        char ch1fileprefix[128];
        char ch1filesuffix[128];
         int ch1filesize; 
         int ch1fileoffset; 
         int ch1doubleprecision;

         int add2channels;            // flag: add 2nd channel to 1st ?                           

        char ch2name[128];            // specifications for channel 2                             
        char ch2filepath[128];
        char ch2fileprefix[128];
        char ch2filesuffix[128];
         int ch2filesize; 
         int ch2fileoffset; 
         int ch2doubleprecision;

        long noiseGPSstart;           // starting time for noise PSD estimation                   
        char noisechannel[128];       // specifications for noise channel                         
        char noisefilepath[128];
        char noisefileprefix[128];
        char noisefilesuffix[128];
         int noisefilesize; 
         int noisefileoffset; 
         int noisedoubleprecision;
      double snr;                    // Save the calculated SNR for each detector
  
      // Elements below this point are determined in `ifoinit()':
      double rightvec[3], leftvec[3], // arm (unit-) vectors            
             normalvec[3],            // normal vector of ifo arm plane  
             orthoarm[3];             // vector orthogonal to right arm and normal vector    
                                      // (usually, but not necessarily identical to 2nd arm) 
      double positionvec[3];          // vector pointing to detector position on earth surface 
      double *raw_noisePSD;
         int PSDsize;
fftw_complex *raw_dataTrafo;
         int FTsize;
         int samplerate;
      double FTstart, deltaFT;
  
      double *noisePSD;               // noise PSD interpolated for the above set of frequencies            
fftw_complex *dataTrafo;              // copy of the dataFT stretch corresponding to above frequencies      
         int lowIndex, highIndex, indexRange;
  
      // Frequency-domain template stuff:
      double *FTin;                   // Fourier transform input                                  
fftw_complex *FTout;                  // FT output (type here identical to `(double) complex')
      double *rawDownsampledWindowedData;     // Copy of raw data, downsampled and windowed
      double *FTwindow;               // Fourier transform input window                           
   fftw_plan FTplan;                  // Fourier transform plan                                   
         int samplesize;              // number of samples (original data)                        
     
};





// Declare functions (prototypes):
void readlocalfile();
void readinputfile(struct runpar *run);
void writeinputfile(struct runpar *run);
void readdatainputfile(struct runpar run, struct interferometer ifo[]);
void setconstants(struct runpar *run);
void set_ifo_data(struct runpar run, struct interferometer ifo[]);

void setrandomtrueparameters(struct runpar *run);
void gettrueparameters(struct parset *par);
void getstartparameters(struct parset *par, struct runpar run);
void getparameterset(struct parset *par, double mc, double eta, double tc, double logdl, double spin, double kappa, 
		     double longi, double sinlati, double phase, double sinthJ0, double phiJ0, double alpha);
void allocparset(struct parset *par, int networksize);
void freeparset(struct parset *par);
void printparset(struct parset par);

void setmcmcseed(struct runpar *run);
void setseed(int *seed);

void mcmc(struct runpar *run, struct interferometer *ifo[]);
void chol(int n, double **A);
void par2arr(struct parset par, double *param);
void arr2par(double *param, struct parset *par);
void par2arrt(struct parset par, double **param); //For the case of parallel tempering
void arr2part(double **param, struct parset *par);
int prior(double *par, int p);

void correlated_mcmc_update(struct interferometer *ifo[], struct parset *state, struct mcmcvariables *mcmc);
void uncorrelated_mcmc_single_update(struct interferometer *ifo[], struct parset *state, struct mcmcvariables *mcmc);
void uncorrelated_mcmc_block_update(struct interferometer *ifo[], struct parset *state, struct mcmcvariables *mcmc);

void write_mcmc_header(struct interferometer *ifo[], struct mcmcvariables mcmc, struct runpar *run);
void write_mcmc_output(struct mcmcvariables mcmc);
void allocate_mcmcvariables(struct mcmcvariables *mcmc);
void free_mcmcvariables(struct mcmcvariables *mcmc);

void update_covariance_matrix(struct mcmcvariables *mcmc);
double anneal_temperature(double temp0, int nburn, int nburn0, int iteri);
void swap_chains(struct mcmcvariables *mcmc);
void write_chain_info(struct mcmcvariables mcmc);




double massratio(double m1, double m2);
void mc2masses(double mc, double eta, double *m1, double *m2);
double chirpmass(double m1, double m2);
double GMST(double GPSsec);
double rightAscension(double longi, double GMST);
double longitude(double rightAscension, double GMST);

double dotproduct(double vec1[3], double vec2[3]);
void facvec(double vec1[3], double fac, double vec2[3]);
void addvec(double vec1[3], double vec2[3], double result[3]);
void normalise(double vec[3]);
void crossproduct(double vec1[3], double vec2[3], double result[3]);
void rotate(double x[3], double angle, double axis[3]);
int righthanded(double x[3], double y[3], double z[3]);
void orthoproject(double x[3], double vec1[3], double vec2[3]);
double angle(double x[3], double y[3]);
void coord2vec(double sinlati, double longi, double x[3]);
void vec2coord(double x[3], double *sinlati, double *longi);

void choleskyinv(double randlibout[], double inverse[]);

void ifoinit(struct interferometer **ifo, int networksize);
void ifodispose(struct interferometer *ifo);
void localpar(struct parset *par, struct interferometer *ifo[], int networksize);
double *filter(int *order, int samplerate, double upperlimit);
double *downsample(double data[], int *datalength, double coef[], int ncoef);
void dataFT(struct interferometer *ifo[], int i, int networksize);
double hann(int j, int N);
double tukey(int j, int N, double r);
void noisePSDestimate(struct interferometer *ifo);
double log_noisePSD(double f, struct interferometer *ifo);
double interpol_log_noisePSD(double f, struct interferometer *ifo);
void writeDataToFiles(struct interferometer *ifo[], int networksize, int mcmcseed);
void writeNoiseToFiles(struct interferometer *ifo[], int networksize, int mcmcseed);
void writeSignalsToFiles(struct interferometer *ifo[], int networksize, int mcmcseed);
void printParameterHeaderToFile(FILE * dump);

void antennaepattern(double altitude, double azimuth, double polarisation,
		     double *Fplus, double *Fcross);
void template(struct parset *par, struct interferometer *ifo[], int ifonr);
void template12(struct parset *par, struct interferometer *ifo[], int ifonr);
		  
//************************************************************************************************************************************************

void template15(struct parset *par, struct interferometer *ifo[], int ifonr);

/*
//void LALHpHc(CoherentGW *waveform, double *hplus, double *hcross, int *l, int length, struct parset *par, struct interferometer *ifo, int ifonr);
void LALHpHc(CoherentGW *waveform, int *l, struct parset *par, struct interferometer *ifo);
//double LALFpFc(CoherentGW *waveform, double *wave, int *l, int length, struct parset *par, int ifonr);
double LALFpFc(CoherentGW *waveform, double *wave, int length, struct parset *par, int ifonr);
void LALfreedom(CoherentGW *waveform);
*/

//************************************************************************************************************************************************

double overlapwithdata(struct parset *par, struct interferometer *ifo[], int ifonr);
double match(struct parset *par, struct interferometer *ifo[], int i, int networksize);
double parmatch(struct parset *par1,struct parset *par2, struct interferometer *ifo[], int networksize);
double paroverlap(struct parset *par1, struct parset *par2, struct interferometer *ifo[], int ifonr);
double vecoverlap(fftw_complex *vec1, fftw_complex *vec2, double * noise, int j1, int j2, double deltaFT);
void signalFFT(fftw_complex * FFTout, struct parset *par, struct interferometer *ifo[], int ifonr);
void computeFishermatrixIFO(struct parset *par, int npar, struct interferometer *ifo[], int networksize, int ifonr, double **matrix);
void computeFishermatrix(struct parset *par, int npar, struct interferometer *ifo[], int networksize, double **matrix);

double ifo_loglikelihood(struct parset *par, struct interferometer *ifo[], int i);
double signaltonoiseratio(struct parset *par, struct interferometer *ifo[], int i);
double net_loglikelihood(struct parset *par, int networksize, struct interferometer *ifo[]);
double logprior(struct parset *par);
double logstartdens(struct parset *par);
double logmultiplier(double mc, double eta, double logdl, double sinlati, double cosiota);
double lgamma(double x);
double logit(double x);
double logitinverse(double x);



#endif
