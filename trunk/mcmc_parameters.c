// mcmc_parameters.c:
// Routines to read/write input files, set constants and set true and null parameters

#include <mcmc.h>





// Read the input file.
// Please keep this routine in sync with writeinputfile() below.
// All parameters that are read in here should become members of the runvar struct and lose their global status
void readinputfile(struct runpar *run)
{
  int i;
  double tmpdbl;
  char bla[500];
  FILE *fin;
  
  if((fin = fopen(run->infilename,"r")) == NULL) {
    printf("   Error reading file: %s, aborting.\n\n\n",run->infilename);
    exit(1);
  } else {
    printf("   Using input file: %s.\n",run->infilename);
  }
  
  
  //Use and l for floats: %lf, %lg, etc, rather than %f, %g
  
  for(i=1;i<=3;i++) { //Read first 3 lines
    fgets(bla,500,fin);
  }
  
  //Basic settings
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin);  sscanf(bla,"%lg",&tmpdbl);
  iter = (int)tmpdbl;
  fgets(bla,500,fin);  sscanf(bla,"%d",&skip);
  fgets(bla,500,fin);  sscanf(bla,"%d",&screenoutput);
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->mcmcseed);
  fgets(bla,500,fin);  sscanf(bla,"%d",&inject);
  for(i=0;i<npar;i++)  fscanf(fin,"%d",&run->setranpar[i]);  //Read the array directly, because sscanf cannot be in a loop...
  fgets(bla,500,fin);  //Read the rest of this line
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->ranparseed);
  fgets(bla,500,fin);  sscanf(bla,"%d",&adapt);
  fgets(bla,500,fin);  sscanf(bla,"%lf",&run->blockfrac);
  for(i=0;i<npar;i++)  fscanf(fin,"%d",&fitpar[i]);  //Read the array directly, because sscanf cannot be in a loop...
  fgets(bla,500,fin);  //Read the rest of this line
  
  
  //Start from offset values:
  fgets(bla,500,fin);  fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin);  sscanf(bla,"%d",&offsetmcmc);
  fgets(bla,500,fin);  sscanf(bla,"%lf",&offsetx);
  for(i=0;i<npar;i++)  fscanf(fin,"%d",&offsetpar[i]);  //Read the array directly, because sscanf cannot be in a loop...
  fgets(bla,500,fin);  //Read the rest of this line
  
  
  //Correlated update proposals:
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin);  sscanf(bla,"%d",&corrupd);
  fgets(bla,500,fin);  sscanf(bla,"%lf",&run->corrfrac);
  fgets(bla,500,fin);  sscanf(bla,"%lg",&tmpdbl);
  ncorr = (int)tmpdbl;
  fgets(bla,500,fin);  sscanf(bla,"%lf",&run->mataccfr);
  fgets(bla,500,fin);  sscanf(bla,"%d",&prmatrixinfo);
  
  
  //Annealing:
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin);  sscanf(bla,"%lf",&temp0);
  fgets(bla,500,fin);  sscanf(bla,"%lg",&tmpdbl);
  nburn = (int)tmpdbl;
  fgets(bla,500,fin);  sscanf(bla,"%lg",&tmpdbl);
  nburn0 = (int)tmpdbl;
  
  //Parallel tempering:
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin);  sscanf(bla,"%d",&partemp);
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->ntemps);
  fgets(bla,500,fin);  sscanf(bla,"%lf",&tempmax);
  fgets(bla,500,fin);  sscanf(bla,"%d",&savehotchains);
  fgets(bla,500,fin);  sscanf(bla,"%d",&prpartempinfo);
  
  //Output:
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin);  sscanf(bla,"%d",&dosnr);
  fgets(bla,500,fin);  sscanf(bla,"%d",&domcmc);
  fgets(bla,500,fin);  sscanf(bla,"%d",&domatch);
  fgets(bla,500,fin);  sscanf(bla,"%d",&intscrout);
  fgets(bla,500,fin);  sscanf(bla,"%d",&writesignal);
  fgets(bla,500,fin);  sscanf(bla,"%d",&printmuch);
  
  //Diverse:
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->networksize);
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->selectdata);
  fgets(bla,500,fin);  sscanf(bla,"%lf",&downsamplefactor);
  fgets(bla,500,fin);  sscanf(bla,"%lf",&run->targetsnr);
  
  //True parameter values:
  fgets(bla,500,fin); fgets(bla,500,fin); fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment lines
  for(i=0;i<npar;i++) fscanf(fin,"%lf",&truepar[i]);  //Read the array directly, because sscanf cannot be in a loop...
  prior_tc_mean = truepar[2];   //prior_tc_mean is used everywhere
  
  //Typical PDF widths:
  fgets(bla,500,fin); fgets(bla,500,fin); fgets(bla,500,fin); //Read the empty and comment line
  for(i=0;i<npar;i++) fscanf(fin,"%lf",&pdfsigs[i]);  //Read the array directly, because sscanf cannot be in a loop...
  
  //Manual temperatures for parallel tempering:
  fgets(bla,500,fin); fgets(bla,500,fin); fgets(bla,500,fin); //Read the empty and comment line
  for(i=0;i<npar;i++) fscanf(fin,"%lf",&run->temps[i]);  //Read the array directly, because sscanf cannot be in a loop...
  
  fclose(fin);
}







// Write a copy of the input file.
// This provides a nicely formatted copy, which may later be used to start a follow-up run.
// Try to keep this in sync with readinputfile() above.
void writeinputfile(struct runpar *run)
{
  int i;
  FILE *fout;
  sprintf(run->outfilename,"%s.%6.6d",run->infilename,run->mcmcseed);
  if((fout = fopen(run->outfilename,"w")) == NULL) {
    printf("   Could not create file: %s, aborting.\n\n\n",run->outfilename);
    exit(1);
  }
  
  fprintf(fout, "  #Input file for spinning MCMC code.  The LINE NUMBER for each parameter should not change!!!\n\n");
  fprintf(fout, "  %-25s  %-18s  %-s\n","#Value:","Variable:","Description:");
  
  
  fprintf(fout, "\n  #Basic settings:\n");
  fprintf(fout, "  %-25.3g  %-18s  %-s\n",    (double)iter,  "iter",           "Total number of iterations to be computed (e.g. 1e7).");
  fprintf(fout, "  %-25d  %-18s  %-s\n",      skip,          "skip",           "Number of iterations to be skipped between stored steps (100 for 1d).");
  fprintf(fout, "  %-25d  %-18s  %-s\n",      screenoutput,  "screenoutput",   "Number of iterations between screen outputs im the MCMC (1000 for 1d).");
  fprintf(fout, "  %-25d  %-18s  %-s\n",      run->mcmcseed, "mcmcseed",       "Random number seed to start the MCMC: 0-let system clock determine seed.");
  fprintf(fout, "  %-25d  %-18s  %-s\n",      inject,        "inject",         "Inject a signal (1) or not (0).");
  fprintf(fout, " ");
  for(i=0;i<npar;i++) fprintf(fout, "%2d",    run->setranpar[i]);
  fprintf(fout, "    %-18s  %-s\n",                          "setranpar[12]",  "Parameters you want to randomise before injecting the signal; 0: use the value in truepar below, 1: randomise.  These are the same parameters as trueval (ie M1, M2, etc!)");
  fprintf(fout, "  %-25d  %-18s  %-s\n",      run->ranparseed,"ranparseed",    "Random number seed for random injection parameters. Don't change between serial chains of the same run!");
  fprintf(fout, "  %-25d  %-18s  %-s\n",      adapt,         "adapt",          "Use adaptation: 0-no, 1-yes.");
  fprintf(fout, "  %-25.2f  %-18s  %-s\n",    run->blockfrac,"blockfrac",      "Fraction of uncorrelated updates that is updated as a block of all parameters (<=0.0: none, >=1.0: all).");
  fprintf(fout, " ");
  for(i=0;i<npar;i++) fprintf(fout, "%2d",    fitpar[i]);
  fprintf(fout, "    %-18s  %-s\n",                          "fitpar[12]",     "Parameters you want to fit for.");
  
  
  fprintf(fout, "\n  #Start from offset values:\n");
  fprintf(fout, "  %-25d  %-18s  %-s\n",      offsetmcmc,    "offsetmcmc",     "Start the MCMC with offset initial parameters: 0-no, 1-yes.  The exact parameters to be offset are determined in offsetpar below.");
  fprintf(fout, "  %-25.1f  %-18s  %-s\n",    offsetx,       "offsetx",        "Start the MCMC with an offset of x times the typical pdf sigma.");
  fprintf(fout, " ");
  for(i=0;i<npar;i++) fprintf(fout, "%2d",    offsetpar[i]);
  fprintf(fout, "    %-18s  %-s\n",                          "offsetpar[12]",  "Parameters you want to start from offset values. At the moment only works if parameter is also 'fit' (i.e. value is 1 in fitpar).");
  
  
  fprintf(fout, "\n  #Correlated update proposals:\n");
  fprintf(fout, "  %-25d  %-18s  %-s\n",      corrupd,       "corrupd",        "Do correlated update proposals: 0-no, 1-yes but update the matrix only once, 2-yes and update the matrix every ncorr iterations.");
  fprintf(fout, "  %-25.1f  %-18s  %-s\n",    run->corrfrac, "corrfrac",       "Fraction of update proposals that is correlated (0.0-1.0, ~0.7 seems ok). corrupd must be 2. Should this replace corrupd?");
  fprintf(fout, "  %-25.3g  %-18s  %-s\n",    (double)ncorr, "ncorr",          "Number of iterations for which the covariance matrix is calculated.");
  fprintf(fout, "  %-25.1f  %-18s  %-s\n",    run->mataccfr, "mataccfr",       "Fraction of elements on the diagonal that must 'improve' in order to accept a new covariance matrix. ???~0.6-0.8 for unimodal, 0.0-0.2 for multimodal???");
  fprintf(fout, "  %-25d  %-18s  %-s\n",      prmatrixinfo,  "prmatrixinfo",   "Print information to screen on proposed matrix updates: 0-none, 1-some (default), 2-add the old and new matrices.");
  
  
  fprintf(fout, "\n  #Annealing:\n");
  fprintf(fout, "  %-25.2f  %-18s  %-s\n",    temp0,         "temp0",          "Starting temperature of the chain, e.g. 100.0. Set 1.0 for no temperature effect.");
  fprintf(fout, "  %-25.3g  %-18s  %-s\n",    (double)nburn, "nburn",          "Number of iterations for the burn-in phase (1e4) at this number, the temperature drops to 1.0.");
  fprintf(fout, "  %-25.3g  %-18s  %-s\n",    (double)nburn0,"nburn0",         "Number of iterations during which temp=temp0 (e.g. 0.1*nburn, should be lower than ~0.9*nburn).");
  
  
  fprintf(fout, "\n  #Parallel tempering:\n");
  fprintf(fout, "  %-25d  %-18s  %-s\n",      partemp,       "partemp",        "Use parallel tempering:  0-no,  1-auto, fixed T ladder,  2-auto, sinusoid T ladder,  3-manual, fixed T ladder,  4-manual, sinusoid T ladder.  For a manual ladder, see near the bottom of the file.");
  fprintf(fout, "  %-25d  %-18s  %-s\n",      run->ntemps,   "ntemps",         "Number of steps in the temperature ladder for parallel tempering, typically 5-10.");
  fprintf(fout, "  %-25.1f  %-18s  %-s\n",    tempmax,       "tempmax",        "Maximum temperature in automatic parallel-tempering ladder (equidistant in log(T)), typically 20-100, e.g. 50.");
  fprintf(fout, "  %-25d  %-18s  %-s\n",      savehotchains, "savehotchains",  "Save hot (T>1) parallel-tempering chains: 0-no (just the T=1 chain), >0-yes; for every saved T=1 point, save every savehotchains-th hot point.");
  fprintf(fout, "  %-25d  %-18s  %-s\n",      prpartempinfo, "prpartempinfo",  "Print information to screen on the temperature chains: 0-none, 1-some ladder info (default), 2-add chain-swap matrix.");
  
  
  fprintf(fout, "\n  #Output:\n");
  fprintf(fout, "  %-25d  %-18s  %-s\n",      dosnr,         "dosnr",          "Calculate the SNR: 0-no, 1-yes.");
  fprintf(fout, "  %-25d  %-18s  %-s\n",      domcmc,        "domcmc",         "Do MCMC: 0-no, 1-yes.");
  fprintf(fout, "  %-25d  %-18s  %-s\n",      domatch,       "domatch",        "Calculate matches: 0-no, 1-yes.");
  fprintf(fout, "  %-25d  %-18s  %-s\n",      intscrout,     "intscrout",      "Print initialisation output to screen: 0-no, 1-yes.");
  fprintf(fout, "  %-25d  %-18s  %-s\n",      writesignal,   "writesignal",    "Write signal, noise, PSDs to file: 0-no, 1-yes.");
  fprintf(fout, "  %-25d  %-18s  %-s\n",      printmuch,     "printmuch",      "Print long stretches of output (1) or not (0).");
  
  
  fprintf(fout, "\n  #Diverse:\n");
  fprintf(fout, "  %-25d  %-18s  %-s\n",    run->networksize,"networksize",    "Set the number of detectors that make up the network: 1: H1, 2: H1L1, 3: H1L1V");
  fprintf(fout, "  %-25d  %-18s  %-s\n",     run->selectdata,"selectdata",     "Select the data set to run on  (set to 0 to print a list of data sets). Make sure you set the true tc and datadir accordingly.");
  fprintf(fout, "  %-25.1f  %-18s  %-s\n", downsamplefactor,"downsamplefactor","Downsample the sampling frequency of the detector (16384 or 20000 Hz) by this factor. Default: 4.0. 10+1.4Mo needs ~16x a<0.1, 8x: a<=0.8, 4x: a>0.8");
  fprintf(fout, "  %-25.1f  %-18s  %-s\n",   run->targetsnr, "targetsnr",      "If > 0: scale the distance such that the network SNR becomes targetsnr");
  //fprintf(fout, "  %-25.1f  %-18s  %-s\n",   cutoff_a,"cutoff_a","Low value of a/M where signal should be cut off, e.g. 7.5.");
  
  
  fprintf(fout, "\n");
  fprintf(fout, "\n  #True parameter values, *not* the exact MCMC parameters and units!  These are used to inject a signal and/or start the MCMC from.\n");
  fprintf(fout, "   M1(Mo)    M2(Mo)            t_c (GPS)   d_L(Mpc)    a_spin    th_SL(d)    R.A.(h)     dec(d)    phic(d)   th_Jo(d)   phi_Jo(d)  alpha(d)  \n");
  for(i=0;i<npar;i++) {
    if(i==2) {
      fprintf(fout, "  %-18.6lf",truepar[i]);
    } else {
      fprintf(fout, "  %-9.4lf",truepar[i]);
    }
  }
  
  fprintf(fout, "\n");
  fprintf(fout, "\n  #Typical PDF widths (used for first correlation matrix and offset run):\n");
  for(i=0;i<npar;i++) fprintf(fout, "  %-7.4f",pdfsigs[i]);
  
  fprintf(fout, "\n");
  fprintf(fout, "\n  #Manual temperature ladder for parallel tempering:\n");
  for(i=0;i<run->ntemps;i++) fprintf(fout, "  %-7.2f",run->temps[i]);
  
  
  fprintf(fout, "\n");
  
  
  /*
  Formats used:
  fprintf(fout, "\n  #:\n");
  fprintf(fout, "  %-25d  %-18s  %-s\n",    ,"","");
  fprintf(fout, "  %-25.1f  %-18s  %-s\n",   ,"","");
  fprintf(fout, "  %-25.1e  %-18s  %-s\n",  ,"","");
  */
  
  fprintf(fout, "\n\n\n");
  fclose(fout);
}





//Read the input file for local (system-dependent) variables: mcmc.local
void readlocalfile()
{
  int i;
  char localfilename[99], bla[500];
  FILE *fin;
  
  sprintf(localfilename,"mcmc.local");
  if((fin = fopen(localfilename,"r")) == NULL) {
    printf("   Error reading file: %s, aborting.\n\n\n",localfilename);
    exit(1);
  }
  
  //Use and l for floats: %lf, %lg, etc, rather than %f, %g
  for(i=1;i<=3;i++) { //Read first 3 lines
    fgets(bla,500,fin);
  }  

  //Data directory:
  fscanf(fin, "%s",datadir);
  
  fclose(fin);
}










// Set the global variables.
// Many of these are now in the input file or unused.
// This routine should eventually disappear.
void setconstants(struct runpar *run)
{
  int i=0,j=0;
  npar           =  12;            // Number of parameters, not *really* a variable... (yet anyway)
  
  i = run->mcmcseed;  //Keeps compiler from complaining:  parameter "run" was never referenced
  j = i;              //Keeps compiler from complaining:  variable "i" was set but never used
  i = j;
  
  /*
  iter           =  1e1;           // Total number of iterations to be computed (~1e7)
  skip           =  10;             // Number of iterations to be skipped between stored steps (100 for 1d)
  screenoutput   =  10;            // Number of iterations between screen outputs im the MCMC (1000 for 1d)
  run->mcmcseed  =  0;             // Random number seed to start the MCMC 0-let system clock determine seed
  adapt          =  1;             // Use adaptation: 0-no, 1-yes
  run->blockfrac   =  0.1;           // Fraction of uncorrelated updates that is updated as a block of all parameters (<=0.0: none, >=1.0: all)
  
  

  //Start from offset values:
  offsetmcmc     =  0;             // Start the MCMC with offset initial parameters: 0-no, 1-yes.  The exact parameters to be offset are determined in offsetpar below.
  offsetx     =  1.e1;          // Start the MCMC with an offset of x times the typical pdf sigma
  //offsetpar[12]        =  {1,1,1,1, 1,1,1,1, 1,1,1,1};  //Parameters you want to start from offset values at the moment only works if parameter is also 'fit' (i.e. value is 1 in fitpar)
  
  for(i=0;i<12;i++) offsetpar[i] = 1;
  
  //Correlated update proposals:
  corrupd        =  2;             // Do correlated update proposals: 0-no, 1-yes but update the matrix only once, 2-yes and update the matrix every ncorr iterations
  ncorr          =  1e2;           // Number of iterations for which the covariance matrix is calculated
  prmatrixinfo   =  1;             // Print information to screen on proposed matrix updates: 0-none, 1-some (default), 2-a lot
  
  //Annealing:
  temp0       =  1.e0;          // Starting temperature of the chain set 1.0 for no temperature effect
  nburn          =  1e5;           // Number of iterations for the burn-in phase (1e4) at this number, the temperature drops to 1.0
  nburn0         =  1e5;           // Number of iterations during which temp=temp0 (e.g. 0.1*nburn, should be lower than ~0.9*nburn)
  
  //Parallel tempering:
  partemp        =  2;             // Use parallel tempering:  0-no,  1-auto, fixed T ladder,  2-auto, sinusoid T ladder,  3-manual, fixed T ladder,  4-manual, sinusoid T ladder
  ntemps         =  3;             // Number of steps in the temperature ladder for parallel tempering
  tempmax     =  5.e1;          // Maximum temperature in automatic parallel-tempering ladder (equidistant in log(T))
  savehotchains  =  0;             // Save all parallel-tempering chains: 0-no, 1-yes
  prpartempinfo  =  2;             // Print information to screen on the temperature chains: 0-none, 1-some ladder info (default), 2-add chain-swap matrix
  
  //Output:
  dosnr          =  1;             // Calculate the SNR: 0-no, 1-yes
  domcmc         =  1;             // Do MCMC: 0-no, 1-yes
  domatch        =  0;             // Calculate matches: 0-no, 1-yes
  intscrout      =  0;             // Print initialisation output to screen: 0-no, 1-yes
  writesignal    =  0;             // Write signal, noise, PSDs to file: 0-no, 1-yes
  printmuch            =  0;             // Print long stretches of output (1) or not (0), global variable, so that you can decide to switch it on/off locally
  
  
  truespin    =  0.8;           // True value of spin
  truetheta   =  55.0;           // True value of theta_SL
  //cutoff_a    =  7.5;           // Low value of a/M where signal should be cut off
  prior_tc_mean = 700009012.345;  //Simulated data, the one MvdS used in the beginning
  //prior_tc_mean = 839366100.345;  //Real (S5) data 1
  //prior_tc_mean = 846226100.345;  //Real (S5) data 2
  downsamplefactor = 16.0;       //Downsample the sampling frequency of the detector (16384 or 20000 Hz) by this factor. Default: 4.0, shouldn't be higher than 8...
  
  //0:mc, 1:eta, 2:tc, 3:logdl,     4:spin, 5:kappa, 6: longi (->RA), 7:sinlati,      8:phase, 9:sinthJ0, 10:phiJ0, 11:alpha
  //0,1,2,3, 4,5,6,7, 8,9,0,1
  //Parameters you want to fit
  //fitpar[12] = {1,1,1,1, 1,1,1,1, 1,1,1,1};  //All
  for(i=0;i<12;i++) fitpar[i] = 1;
  //fitpar[12] = {0,0,1,1, 0,0,1,1, 0,0,0,0};  //t_c, dist and sky position
  //fitpar[12] = {0,0,1,1, 0,0,1,1, 1,1,1,1};  //Extrinsic parameters
  //fitpar[12] = {1,1,0,0, 1,1,0,0, 0,0,0,0};  //Intrinsic parameters
  
  //pdfsigs[12] = {0.0233,0.0252,0.008,0.260,  0.1129,0.16,0.57,0.139,  1.5,0.135,0.282,1.5};  //Typical width of a PDF, for 1D, correlated updates, a=0.8: 90% range
  pdfsigs[0]  = 0.0233;
  pdfsigs[1]  = 0.0252;
  pdfsigs[2]  = 0.008;
  pdfsigs[3]  = 0.260;
  pdfsigs[4]  = 0.1129;
  pdfsigs[5]  = 0.16;
  pdfsigs[6]  = 0.57;
  pdfsigs[7]  = 0.139;
  pdfsigs[8]  = 1.5;
  pdfsigs[9]  = 0.135;
  pdfsigs[10] = 0.282;
  pdfsigs[11] = 1.5;
  //pdfsigs[12] = {0.0233,0.0252,0.05,0.30,  0.1129,0.16,1.3,0.4,  1.5,0.135,0.282,1.5};  //Test run with only tc, d, RA and dec.
  
  //datadir[] = "/home/sluys/work/GW/programs/MCMC/data"; // Directory where the datafiles are located: Mac, Fugu
  sprintf(datadir,"/home/sluys/work/GW/programs/MCMC/data");
  //const char datadir[] = "/guest/sluys/work/GW/programs/MCMC/data"; // Directory where the datafiles are located: Tsunami
  //const char datadir[] = "/attic/sluys/work/GW/programs/MCMC/data"; // Directory where the datafiles are located: Typhoon

  
  */
  
  /*
  //Not in MCMC units!
  truepar[0]  = 10.0;              // M1
  truepar[1]  = 1.4;               // M2
  truepar[2]  = 700009012.345000;  // tc
  truepar[3]  = 13.0;              // d_L
  truepar[4]  = 0.8;               // spin
  truepar[5]  = 55.0;              // theta (deg)
  truepar[6]  = 91.0;              // 'Greenwich hour angle' (saved as RA) (deg)
  truepar[7]  = 40.0;              // declination (deg)
  truepar[8]  = 11.459156;         // phase (deg)
  truepar[9]  = 15.0;              // theta_J0 (deg)
  truepar[10] = 125.0;             // phiJ0 (deg)
  truepar[11] = 51.566202;         // alpha (deg)
  */
  
  //--- Program control: --------------------------------------------------------------------------------------
  tempi = 0; //A global variable that determines the current chain (temperature) in the temperature ladder
  
  
  
  
  
  /*  note that since `complex.h' is included BEFORE `fftw3.h',   */
  /*  the data type `fftw_complex' equals `double complex'        */
  /*  (see: http://www.fftw.org/fftw3_doc/Complex-numbers.html).  */
  
  /*-- define some physical constants --*/
  /* Nelson Christensen wrote (18/03/'05):                                */
  /* > Christian and Mark                                                 */
  /* > For the data we are using we need to use the following constants:  */
  Ms   = 1.9889194662e30;   /* solar mass (kg)                                            */ 
  Mpc  = 3.08568025e22;     /* metres in a Mpc  (LAL: 3.0856775807e22)                    */
  G    = 6.67259e-11;       /* 6.674215e-11; */ /* gravity constant (SI)                  */
  c    = 299792458.0;       /* speed of light (m/s)                                       */
  Mpcs = 1.029272137e14;    /* seconds in a Mpc  (Mpc/c)                                  */
  pi   = 3.141592653589793;   // pi
  tpi  = 6.283185307179586;   // 2 pi
  mtpi = 6.283185307179586e6; // Large multiple of 2 pi
  
  
  
  /*-- set some prior parameters: --*/
  //prior_tc_margin  = 0.005;/*0.005;*/  /* margin around `prior_tc_mean' (seconds, see below) */
  //prior_mass_lower =   1.0;/*0.9;*/  /* Ms  */
  //prior_mass_upper =  10.0; /*2.1;*/
  //prior_dist_90    =  50.0; /* distance at which a 2-2Ms inspiral has 90% detection prob. */
  //prior_dist_10    =  60.0; /* distance at which a 2-2Ms inspiral has 10% detection prob. */
  
  /*--- logfile for MCMC chain: --------------------------------------------------------------------------*/
  parallelchains =    2; /* number of parallel chains                                        */
  impodraws      =   10; /* number of importance-resampled starting draws                    */
  logpostdiff    = 20.0; /* max. difference of log-post.-densities allowed in startsample    */
  anneal         =   10; /* number of annealing iterations                                   */
  annealtemp     =  1.0; /* `temperature' from which (variance) annealing starts             */
  covest         = 29000; /* starting iteration for covariance estimation (9900000)          */
  covfix         = 29500; /* number of iterations until covariance is fixed (no more updates) (9950000) */
  covskip        = 100; /* iterations to be skipped between covariance updates (100)        */
  initweight     = 1000; /* weight (in iterations) of initial covariance in later updates (1000000)   */
  studentDF      =  3.0; /* degrees of freedom for Student-t jumps (0.0 indicates Normal)    */
  modifiedStudent =   1; /* modify the (student-t) proposal                                  */
  propscale      = 0.2 ; /* scale of proposal covariance relative to assumed posterior cov.  */
  tukeywin       = 0.05; /* parameter for Tukey-window used in dataFT (non-flat fraction).   Use 0.15 for Virgo data  */
  randomseed1    = 3335; /* random seed                                                      */
  randomseed2    = 4449; /*   "     "                                                        */
  unifMcJump     = 0.10; /* probability for uniform chirp mass  proposals                    */
  unifTcJump     = 0.10; /* probability for uniform coalescence time proposals               */
  unifDirectJump = 0.10; /* probability for uniform direction proposals                      */
  invDirectJump  = 0.01; /* probability for inverse direction proposals                      */
  wideDirectJump = 0.20; /* probability for 10 times as wide jumps                           */
  unifOrientJump = 0.10; /* probability for uniform orientation proposals                    */
  invPhaseJump   = 0.05; /* probability for inverse phase proposals  (phi + pi)              */
  invIncliJump   = 0.05; /* probability for inverse inclination proposals  (pi - iota)       */
  iter_per_min   = 60.0; /* (a guess in order to estimate the computation time)              */
  mutationprob   = 0.75; /* mutation probability (remaining moves are crossovers)            */
  
  //In mcmc.input now
  //inject            = 1; /* inject `homemade' signal?                                        */
  
  annealfact = 2.0;      /* temperature increase factor for subsequent chains                */
  
}



void setrandomtrueparameters(struct runpar *run)  //Get random values for the 'true' parameters for the 12-parameter spinning template.  
// *** This changes the injected signal!!! ***
{
  int i=0;
  gsl_rng *ran;
  ran = gsl_rng_alloc(gsl_rng_mt19937);  // GSL random-number seed
  if(1==2) {  //Select a random seed, *** ONLY FOR TESTING ***
    printf("\n  *** SELECTING RANDOM SEED ***  This should only be done while testing!!! setrandumtrueparameters() \n\n");
    run->ranparseed = 0;
    setseed(&run->ranparseed);
    //printf("  Seed: %d\n", run->ranparseed);
  }
  gsl_rng_set(ran, run->ranparseed);     // Set seed
  
  //Lower and upper boundaries:
  double *lb,*ub,db,dt;
  lb = (double*)calloc(12,sizeof(double));
  ub = (double*)calloc(12,sizeof(double));
  lb[0] = 5.0;       //M1 (Mo)
  ub[0] = 15.0;
  lb[1] = 1.2;       //M2 (Mo)
  ub[1] = 1.6;
  dt = 0.5; //This is dt/2
  lb[2] = prior_tc_mean - dt; //t_c
  ub[2] = prior_tc_mean + dt;
  lb[3] = 10.0;      //d_L (Mpc)  !Linear!
  ub[3] = 30.0;
  lb[4] = 1.e-10;    //a_spin (0-1)
  ub[4] = 0.999999;
  lb[5] = 0.001;     //th_SL (deg) (not used)
  ub[5] = 179.999;
  lb[6] = 0.0;       //RA (h)
  ub[6] = 24.0;
  lb[7] = -89.999;   //dec (deg) (not used)
  ub[7] = 89.999;
  lb[8] = 0.0;       //phi_c (deg)
  ub[8] = 360.0;
  lb[9] = -89.999;   //theta_J0 (deg) (not used)
  ub[9] = 89.999;
  lb[10] = 0.0;      //phi_Jo (deg)
  ub[10] = 360.0;
  lb[11] = 0.0;      //alpha_c (deg)
  ub[11] = 360.0;
  
  for(i=0;i<npar;i++) {
    db = ub[i]-lb[i];
    if(run->setranpar[i]==1) truepar[i] = gsl_rng_uniform(ran)*db + lb[i];
    if(i==5 && run->setranpar[i]==1) truepar[i] = acos(gsl_rng_uniform(ran)*2.0 - 1.0)*r2d;  //kappa -> th_SL
    if((i==7 || i==9)  && run->setranpar[i]==1) truepar[i] = asin(gsl_rng_uniform(ran)*2.0 - 1.0)*r2d;  //sin(dec)->dec, sin(th_J0)->th_J0
    //printf("  %d  %lf  %lf  %lf  %lf\n",i,lb[i],ub[i],db,truepar[i]);
  }
  
  free(lb);
  free(ub);
  gsl_rng_free(ran);
}




void gettrueparameters(struct parset *par)  //Set the parameters for the 12-parameter spinning template to the 'true values'
{
  par->m1       = truepar[0];                    // M1 (10.0)
  par->m2       = truepar[1];                    // M2  (1.4)
  par->m        = par->m1+par->m2;
  par->mu       = par->m1*par->m2/par->m;
  par->eta      = par->mu/par->m;                // mass ratio                
  par->mc       = par->m*pow(par->eta,0.6);      // chirp mass. in Mo         
  par->tc       = truepar[2];                    // coalescence time
  par->logdl    = log(truepar[3]);               // log-distance (Mpc) (17.5)             
  
  par->spin     = truepar[4];                    // magnitude of total spin   (0.1)
  par->kappa    = cos(truepar[5]*d2r);           // L^.S^, cos of angle between L^ & S^  (0.819152)
  //par->longi    = truepar[6]*d2r;                // longitude (~'Greenwich hour angle', saved as RA)     (120)
  par->longi    = fmod(longitude(truepar[6]*h2r,GMST(par->tc))+mtpi,tpi);  //The parameter in the input and output is RA; the MCMC parameter is 'longi' ~ Greenwich hour angle
  par->sinlati  = sin(truepar[7]*d2r);           // sin latitude (sin(delta))  (40)     
  
  par->phase    = truepar[8]*d2r;                // orbital phase   (phi_c)   (0.2)
  par->sinthJ0  = sin(truepar[9]*d2r);           // sin Theta_J0 ~ latitude, pi/2 = NP    (15)
  par->phiJ0    = truepar[10]*d2r;               // Phi_J0 ~ azimuthal            (125)
  par->alpha    = truepar[11]*d2r;               // Alpha_c                       (0.9 rad = 51.566202deg)
  
  par->loctc    = NULL;
  par->localti  = NULL;
  par->locazi   = NULL;
  par->locpolar = NULL;
}




void getnullparameters(struct parset *par)  //Set the parameters for the 12-parameter spinning template to 'null values', to simulate absence of a signal
{
  par->mc       = 0.01;
  par->eta      = 0.0001;
  par->tc       = prior_tc_mean;
  par->logdl    = 10000.0;
  par->sinlati  = 0.0;
  par->longi    = 0.0;
  par->phase    = 0.0;
  par->spin     = 0.000001;              // magnitude of total spin   
  par->kappa    = 0.0;                   // L^.S^, cos of angle between L^ & S^
  par->sinthJ0  = 0.000001;              // sin Theta_J0 ~ latitude, pi/2=NP    
  par->phiJ0    = 0.0;                   // Phi_J0 ~ azimuthal        
  par->alpha    = 0.0;                   // Alpha_c                   
  
  par->loctc    = NULL;
  par->localti  = NULL;
  par->locazi   = NULL;
  par->locpolar = NULL;
}








