// mcmc_parameters.c:
// Routines to read/write input files, set constants and set true and null parameters

#include <mcmc.h>





// Read the main input file.
// Please keep this routine in sync with writeinputfile() below.
// All parameters that are read in here should become members of the runvar struct and lose their global status
void readMainInputfile(struct runpar *run)
{
  int i;
  char bla[500];
  FILE *fin;
  
  if((fin = fopen(run->infilename,"r")) == NULL) {
    printf("   Error reading main input file: %s, aborting.\n\n\n",run->infilename);
    exit(1);
  } else {
    printf("   Using main input file: %s.\n",run->infilename);
  }
  
  
  //Use and l for floats: %lf, %lg, etc, rather than %f, %g
  
  for(i=1;i<=3;i++) fgets(bla,500,fin);  //Read first 3 lines
  
  //Operation and output:
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin);  sscanf(bla,"%d",&doSNR);
  fgets(bla,500,fin);  sscanf(bla,"%d",&doMCMC);
  fgets(bla,500,fin);  sscanf(bla,"%d",&doMatch);
  fgets(bla,500,fin);  sscanf(bla,"%d",&intscrout);
  fgets(bla,500,fin);  sscanf(bla,"%d",&writeSignal);
  fgets(bla,500,fin);  sscanf(bla,"%d",&printMuch);
  
  
  
  
  //Software injection:
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin);  sscanf(bla,"%d",&inject);
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->mcmcWaveform);
  run->injectionWaveform = run->mcmcWaveform;  //For now
  
  if(inject>=1) {
    if(run->injectionWaveform==1) {
      printf("   Using Apostolatos, 1.5PN, 12-parameter waveform for the software injection.\n");
      run->nMCMCpar=12;
    } else if(run->injectionWaveform==2) {
      printf("   Using LAL, 3.5PN, 12-parameter waveform for the software injection.\n");
      run->nMCMCpar=12;
    } else if(run->injectionWaveform==3) {
      printf("   Using LAL, 3.5PN, 15-parameter waveform for the software injection.\n");
      run->nMCMCpar=15;
    } else {
      printf("   Unknown waveform chosen as MCMC template: %d.   Available waveforms are:\n",run->injectionWaveform);
      printf("     1: Apostolatos, simple precession, 12 parameters\n");
      printf("     2: LAL, single spin, 12 parameters\n");
      printf("     3: LAL, double spin, 15 parameters\n");
      printf("   Please set injectionWaveform in mcmc.input to one of these values.\n\n");
      exit(1);
    }
  }
  npar = run->nMCMCpar;
  run->nInjectPar = run->nMCMCpar;
  run->ranInjPar  = (int*)calloc(run->nMCMCpar,sizeof(int)); //Now that run->nMCMCpar is known, initialise the array.... 
  
  fgets(bla,500,fin);  sscanf(bla,"%lf",&run->injectionSNR);
  for(i=0;i<run->nMCMCpar;i++)  fscanf(fin,"%d",&run->ranInjPar[i]);  //Read the array directly, because sscanf cannot be in a loop...
  fgets(bla,500,fin);  //Read the rest of this line
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->ranparseed);
  


  
  //True parameter values:
  fgets(bla,500,fin); fgets(bla,500,fin); fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment lines
  for(i=0;i<run->nMCMCpar;i++) fscanf(fin,"%lf",&truepar[i]);  //Read the array directly, because sscanf cannot be in a loop...
  prior_tc_mean = truepar[2];   //prior_tc_mean is used everywhere
  //Starting parameter values:
  for(i=0;i<run->nMCMCpar;i++) fscanf(fin,"%lf",&run->startpar[i]);  //Read the array directly, because sscanf cannot be in a loop...
  //if(offsetmcmc==0 || offsetmcmc==1) for(i=0;i<run->nMCMCpar;i++) run->startpar[i] = truepar[i];  //Set the starting parameters equal to the true, injection parameters
  
  //Secondary input files:
  fgets(bla,500,fin); fgets(bla,500,fin); fgets(bla,500,fin); fgets(bla,500,fin); //Read the empty and comment line
  fgets(bla,500,fin); sscanf(bla,"%s",run->MCMCinfilename);
  fgets(bla,500,fin); sscanf(bla,"%s",run->datainfilename);
  fgets(bla,500,fin); sscanf(bla,"%s",run->parameterinfilename);
  
  fclose(fin);
}







// Write a copy of the input file.
// This provides a nicely formatted copy, which may later be used to start a follow-up run.
// Try to keep this in sync with readinputfile() above.
void writeMainInputfile(struct runpar *run)
{
  int i;
  FILE *fout;
  sprintf(run->outfilename,"%s.%6.6d",run->infilename,run->MCMCseed);
  if((fout = fopen(run->outfilename,"w")) == NULL) {
    printf("   Could not create file: %s, aborting.\n\n\n",run->outfilename);
    exit(1);
  }
  
  fprintf(fout, "  #Input file for spinning MCMC code.  The LINE NUMBER for each parameter should not change!!!\n\n");
  fprintf(fout, "  %-39s  %-18s  %-s\n","#Value:","Variable:","Description:");
  
  
  fprintf(fout, "  \n  #Basic settings:\n");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      run->mcmcWaveform, "mcmcWaveform", "Waveform version used as MCMC template:  1 for 1.5PN 12-parameter Apostolatos, 2 for 3.5PN 12-parameter LAL, 3 for 3.5PN 15-parameter LAL");
  fprintf(fout, "  %-39.3g  %-18s  %-s\n",    (double)iter,  "iter",           "Total number of iterations to be computed (e.g. 1e7).");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      thinOutput,    "thinOutput",           "Number of iterations to be skipped between stored steps (100 for 1d).");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      thinScreenOutput,  "thinScreenOutput",   "Number of iterations between screen outputs im the MCMC (1000 for 1d).");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      run->MCMCseed, "MCMCseed",       "Random number seed to start the MCMC: 0-let system clock determine seed.");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      inject,        "inject",         "Inject a signal (1) or not (0).");
  fprintf(fout, " ");
  for(i=0;i<run->nMCMCpar;i++) fprintf(fout, "%2d",    run->ranInjPar[i]);
  for(i=0;i<max(19-run->nMCMCpar,0);i++) fprintf(fout, "  ");  //Line up the next colum nicely, for up to 20 parameters
  fprintf(fout, "    %-18s  %-s\n",                          "ranInjPar[]",  "Parameters you want to randomise before injecting the signal; 0: use the value in truepar below, 1: randomise.  These are the same parameters as trueval (ie M1, M2, etc!)");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      run->ranparseed,"ranparseed",    "Random number seed for random injection parameters. Don't change between serial chains of the same run!");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      adapt,         "adapt",          "Use adaptation: 0-no, 1-yes.");
  fprintf(fout, "  %-39.2f  %-18s  %-s\n",    run->blockfrac,"blockfrac",      "Fraction of uncorrelated updates that is updated as a block of all parameters (<=0.0: none, >=1.0: all).");
  fprintf(fout, " ");
  
  
  fprintf(fout, "  \n  #Start from offset values:\n");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      offsetmcmc,    "offsetmcmc",     "Start the MCMC with offset initial parameters: 0-no: use injection parameters, 1-yes: randomly around the injected parameters, 2-yes: at the starting parameters, 3-yes: randomly around the starting parameters.  The individual parameters to be offset are determined in offsetpar below.");
  fprintf(fout, "  %-39.1f  %-18s  %-s\n",    offsetx,       "offsetx",        "Start the MCMC with an offset of x times the typical pdf sigma.");
  fprintf(fout, " ");
  for(i=0;i<run->nMCMCpar;i++) fprintf(fout, "%2d",    offsetpar[i]);
  for(i=0;i<max(19-run->nMCMCpar,0);i++) fprintf(fout, "  ");  //Line up the next colum nicely, for up to 20 parameters
  fprintf(fout, "    %-18s  %-s\n",                          "offsetpar[]",  "Parameters you want to start from random offset values. At the moment only works if parameter is also 'fit' (i.e. value is 0 in parFix).");
  
  
  fprintf(fout, "  \n  #Correlated update proposals:\n");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      corrupd,       "corrupd",        "Do correlated update proposals: 0-no, 1-yes but update the matrix only once, 2-yes and update the matrix every ncorr iterations.");
  fprintf(fout, "  %-39.1f  %-18s  %-s\n",    run->corrfrac, "corrfrac",       "Fraction of update proposals that is correlated (0.0-1.0, ~0.7 seems ok). corrupd must be 2. Should this replace corrupd?");
  fprintf(fout, "  %-39.3g  %-18s  %-s\n",    (double)ncorr, "ncorr",          "Number of iterations for which the covariance matrix is calculated.");
  fprintf(fout, "  %-39.1f  %-18s  %-s\n",    run->mataccfr, "mataccfr",       "Fraction of elements on the diagonal that must 'improve' in order to accept a new covariance matrix. ???~0.6-0.8 for unimodal, 0.0-0.2 for multimodal???");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      prmatrixinfo,  "prmatrixinfo",   "Print information to screen on proposed matrix updates: 0-none, 1-some (default), 2-add the old and new matrices.");
  
  
  fprintf(fout, "  \n  #Annealing:\n");
  fprintf(fout, "  %-39.2f  %-18s  %-s\n",    temp0,         "temp0",          "Starting temperature of the chain, e.g. 100.0. Set 1.0 for no temperature effect.");
  fprintf(fout, "  %-39.3g  %-18s  %-s\n",    (double)nburn, "nburn",          "Number of iterations for the burn-in phase (1e4) at this number, the temperature drops to 1.0.");
  fprintf(fout, "  %-39.3g  %-18s  %-s\n",    (double)nburn0,"nburn0",         "Number of iterations during which temp=temp0 (e.g. 0.1*nburn, should be lower than ~0.9*nburn).");
  
  
  fprintf(fout, "  \n  #Parallel tempering:\n");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      partemp,       "partemp",        "Use parallel tempering:  0-no,  1-auto, fixed T ladder,  2-auto, sinusoid T ladder,  3-manual, fixed T ladder,  4-manual, sinusoid T ladder.  For a manual ladder, see near the bottom of the file.");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      run->ntemps,   "ntemps",         "Number of steps in the temperature ladder for parallel tempering, typically 5-10.");
  fprintf(fout, "  %-39.1f  %-18s  %-s\n",    tempmax,       "tempmax",        "Maximum temperature in automatic parallel-tempering ladder (equidistant in log(T)), typically 20-100, e.g. 50.");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      savehotchains, "savehotchains",  "Save hot (T>1) parallel-tempering chains: 0-no (just the T=1 chain), >0-yes; for every saved T=1 point, save every savehotchains-th hot point.");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      prpartempinfo, "prpartempinfo",  "Print information to screen on the temperature chains: 0-none, 1-some ladder info (default), 2-add chain-swap matrix.");
  
  
  fprintf(fout, "  \n  #Output:\n");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      doSNR,         "doSNR",          "Calculate the SNR: 0-no, 1-yes.");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      doMCMC,        "doMCMC",         "Do MCMC: 0-no, 1-yes.");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      doMatch,       "doMatch",        "Calculate matches: 0-no, 1-yes.");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      intscrout,     "intscrout",      "Print initialisation output to screen: 0-no, 1-yes.");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      writeSignal,   "writeSignal",    "Write signal, noise, PSDs to file: 0-no, 1-yes.");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      printMuch,     "printMuch",      "Print long stretches of output (1) or not (0).");
  

  fprintf(fout, "  \n  #Data handling:\n");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      run->selectdata,        "selectdata",       "Select the data set to run on:  0 use input file  (set to -1 to print a list of data sets). Make sure you set the true tc and datadir accordingly.");
  fprintf(fout, "  %-39d  %-18s  %-s\n",      downsamplefactor,       "downsamplefactor", "Downsample the sampling frequency of the detector (16-20kHz for the detectors, 4kHz for NINJA) by this factor.  Default (for detectors): 4.0. 10+1.4Mo needs ~16x a<0.1, 8x: a<=0.8, 4x: a>0.8");
  fprintf(fout, "  %-39.1f  %-18s  %-s\n",    run->databeforetc,      "databeforetc",     "The stretch of data in seconds before presumed coalescence that is read in as part of the data segment");
  fprintf(fout, "  %-39.1f  %-18s  %-s\n",    run->dataaftertc,       "dataaftertc",      "The stretch of data in seconds after presumed coalescence that is read in as part of the data segment");
  fprintf(fout, "  %-39.1f  %-18s  %-s\n",    run->lowfrequencycut,   "lowfrequencycut",  "Templates and overlap integration start at this frequency");
  fprintf(fout, "  %-39.1f  %-18s  %-s\n",    run->highfrequencycut,  "highfrequencycut", "Overlap integration ends at this frequency");

  
  fprintf(fout, "  \n  #Diverse:\n");
  fprintf(fout, "  %-39d  %-18s  %-s\n",    run->networksize,"networksize",    "Set the number of detectors that make up the network");
  fprintf(fout, "  ");
  for(i=0;i<run->networksize;i++) fprintf(fout, "%-3d",run->selectifos[i]);
  for(i=0;i<38-3*(run->networksize-1);i++) fprintf(fout, " ");
  fprintf(fout, "%-18s  %-s\n", "selectifos",    "Select the IFOs to use  1: H1, 2: L1, 3: V");
  fprintf(fout, "  %-39.1f  %-18s  %-s\n",   run->injectionSNR, "injectionSNR",      "If > 0: scale the distance such that the network SNR becomes injectionSNR");
  


  
  //fprintf(fout, "  %-39.1f  %-18s  %-s\n",   cutoff_a,"cutoff_a","Low value of a/M where signal should be cut off, e.g. 7.5.");
  
  
  fprintf(fout, "  \n");
  fprintf(fout, "  \n  #Injection (first line) and starting (second line) parameter values, these are the exact MCMC parameters and units\n");
  if(run->mcmcWaveform==1) {
    fprintf(fout, "      Mc        eta             t_c(GPS)  log(d_L)     a_spin cos(th_SL)       R.A.   sin(dec)      phi_c   sinth_Jo      phi_Jo     alpha     \n");
  } else if(run->mcmcWaveform==2 || run->mcmcWaveform==3) {
    fprintf(fout, "      Mc        eta             t_0(GPS)   log(dL)       R.A.     sinDec     cos(i)      phi_c        psi    aspin_1 cos(theta1)     phi_1    aspin_2 cos(theta2)     phi_2  \n");
  } else {
    fprintf(fout, "  \n");
  }
  for(i=0;i<run->nMCMCpar;i++) {
    if(i==2) {
      fprintf(fout, "  %-18.6lf",truepar[i]);
    } else {
      fprintf(fout, "  %-9.4lf",truepar[i]);
    }
  }
  fprintf(fout, "  \n");
  for(i=0;i<run->nMCMCpar;i++) {
    if(i==2) {
      fprintf(fout, "  %-18.6lf",run->startpar[i]);
    } else {
      fprintf(fout, "  %-9.4lf",run->startpar[i]);
    }
  }
  fprintf(fout, "\n");
  
  fprintf(fout, "  \n  #Typical PDF widths (used for first correlation matrix and offset run):\n");
  for(i=0;i<run->nMCMCpar;i++) fprintf(fout, "  %-7.4f",pdfsigs[i]);
  
  fprintf(fout, "\n");
  fprintf(fout, "  \n  #Manual temperature ladder for parallel tempering:\n");
  for(i=0;i<run->ntemps;i++) fprintf(fout, "  %-7.2f",run->temps[i]);
  
  
  fprintf(fout, "\n");
  
  fprintf(fout, "\n");
  fprintf(fout, "  \n  #Secondary input files:\n");
  fprintf(fout, "  %-39s  %-18s  %-s\n",      run->datainfilename,       "datainfilename",        "File name of the data/noise input file, e.g. mcmc.data");
  
  /*
  Formats used:
  fprintf(fout, "  \n  #:\n");
  fprintf(fout, "  %-39d  %-18s  %-s\n",    ,"","");
  fprintf(fout, "  %-39.1f  %-18s  %-s\n",   ,"","");
  fprintf(fout, "  %-39.1e  %-18s  %-s\n",  ,"","");
  */
  
  fprintf(fout, "\n\n\n");
  fclose(fout);
}





//Read the input file for local (system-dependent) variables: mcmc.local
void readLocalInputfile()
{
  int i;
  char localfilename[99], bla[500];
  FILE *fin;
  
  sprintf(localfilename,"mcmc.local");
  if((fin = fopen(localfilename,"r")) == NULL) {
    printf("   Error reading local file: %s, aborting.\n\n\n",localfilename);
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







// Read the MCMC input file.
// All parameters that are read in here should become members of the runvar struct and lose their global status
void readMCMCinputfile(struct runpar *run)
{
  int i;
  double tmpdbl;
  char bla[500];
  FILE *fin;
  
  if((fin = fopen(run->MCMCinfilename,"r")) == NULL) {
    printf("   Error reading MCMC input file: %s, aborting.\n\n\n",run->MCMCinfilename);
    exit(1);
  } else {
    printf("   Using MCMC input file: %s.\n",run->MCMCinfilename);
  }
  
  
  //Use and l for floats: %lf, %lg, etc, rather than %f, %g
  
  for(i=1;i<=3;i++) { //Read first 3 lines
    fgets(bla,500,fin);
  }
  
  //Basic settings
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->mcmcWaveform);
  run->injectionWaveform = run->mcmcWaveform;  //For now
  
  
  if(run->mcmcWaveform==1) {
    printf("   Using Apostolatos, 1.5PN, 12-parameter waveform as the MCMC template.\n");
    run->nMCMCpar=12;
  } else if(run->mcmcWaveform==2) {
    printf("   Using LAL, 3.5PN, 12-parameter waveform as the MCMC template.\n");
    run->nMCMCpar=12;
  } else if(run->mcmcWaveform==3) {
    printf("   Using LAL, 3.5PN, 15-parameter waveform as the MCMC template.\n");
    run->nMCMCpar=15;
  } else {
    printf("   Unknown waveform chosen as MCMC template: %d.   Available waveforms are:\n",run->mcmcWaveform);
    printf("     1: Apostolatos, simple precession, 12 parameters\n");
    printf("     2: LAL, single spin, 12 parameters\n");
    printf("     3: LAL, double spin, 15 parameters\n");
    printf("   Please set mcmcWaveform in mcmc.input to one of these values.\n\n");
    exit(1);
  }
  npar = run->nMCMCpar;
  run->nInjectPar = run->nMCMCpar;
  
  
  run->ranInjPar  = (int*)calloc(run->nMCMCpar,sizeof(int)); //Now that run->nMCMCpar is known, initialise the array.... 
  
  fgets(bla,500,fin);  sscanf(bla,"%lg",&tmpdbl);    
  iter = (int)tmpdbl;
  fgets(bla,500,fin);  sscanf(bla,"%d",&thinOutput);
  fgets(bla,500,fin);  sscanf(bla,"%d",&thinScreenOutput);
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->MCMCseed);
  fgets(bla,500,fin);  sscanf(bla,"%d",&adapt);
  fgets(bla,500,fin);  sscanf(bla,"%lf",&run->blockfrac);
  
  

  //Start from offset values:
  fgets(bla,500,fin);  fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin);  sscanf(bla,"%d",&offsetmcmc);
  fgets(bla,500,fin);  sscanf(bla,"%lf",&offsetx);
  for(i=0;i<run->nMCMCpar;i++)  fscanf(fin,"%d",&offsetpar[i]);  //Read the array directly, because sscanf cannot be in a loop...
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
  
  //Typical PDF widths:
  fgets(bla,500,fin); fgets(bla,500,fin); //Read the empty and comment line
  for(i=0;i<run->nMCMCpar;i++) fscanf(fin,"%lf",&pdfsigs[i]);  //Read the array directly, because sscanf cannot be in a loop...

  
  //Manual temperature ladder for parallel tempering:
  fgets(bla,500,fin); fgets(bla,500,fin); fgets(bla,500,fin); //Read the empty and comment line
  for(i=0;i<run->ntemps;i++) fscanf(fin,"%lf",&run->temps[i]);  //Read the array directly, because sscanf cannot be in a loop...
  
  fclose(fin);
}
//End of void readMCMCinputfile(struct runpar *run)












// Read the data input file.
void readDataInputfile(struct runpar *run, struct interferometer ifo[])
{
  int i=0,j=0;
  double lati,longi,leftarm,rightarm;
  char bla[500], subdir[500];
  FILE *fin;
  
  if((fin = fopen(run->datainfilename,"r")) == NULL) {
    printf("   Error reading data file: %s, aborting.\n\n\n",run->datainfilename);
    exit(1);
  }
  
  
  //Use and l for floats: %lf, %lg, etc, rather than %f, %g
  
  for(j=1;j<=3;j++) fgets(bla,500,fin);  //Read first 3 lines
  fgets(run->datasetName,80,fin);  fgets(bla,500,fin);  //Read name of the data set used, and then the rest of the line
  
  //Detector network:
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->networksize);
  for(i=0;i<run->networksize;i++) fscanf(fin,"%d",&run->selectifos[i]);  //Read the array directly, because sscanf cannot be in a loop...
  fgets(bla,500,fin);  //Read the rest of the line
  
  
  //Data handling:
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin); sscanf(bla,"%d",&downsamplefactor);
  fgets(bla,500,fin); sscanf(bla,"%lf",&run->databeforetc);
  fgets(bla,500,fin); sscanf(bla,"%lf",&run->dataaftertc);
  fgets(bla,500,fin); sscanf(bla,"%lf",&run->lowfrequencycut);
  fgets(bla,500,fin); sscanf(bla,"%lf",&run->highfrequencycut);
  fgets(bla,500,fin); sscanf(bla,"%lf",&run->tukeywin);
  
  
  //Read input for PSD estimation:
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment lines
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->PSDsegmentNumber);
  fgets(bla,500,fin);  sscanf(bla,"%lf",&run->PSDsegmentLength);
  
  
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment lines
  for(i=0;i<run->networksize;i++){
    fgets(bla,500,fin); fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment lines
    
    fgets(bla,500,fin);  sscanf(bla,"%s",ifo[i].name);
    fgets(bla,500,fin);  sscanf(bla,"%lf",&lati);
    fgets(bla,500,fin);  sscanf(bla,"%lf",&longi);
    fgets(bla,500,fin);  sscanf(bla,"%lf",&rightarm);
    fgets(bla,500,fin);  sscanf(bla,"%lf",&leftarm);
    
    ifo[i].lati      = lati     /180.0*pi;
    ifo[i].longi     = longi    /180.0*pi;
    ifo[i].rightarm  = rightarm /180.0*pi;
    ifo[i].leftarm   = leftarm  /180.0*pi;
    
    fgets(bla,500,fin);  //Read the empty line
    
    fgets(bla,500,fin);  sscanf(bla,"%s",ifo[i].ch1name);
    //fgets(bla,500,fin);  sscanf(bla,"%s",&ifo[i].ch1filepath);
    fgets(bla,500,fin);  sscanf(bla,"%s",subdir);
    sprintf(ifo[i].ch1filepath,"%s%s%s",datadir,"/",subdir);
    fgets(bla,500,fin);  sscanf(bla,"%s",ifo[i].ch1fileprefix);
    fgets(bla,500,fin);  sscanf(bla,"%s",ifo[i].ch1filesuffix);
    fgets(bla,500,fin);  sscanf(bla,"%d",&ifo[i].ch1filesize);
    fgets(bla,500,fin);  sscanf(bla,"%d",&ifo[i].ch1fileoffset);
    fgets(bla,500,fin);  sscanf(bla,"%d",&ifo[i].ch1doubleprecision);
    fgets(bla,500,fin);  sscanf(bla,"%d",&ifo[i].add2channels);
    
    fgets(bla,500,fin);  //Read the empty line
    
    fgets(bla,500,fin);  sscanf(bla,"%ld",&ifo[i].noiseGPSstart);
    fgets(bla,500,fin);  sscanf(bla,"%s",ifo[i].noisechannel);
    //fgets(bla,500,fin);  sscanf(bla,"%s",&ifo[i].noisefilepath);
    fgets(bla,500,fin);  sscanf(bla,"%s",subdir);
    sprintf(ifo[i].noisefilepath,"%s%s%s",datadir,"/",subdir);
    fgets(bla,500,fin);  sscanf(bla,"%s",ifo[i].noisefileprefix);
    fgets(bla,500,fin);  sscanf(bla,"%s",ifo[i].noisefilesuffix);
    fgets(bla,500,fin);  sscanf(bla,"%d",&ifo[i].noisefilesize);
    fgets(bla,500,fin);  sscanf(bla,"%d",&ifo[i].noisefileoffset);
    fgets(bla,500,fin);  sscanf(bla,"%d",&ifo[i].noisedoubleprecision);
    
  }
  fclose(fin);
  
}





// All parameters that are read in here should be members of the runvar struct
void readParameterInputfile(struct runpar *run)
{
  int i;
  char bla[500];
  FILE *fin;
  
  if((fin = fopen(run->parameterinfilename,"r")) == NULL) {
    printf("   Error reading parameter input file: %s, aborting.\n\n\n",run->parameterinfilename);
    exit(1);
  } else {
    printf("   Using parameter input file: %s.\n",run->parameterinfilename);
  }
  
  
  //Use and l for floats: %lf, %lg, etc, rather than %f, %g
  
  for(i=1;i<=2;i++) fgets(bla,500,fin);  //Read first 2 lines
  
  //Priors:
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->priorSet);
  
  
  
  //Parameters:
  for(i=1;i<=5;i++) fgets(bla,500,fin);  //Read empty and comment lines
  
  for(i=0;i<run->nMCMCpar;i++) {
    fscanf(fin,"%d %d %lf %d %d %lf %d %lf %lf",&run->parNumber[i],&run->parID[i],&run->parBestVal[i],&run->parFix[i],&run->parStartMCMC[i],&run->parSigma[i],&run->priorType[i],&run->priorBoundLow[i],&run->priorBoundUp[i]);
    fgets(bla,500,fin);  //Read rest of the line
    
    //printf("%d:  %d %d %lf %d %lf %d %lf %lf\n",i,run->parNumber[i],run->parID[i],run->parBestVal[i],run->parStartMCMC[i],run->parSigma[i],
    //   run->priorType[i],run->priorBoundLow[i],run->priorBoundUp[i]);
    
    
    if(run->parNumber[i] != i+1) {
      printf("   Error reading parameter input file %s:  parameter %d has number %d.\n   Aborting...\n\n",run->parameterinfilename,i+1,run->parNumber[i]);
      exit(1);
    }
    
    
    run->parRevID[run->parID[i]] = i;  //Reverse parameter ID
    
    
    
    // Get the desired boundary conditions:
    switch (run->priorType[i]) {
    case 11 :
      break;
    case 12 :
      if(run->priorBoundLow[i] > 0.0 || run->priorBoundUp[i] < 0.0) {
	printf("\n\n   Error reading parameter input file %s, parameter %d (%s):\n     for priorType = 12, priorBoundLow and priorBoundUp must be <= 0 and >= 0 respectively.\n   Aborting...\n\n",
	       run->parameterinfilename,run->parNumber[i],run->parAbrev[run->parID[i]]);
	exit(1);
      }
      run->priorBoundLow[i] = run->parBestVal[i] + run->priorBoundLow[i];
      run->priorBoundUp[i]  = run->parBestVal[i] + run->priorBoundUp[i];
      break;
    case 13 :
      if(run->priorBoundLow[i] > 1.0 || run->priorBoundUp[i] < 1.0) {
	printf("\n\n   Error reading parameter input file %s, parameter %d (%s):\n     for priorType = 13, priorBoundLow and priorBoundUp must be <= 1 and >= 1 respectively.\n   Aborting...\n\n",
	       run->parameterinfilename,run->parNumber[i],run->parAbrev[run->parID[i]]);
	exit(1);
      }
      run->priorBoundLow[i] = run->parBestVal[i] * run->priorBoundLow[i];
      run->priorBoundUp[i]  = run->parBestVal[i] * run->priorBoundUp[i];
      break;
    case 21 : 
      run->priorBoundLow[i] = 0.0;
      run->priorBoundUp[i]  = tpi;
      break;
    case 22 :
      run->priorBoundLow[i] = 0.0;
      run->priorBoundUp[i]  = pi;
      break;
    default :
      printf("\n\n   Error reading parameter input file %s, parameter %d (%s):\n     %d is not a valid option for priorType.\n   Aborting...\n\n",
	     run->parameterinfilename,run->parNumber[i],run->parAbrev[run->parID[i]],run->priorType[i]);
      exit(1);
    } //End switch
    
    
    // Check whether value for fix is valid
    if(run->parFix[i] < 0 || run->parFix[i] > 2) {
      printf("\n\n   Error reading parameter input file %s, parameter %d (%s):\n     %d is not a valid option for parFix.\n   Aborting...\n\n",
	     run->parameterinfilename,run->parNumber[i],run->parAbrev[run->parID[i]],run->parFix[i]);
	exit(1);
    }      
    
    // Check whether value for start is valid
    if(run->parStartMCMC[i] < 1 || run->parStartMCMC[i] > 5) {
      printf("\n\n   Error reading parameter input file %s, parameter %d (%s):\n     %d is not a valid option for parStartMCMC.\n   Aborting...\n\n",
	     run->parameterinfilename,run->parNumber[i],run->parAbrev[run->parID[i]],run->parStartMCMC[i]);
	exit(1);
    }      
    
    //Check whether the lower prior boundary < the upper
    if(run->priorBoundLow[i] >= run->priorBoundUp[i]) {
      printf("\n\n   Error reading parameter input file %s, parameter %d (%s):\n     the lower boundary of the prior is larger than or equal to the upper boundary (%lf vs. %lf).\n   Aborting...\n\n",
	     run->parameterinfilename,run->parNumber[i],run->parAbrev[run->parID[i]],run->priorBoundLow[i],run->priorBoundUp[i]);
      exit(1);
    }
    
    //Check whether  lower prior boundary <= best value <= upper boundary
    if(run->parBestVal[i] < run->priorBoundLow[i] || run->parBestVal[i] > run->priorBoundUp[i]) {
      printf("\n\n   Error reading parameter input file %s, parameter %d (%s):\n     the best value lies outside the prior range.\n   Aborting...\n\n",
	     run->parameterinfilename,run->parNumber[i],run->parAbrev[run->parID[i]]);
      exit(1);
    }
    
  } //End for
  
  fclose(fin);
  
  
  
  //Print MCMC parameters and prior ranges to screen:
  char FixStr[3][99];
  strcpy(FixStr[0],"No, let it free");
  strcpy(FixStr[1],"Yes, to best value");
  strcpy(FixStr[2],"Yes, to injection");
  
  char StartStr[6][99];
  strcpy(StartStr[1],"From best value");
  strcpy(StartStr[2],"Randomly near best value");
  strcpy(StartStr[3],"From injection value");
  strcpy(StartStr[4],"Randomly near injection");
  strcpy(StartStr[5],"Randomly from prior");
  
  printf("\n      Nr: Name:                Best value:     Prior:     min:            max:    Fix parameter?        Start chain:\n");
  for(i=0;i<run->nMCMCpar;i++) {
    printf("      %2d  %-11s     %15.4lf     %15.4lf %15.4lf     %-20s  %-25s\n",run->parNumber[i],run->parAbrev[run->parID[i]],run->parBestVal[i],
	   run->priorBoundLow[i],run->priorBoundUp[i],  FixStr[run->parFix[i]],StartStr[run->parStartMCMC[i]]);
  }
  printf("\n");
  
}



void setParameterNames(struct runpar * run)
{
  //Set 01: time
  strcpy(run->parAbrev[11], "t_c");
  run->parDef[11] = 1;
  strcpy(run->parAbrev[12], "t_40");
  run->parDef[12] = 1;
  
  //Set 02: distance
  strcpy(run->parAbrev[21], "d^3");
  run->parDef[21] = 1;
  strcpy(run->parAbrev[22], "log(d)");
  run->parDef[22] = 1;
  
  //Set 03: sky position
  strcpy(run->parAbrev[31], "R.A.");
  run->parDef[31] = 1;
  strcpy(run->parAbrev[32], "sin(dec)");
  run->parDef[32] = 1;
  
  //Set 04: phase
  strcpy(run->parAbrev[41], "phi_orb");
  run->parDef[41] = 1;
  
  //Set 05: orientation
  strcpy(run->parAbrev[51], "cos(i)");
  run->parDef[51] = 1;
  strcpy(run->parAbrev[52], "psi");
  run->parDef[52] = 1;
  strcpy(run->parAbrev[53], "sin th_J0");
  run->parDef[53] = 1;
  strcpy(run->parAbrev[54], "phi_J0");
  run->parDef[54] = 1;
  
  //Set 06: mass
  strcpy(run->parAbrev[61], "Mc");
  run->parDef[61] = 1;
  strcpy(run->parAbrev[62], "eta");
  run->parDef[62] = 1;
  strcpy(run->parAbrev[63], "M1");
  run->parDef[63] = 1;
  strcpy(run->parAbrev[64], "M2");
  run->parDef[64] = 1;
  
  //Set 07: spin
  strcpy(run->parAbrev[71], "a_spin1");
  run->parDef[71] = 1;
  strcpy(run->parAbrev[72], "cs th_sp1");
  run->parDef[72] = 1;
  strcpy(run->parAbrev[73], "phi_spin1");
  run->parDef[73] = 1;
  strcpy(run->parAbrev[74], "a_spin2");
  run->parDef[74] = 1;
  strcpy(run->parAbrev[75], "cs th_sp2");
  run->parDef[75] = 1;
  strcpy(run->parAbrev[76], "phi_spin2");
  run->parDef[76] = 1;
  
  //Set 08: merger, ringdown
  //strcpy(run->parAbrev[], "");
  //run->parDef[] = 1;
  
  //Set:
  //strcpy(run->parAbrev[], "");
  //run->parDef[] = 1;
  
}






// Set the global variables.
// Many of these are now in the input file or unused.
// This routine should eventually contain mathematical and (astro)physical constants only
void setconstants()
{
  tempi = 0; //A global variable that determines the current chain (temperature) in the temperature ladder
  
  // Mathematical constants:
  pi   = 3.141592653589793;   // pi
  tpi  = 6.283185307179586;   // 2 pi
  mtpi = 6.283185307179586e6; // Large multiple of 2 pi (2 megapi)
  
  // Define some physical constants:
  G    = 6.67259e-11;         // 6.674215e-11; */ /* gravity constant (SI)
  c    = 299792458.0;         // speed of light (m/s)
  
  Ms   = 1.9889194662e30;     // solar mass (kg)
  Mpc  = 3.08568025e22;       // metres in a Mpc  (LAL: 3.0856775807e22)
  Mpcs = 1.029272137e14;      // seconds in a Mpc  (Mpc/c)
}








void gettrueparameters(struct parset *par, int nTruePar)  //Set the parameters for the 12-parameter spinning template to the 'true values'
{
  int i=0;
  for(i=0;i<nTruePar;i++) {
    par->par[i]      = truepar[i];
  }
  
  par->tc       = truepar[2];                    // coalescence time
  //par->longi    = fmod(longitude(truepar[6],GMST(par->tc))+mtpi,tpi);  //The parameter in the input and output is RA; the MCMC parameter is 'longi' ~ Greenwich hour angle
  par->longi    = truepar[6];                    //The parameter in the input and output is RA; the MCMC parameter is 'longi' ~ Greenwich hour angle
  par->sinlati  = truepar[7];           // sin latitude (sin(delta))  (40)     
  par->sinthJ0  = truepar[9];           // sin Theta_J0 ~ latitude, pi/2 = NP    (15)
  //par->phiJ0    = fmod(longitude(truepar[10],GMST(par->tc))+mtpi,tpi);  //The parameter in the input and output is a 'RA'; the MCMC parameter is a 'longi' ~ Greenwich hour angle
  par->phiJ0    = truepar[10];               // Phi_J0 ~ azimuthal            (125)
  
  /*
  parSTOPGREPFROMFINDINGTHIS->m1       = truepar[0];                    // M1 (10.0)
  parSTOPGREPFROMFINDINGTHIS->m2       = truepar[1];                    // M2  (1.4)
  parSTOPGREPFROMFINDINGTHIS->m        = parSTOPGREPFROMFINDINGTHIS->m1+parSTOPGREPFROMFINDINGTHIS->m2;
  parSTOPGREPFROMFINDINGTHIS->mu       = parSTOPGREPFROMFINDINGTHIS->m1*parSTOPGREPFROMFINDINGTHIS->m2/parSTOPGREPFROMFINDINGTHIS->m;
  
  parSTOPGREPFROMFINDINGTHIS->eta      = parSTOPGREPFROMFINDINGTHIS->mu/parSTOPGREPFROMFINDINGTHIS->m;                // mass ratio                
  parSTOPGREPFROMFINDINGTHIS->mc       = parSTOPGREPFROMFINDINGTHIS->m*pow(parSTOPGREPFROMFINDINGTHIS->eta,0.6);      // chirp mass. in Mo         
  parSTOPGREPFROMFINDINGTHIS->logdl    = log(truepar[3]);               // log-distance (Mpc) (17.5)             
  
  parSTOPGREPFROMFINDINGTHIS->spin     = truepar[4];                    // magnitude of total spin   (0.1)
  parSTOPGREPFROMFINDINGTHIS->kappa    = cos(truepar[5]*d2r);           // L^.S^, cos of angle between L^ & S^  (0.819152)
  
  parSTOPGREPFROMFINDINGTHIS->phase    = truepar[8]*d2r;                // orbital phase   (phi_c)   (0.2)
  parSTOPGREPFROMFINDINGTHIS->alpha    = truepar[11]*d2r;               // Alpha_c                       (0.9 rad = 51.566202deg)
  */
  
  par->loctc    = NULL;
  par->localti  = NULL;
  par->locazi   = NULL;
  par->locpolar = NULL;
}


void getstartparameters(struct parset *par, struct runpar run)  //Set the parameters for the 12-parameter spinning template to the starting values for the MCMC chain
{
  
  int i=0;
  for(i=0;i<run.nMCMCpar;i++) {
    par->par[i]      = run.startpar[i];
  }
  
  /*
  parSTOPGREPFROMFINDINGTHIS->m1       = run.startpar[0];                    // M1 (10.0)
  parSTOPGREPFROMFINDINGTHIS->m2       = run.startpar[1];                    // M2  (1.4)
  parSTOPGREPFROMFINDINGTHIS->m        = parSTOPGREPFROMFINDINGTHIS->m1+parSTOPGREPFROMFINDINGTHIS->m2;
  parSTOPGREPFROMFINDINGTHIS->mu       = parSTOPGREPFROMFINDINGTHIS->m1*parSTOPGREPFROMFINDINGTHIS->m2/parSTOPGREPFROMFINDINGTHIS->m;
  parSTOPGREPFROMFINDINGTHIS->eta      = parSTOPGREPFROMFINDINGTHIS->mu/parSTOPGREPFROMFINDINGTHIS->m;                // mass ratio                
  parSTOPGREPFROMFINDINGTHIS->mc       = parSTOPGREPFROMFINDINGTHIS->m*pow(parSTOPGREPFROMFINDINGTHIS->eta,0.6);      // chirp mass. in Mo         
  parSTOPGREPFROMFINDINGTHIS->tc       = run.startpar[2];                    // coalescence time
  parSTOPGREPFROMFINDINGTHIS->logdl    = log(run.startpar[3]);               // log-distance (Mpc) (17.5)             
  
  parSTOPGREPFROMFINDINGTHIS->spin     = run.startpar[4];                    // magnitude of total spin   (0.1)
  parSTOPGREPFROMFINDINGTHIS->kappa    = cos(run.startpar[5]*d2r);           // L^.S^, cos of angle between L^ & S^  (0.819152)
  parSTOPGREPFROMFINDINGTHIS->longi    = fmod(longitude(run.startpar[6]*h2r,GMST(parSTOPGREPFROMFINDINGTHIS->tc))+mtpi,tpi);  //The parameter in the input and output is RA; the MCMC parameter is 'longi' ~ Greenwich hour angle
  parSTOPGREPFROMFINDINGTHIS->sinlati  = sin(run.startpar[7]*d2r);           // sin latitude (sin(delta))  (40)     
  
  parSTOPGREPFROMFINDINGTHIS->phase    = run.startpar[8]*d2r;                // orbital phase   (phi_c)   (0.2)
  parSTOPGREPFROMFINDINGTHIS->sinthJ0  = sin(run.startpar[9]*d2r);           // sin Theta_J0 ~ latitude, pi/2 = NP    (15)
  parSTOPGREPFROMFINDINGTHIS->phiJ0    = run.startpar[10]*d2r;               // Phi_J0 ~ azimuthal            (125)
  parSTOPGREPFROMFINDINGTHIS->alpha    = run.startpar[11]*d2r;               // Alpha_c                       (0.9 rad = 51.566202deg)
  */
  
  par->loctc    = NULL;
  par->localti  = NULL;
  par->locazi   = NULL;
  par->locpolar = NULL;
}


//Allocate memory for the vectors in the struct parset
void allocparset(struct parset *par, int networksize)
{
  par->loctc    = NULL;
  par->localti  = NULL;
  par->locazi   = NULL;
  par->locpolar = NULL;
  
  par->loctc    = (double*)calloc(networksize,sizeof(double));
  par->localti  = (double*)calloc(networksize,sizeof(double));
  par->locazi   = (double*)calloc(networksize,sizeof(double));
  par->locpolar = (double*)calloc(networksize,sizeof(double));
}


//Deallocate the vectors in the struct parset
void freeparset(struct parset *par)
{
  free(par->loctc);         par->loctc        = NULL;
  free(par->localti);       par->localti      = NULL;
  free(par->locazi);        par->locazi       = NULL;
  free(par->locpolar);      par->locpolar     = NULL;
}







