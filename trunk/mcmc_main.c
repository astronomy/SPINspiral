// mcmc_main.c
// Main routine of spinning MCMC code



#include <mcmc.h>


// Main program:
int main(int argc, char * argv[])
{
  // Interferometers are managed via the `database'; the `network' is a vector of pointers to the database (see below).
  // The interferometers that are actually used need to be initialised via the `ifoinit()'-function in order to determine noise PSD, signal FT &c.
  
  printf("\n\n   Starting MCMC code...\n");
  int i;
  double snr;
  
  useoldmcmcoutputformat = 0; //Set to 1 if you want to ... exactly!
  
  //Initialise stuff for the run
  struct runpar run;
  setconstants(&run);    //Set the global constants (which are variable in C). This routine should eventually disappear.
  run.setranpar  = (int*)calloc(npar,sizeof(int));
  sprintf(run.infilename,"mcmc.input"); //Default input filename
  if(argc > 1) sprintf(run.infilename,argv[1]);
  
  readlocalfile();                //Read system-dependent data, e.g. path to data files
  readinputfile(&run);            //Read data for this run from input.mcmc
  //setmcmcseed(&run);              //Set mcmcseed if 0, otherwise keep the current value
  setseed(&run.mcmcseed);         //Set mcmcseed if 0, otherwise keep the current value; more general routine
  setrandomtrueparameters(&run);  //Randomise the injection parameters where wanted
  writeinputfile(&run);           //Write run data to nicely formatted input.mcmc.<mcmcseed>
  
  //Set up the data for the IFOs you may want to use (H1,L1 + VIRGO by default)
  struct interferometer database[3];
  set_ifo_data(run, database);  
  
  //Define interferometer network which IFOs.  The first run.networksize are actually used
  struct interferometer *network[3] = {&database[0], &database[1], &database[2]};
  int networksize = run.networksize;
  
  //Initialise interferometers, read and prepare data, inject signal (takes some time)
  if(networksize == 1) {
    printf("   Initialising 1 IFO, reading noise and data...\n");
  } else {
    printf("   Initialising %d IFOs, reading noise and data files...\n",networksize);
  }
  ifoinit(network, networksize);
  if(inject) {
    if(run.targetsnr < 0.001) printf("   A signal with the 'true' parameter values was injected.\n");
  } else {
    printf("   No signal was injected.\n");
  }
  
  
  //Get a parameter set to calculate SNR or write the wavefrom to disc
  struct parset dummypar;
  gettrueparameters(&dummypar);
  dummypar.loctc    = (double*)calloc(networksize,sizeof(double));
  dummypar.localti  = (double*)calloc(networksize,sizeof(double));
  dummypar.locazi   = (double*)calloc(networksize,sizeof(double));
  dummypar.locpolar = (double*)calloc(networksize,sizeof(double));
  localpar(&dummypar, network, networksize);
  
  
  //Calculate SNR
  run.netsnr = 0.0;
  if(dosnr==1) {
    for (i=0; i<networksize; ++i) {
      snr = signaltonoiseratio(&dummypar, network, i);
      network[i]->snr = snr;
      run.netsnr += snr*snr;
    }
    run.netsnr = sqrt(run.netsnr);
  }
  
  //Get the desired SNR by scaling the distance
  if(run.targetsnr > 0.001 && inject==1) {
    truepar[3] *= (run.netsnr/run.targetsnr);  //Use total network SNR
    //truepar[3] *= (run.netsnr/(run.targetsnr*sqrt((double)networksize)));  //Use geometric average SNR
    printf("   Setting distance to %lf Mpc to get a network SNR of %lf.\n",truepar[3],run.targetsnr);
    gettrueparameters(&dummypar);
    dummypar.loctc    = (double*)calloc(networksize,sizeof(double));
    dummypar.localti  = (double*)calloc(networksize,sizeof(double));
    dummypar.locazi   = (double*)calloc(networksize,sizeof(double));
    dummypar.locpolar = (double*)calloc(networksize,sizeof(double));
    localpar(&dummypar, network, networksize);
    
    //Recalculate SNR
    run.netsnr = 0.0;
    if(dosnr==1) {
      for (i=0; i<networksize; ++i) {
	snr = signaltonoiseratio(&dummypar, network, i);
	network[i]->snr = snr;
	run.netsnr += snr*snr;
      }
      run.netsnr = sqrt(run.netsnr);
    }
    
    for (i=0; i<networksize; ++i)
      ifodispose(network[i]);
    //Reinitialise interferometers, read and prepare data, inject signal (takes some time)
    if(networksize == 1) {
      printf("   Reinitialising 1 IFO, reading data...\n");
    } else {
      printf("   Reinitialising %d IFOs, reading datafiles...\n",networksize);
    }
    ifoinit(network, networksize);
    printf("   A signal with the 'true' parameter values was injected.\n");
  }
  
  
  
  //Calculate 'null-likelihood'
  struct parset nullpar;
  getnullparameters(&nullpar);
  nullpar.loctc    = (double*)calloc(networksize,sizeof(double));
  nullpar.localti  = (double*)calloc(networksize,sizeof(double));
  nullpar.locazi   = (double*)calloc(networksize,sizeof(double));
  nullpar.locpolar = (double*)calloc(networksize,sizeof(double));
  localpar(&nullpar, network, networksize);
  run.logL0 = net_loglikelihood(&nullpar, networksize, network);
  if(inject == 0) run.logL0 *= 1.01;  //If no signal is injected, presumably there is one present in the data; enlarge the range that log(L) can take by owering Lo (since L>Lo is forced)
  
  
  
  //Write the signal and its FFT to disc
  if(writesignal) {
    for (i=0; i<networksize; ++i) {
      //printmuch=1;
      writesignaltodisc(&dummypar, network, i);
      //printmuch=0;
    }
  }
  writesignal=0;
  
  
  //Write some injection-signal parameters to screen (used to be in the MCMC routine)
  printf("\n\n");
  printf("   Global     :    Source position:  RA: %5.2lfh, dec: %6.2lfd;  J0 points to:  RA: %5.2lfh, dec: %6.2lfd;   inclination J0: %5.2lfd  \n",rightAscension(dummypar.longi,GMST(dummypar.tc))*r2h,asin(dummypar.sinlati)*r2d,rightAscension(dummypar.phiJ0,GMST(dummypar.tc))*r2h,asin(dummypar.sinthJ0)*r2d,(pi/2.0-acos(dummypar.NdJ))*r2d);
  for(i=0;i<networksize;i++) printf("   %-11s:    theta: %5.1lfd,  phi: %5.1lfd;   azimuth: %5.1lfd,  altitude: %5.1lfd\n",network[i]->name,dummypar.localti[i]*r2d,dummypar.locazi[i]*r2d,fmod(pi-(dummypar.locazi[i]+network[i]->rightarm)+mtpi,tpi)*r2d,(pi/2.0-dummypar.localti[i])*r2d);
  
  printf("\n  %10s  %10s  %6s  %20s  %6s  ","niter","nburn","seed","null likelihood","ndet");
  for(i=0;i<networksize;i++) printf("%16s%4s  ",network[i]->name,"SNR");
  printf("%20s  ","Network SNR");
  printf("\n  %10d  %10d  %6d  %20.10lf  %6d  ",iter,nburn,run.mcmcseed,run.logL0,networksize);
  for(i=0;i<networksize;i++) printf("%20.10lf  ",network[i]->snr);
  printf("%20.10lf\n\n",run.netsnr);
  printf("    %8s  %8s  %17s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n", "Mc","eta","tc","logdL","spin","kappa","longi","sinlati","phase","sinthJ0","phiJ0","alpha");
  printf("    %8.5f  %8.5f  %17.6lf  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n\n", dummypar.mc,dummypar.eta,dummypar.tc,dummypar.logdl,dummypar.spin,dummypar.kappa,dummypar.longi,dummypar.sinlati,dummypar.phase,dummypar.sinthJ0,dummypar.phiJ0,dummypar.alpha);
  printf("    %8s  %8s  %17s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n", "M1","M2","tc","d_L","spin","th_SL","RA","Dec","phase","th_J0","phi_J0","alpha");
  printf("    %8.5f  %8.5f  %17.6lf  %8.2f  %8.5f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n\n", dummypar.m1,dummypar.m2,dummypar.tc,exp(dummypar.logdl),dummypar.spin,acos(dummypar.kappa)*r2d,rightAscension(dummypar.longi,GMST(dummypar.tc))*r2h,asin(dummypar.sinlati)*r2d,dummypar.phase*r2d,asin(dummypar.sinthJ0)*r2d,dummypar.phiJ0*r2d,dummypar.alpha*r2d);
  
  
  
  //Do MCMC
  if(domcmc==1) {
    //printmuch=1;
    mcmc(&run, network);
    //printmuch=0;
  }
  
  
  
  //Calculate matches between two signals
  if(domatch==1) {
    printf("\n");
    gettrueparameters(&dummypar);
    dummypar.loctc    = (double*)calloc(networksize,sizeof(double));
    dummypar.localti  = (double*)calloc(networksize,sizeof(double));
    dummypar.locazi   = (double*)calloc(networksize,sizeof(double));
    dummypar.locpolar = (double*)calloc(networksize,sizeof(double));
    
    FILE *fout;
    fout = fopen("tc.dat","w");
    double fac=0.0;
    double matchpar = dummypar.tc,matchres=0.0;
    for(fac=-0.002;fac<0.002;fac+=0.00005) {
      dummypar.tc = matchpar+fac;
      for(i=0;i<networksize;i++) {
        localpar(&dummypar, network, networksize);
        matchres = match(&dummypar,network,i,networksize);
        printf("%10.6f  %10.6f\n",fac,matchres);
        fprintf(fout,"%10.6f  %10.6f\n",fac,matchres);
      }
    }
    fclose(fout);
  }
  
  
  //Get rid of allocated memory and quit
  for (i=0; i<networksize; ++i) ifodispose(network[i]);
  
  free(run.setranpar);
  
  free(nullpar.loctc);
  free(nullpar.localti);
  free(nullpar.locazi);
  free(nullpar.locpolar);
  
  free(dummypar.loctc);
  free(dummypar.localti);
  free(dummypar.locazi);
  free(dummypar.locpolar);
  
  printf("\n   MCMC code done.\n\n\n");
  return 0;
}







//Deallocate the struct parset
void pardispose(struct parset *par)
{
  free(par->loctc);         par->loctc        = NULL;
  free(par->localti);       par->localti      = NULL;
  free(par->locazi);        par->locazi       = NULL;
  free(par->locpolar);      par->locpolar     = NULL;
}


void setmcmcseed(struct runpar *run)
//If run->mcmcseed==0, set it using the system clock
{
  struct timeval time;
  struct timezone tz;
  gettimeofday(&time, &tz);
  if(run->mcmcseed==0) run->mcmcseed = time.tv_usec;
}


void setseed(int *seed)
//If seed==0, set (randomise) it using the system clock
{
  struct timeval time;
  struct timezone tz;
  gettimeofday(&time, &tz);
  if(*seed==0) *seed = time.tv_usec;
}







