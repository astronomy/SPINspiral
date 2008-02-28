// waveform_main.c:
// Main routine to generate and save waveform



#include <mcmc.h>


// Main program:
int main(int argc, char * argv[])
{
  // Interferometers are managed via the `database'; the `network' is a vector of pointers to the database (see below).
  // The interferometers that are actually used need to be initialised via the `ifoinit()'-function in order to determine noise PSD, signal FT &c.
  
  printf("\n\n   Starting MCMC code...\n");
  int i;
  double snr;
  
  //Initialise stuff for the run
  struct runpar run;
  sprintf(run.infilename,"mcmc.input");
  if(argc > 1) sprintf(run.infilename,argv[1]);

  setconstants(&run);    //Set the global constants (which are variable in C). This routine should eventuelly disappear.
  readlocalfile();       //Read system-dependent data, e.g. path to data files
  readinputfile(&run);   //Read data for this run from input.mcmc
  setmcmcseed(&run);     //Set mcmcseed (if 0), or use the current value
  writeinputfile(&run);  //Write run data to nicely formatted input.mcmc.<mcmcseed>
  
  dosnr = 1;
  domcmc = 0;
  domatch = 0;
  writesignal = 0;
  
  //Set up the data of the IFOs you may want to use (H1,L1 + VIRGO by default)
  struct interferometer database[3];
  set_ifo_data(database);
  
  
  //Define interferometer network; how many and which IFOs
  //const int networksize = 1;
  //struct interferometer *network[1] = {&database[0]};
  const int networksize = 2;
  struct interferometer *network[2] = {&database[0], &database[1]};
  //const int networksize = 3;
  //struct interferometer *network[3] = {&database[0], &database[1], &database[2]};
  
  
  //Initialise interferometers, read and prepare data, inject signal (takes some time)
  if(networksize == 1) {
    printf("   Initialising 1 IFO, reading data...\n");
  } else {
    printf("   Initialising %d IFOs, reading datafiles...\n",networksize);
  }
  ifoinit(network, networksize);
  if(inject) {
    printf("   A signal with the 'true' parameter values was injected.\n");
  } else {
    printf("   No signal was injected.\n");
  }
  
  
  //Calculate 'null-likelihood'
  struct parset nullpar;
  setnullparameters(&nullpar);
  nullpar.loctc    = (double*)calloc(networksize,sizeof(double));
  nullpar.localti  = (double*)calloc(networksize,sizeof(double));
  nullpar.locazi   = (double*)calloc(networksize,sizeof(double));
  nullpar.locpolar = (double*)calloc(networksize,sizeof(double));
  localpar(&nullpar, network, networksize);
  NullLikelihood = net_loglikelihood(&nullpar, networksize, network);
  if(inject == 0) NullLikelihood *= 1.01;  //If no signal is injected, presumably there is one present in the data; enlarge the range that log(L) can take by owering Lo (since L>Lo is forced)
  
  //Get a parameter set to calculate SNR or write the wavefrom to disc
  struct parset dummypar;
  settrueparameters(&dummypar);
  dummypar.loctc    = (double*)calloc(networksize,sizeof(double));
  dummypar.localti  = (double*)calloc(networksize,sizeof(double));
  dummypar.locazi   = (double*)calloc(networksize,sizeof(double));
  dummypar.locpolar = (double*)calloc(networksize,sizeof(double));
  localpar(&dummypar, network, networksize);
  
  
  //Calculate SNR
  if(dosnr==1) {
    for (i=0; i<networksize; ++i) {
      //printmuch=1;
      snr = signaltonoiseratio(&dummypar, network, i);
      //printmuch=0;
      //printf("  %10.2f\n",snr);
      network[i]->snr = snr;
      //printf("  Det: %d,  GMST: %20.10lf,  Azimuth: %20.10lf,  Altitude: %20.10lf  \n",i,GMST(prior_tc_mean),dummypar.locazi[i],dummypar.localti[i]);
    }
  }
  
  //Write the signal and its FFT to disc
  writesignal = 1;
  if(writesignal) {
    for (i=0; i<networksize; ++i) {
      //printmuch=1;
      writesignaltodisc(&dummypar, network, i);
      //printmuch=0;
    }
  }
  writesignal=0;
  
  /*
  //Do MCMC
  if(domcmc==1) {
    //printmuch=1;
    mcmc(&run, networksize, network);
    //printmuch=0;
  } else {
  */
  printf("%10s  %10s  %6s  %20s  %6s  ","niter","nburn","seed","null likelihood","ndet");
  for(i=0;i<networksize;i++) {
    printf("%16s%4s  ",network[i]->name,"SNR");
  }
  printf("\n%10d  %10d  %6d  %20.10lf  %6d  ",iter,nburn,run.mcmcseed,NullLikelihood,networksize);
  for(i=0;i<networksize;i++) {
    printf("%20.10lf  ",network[i]->snr);
  }
  printf("\n\n%8s  %8s  %17s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n",
	 "Mc","eta","tc","logdL","sinlati","longi","phase","spin","kappa","sinthJ0","phiJ0","alpha");
  printf("%8.5f  %8.5f  %17.6lf  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n\n",
	 dummypar.mc,dummypar.eta,dummypar.tc,dummypar.logdl,dummypar.sinlati,dummypar.longi,dummypar.phase,dummypar.spin,dummypar.kappa,dummypar.sinthJ0,dummypar.phiJ0,dummypar.alpha);
  //}
  
  //Calculate matches between two signals
  if(domatch==1) {
    printf("\n");
    settrueparameters(&dummypar);
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
  for (i=0; i<networksize; ++i)
    ifodispose(network[i]);
  
  
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







