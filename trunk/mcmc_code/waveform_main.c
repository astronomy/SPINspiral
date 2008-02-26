// waveform_main.c:
// Main routine to generate and save waveform



#include <mcmc.h>


/* main program: */
main()
{
  // Interferometers are managed via the `database'; the `network' is a vector of pointers to the database (see below).
  // The interferometers that are actually used need to be initialised via the `ifoinit()'-function in order to determine noise PSD, signal FT &c.
  
  printf("\n\n   Starting MCMC code...\n");
  int i;
  double snr;
  
  //Initialise stuff for the run
  setconstants();    //Set the global constants (which are variable in C)
  readinputfile();   //Read data for this run from input.mcmc
  setmcmcseed();     //Set mcmcseed (if 0), or use the current value
  writeinputfile();  //Write run data to nicely formatted input.mcmc.<mcmcseed>
  
  dosnr = 1;
  domcmc = 0;
  domatch = 0;
  writesignal = 0;
  
  /*--- SET UP INTERFEROMETER `DATABASE': ---*/
  struct interferometer database[3];
  set_ifo_data(database);  //Set all the data of the IFOs you may want to use (H1,L1 + VIRGO by default)
  
  //--- DEFINE INTERFEROMETER NETWORK: --- (0: H1, 1: L1, 2: Virgo)
  const int networksize = 1;
  struct interferometer *network[1] = {&database[0]};
  //const int networksize = 2;
  //struct interferometer *network[2] = {&database[0], &database[1]};
  //const int networksize = 3;
  //struct interferometer *network[3] = {&database[0], &database[1], &database[2]};
  
  
  //Initialise interferometers, read and prepare data (takes some time)
  printf("   Initialising IFOs, reading datafiles...\n");
  ifoinit(network, networksize);
  
  
  //Calculate 'null-likelihood'
  struct parset nullpar;
  setnullparameters(&nullpar);
  nullpar.loctc    = (double*)calloc(networksize,sizeof(double));
  nullpar.localti  = (double*)calloc(networksize,sizeof(double));
  nullpar.locazi   = (double*)calloc(networksize,sizeof(double));
  nullpar.locpolar = (double*)calloc(networksize,sizeof(double));
  localpar(&nullpar, network, networksize);
  NullLikelihood = net_loglikelihood(&nullpar, networksize, network);
  
  //Get a parameter set to calculate SNR or write the wavefrom to disc
  struct parset dummypar;
  settrueparameters(&dummypar);
  dummypar.loctc    = (double*)calloc(networksize,sizeof(double));
  dummypar.localti  = (double*)calloc(networksize,sizeof(double));
  dummypar.locazi   = (double*)calloc(networksize,sizeof(double));
  dummypar.locpolar = (double*)calloc(networksize,sizeof(double));
  localpar(&dummypar, network, networksize);
  
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
  if(domcmc==1) {
    //printmuch=1;
    mcmc(networksize, network);
    //printmuch=0;
  } else {
  */
    printf("%10s  %10s  %6s  %20s  %6s  ","niter","nburn","seed","null likelihood","ndet");
    for(i=0;i<networksize;i++) {
      printf("%16s%4s  ",network[i]->name,"SNR");
    }
    printf("\n%10d  %10d  %6d  %20.10lf  %6d  ",iter,nburn,mcmcseed,NullLikelihood,networksize);
    for(i=0;i<networksize;i++) {
      printf("%20.10lf  ",network[i]->snr);
    }
    printf("\n\n%8s  %8s  %17s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n",
	   "Mc","eta","tc","logdL","sinlati","longi","phase","spin","kappa","sinthJ0","phiJ0","alpha");
    printf("%8.5f  %8.5f  %17.6lf  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n\n",
	   dummypar.mc,dummypar.eta,dummypar.tc,dummypar.logdl,dummypar.sinlati,dummypar.longi,dummypar.phase,dummypar.spin,dummypar.kappa,dummypar.sinthJ0,dummypar.phiJ0,dummypar.alpha);
    //}
  
  //  struct parset dummypar;
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








void pardispose(struct parset *par)
{
  free(par->loctc);         par->loctc        = NULL;
  free(par->localti);       par->localti      = NULL;
  free(par->locazi);        par->locazi       = NULL;
  free(par->locpolar);      par->locpolar     = NULL;
}




void setmcmcseed()
//If mcmcseed==0, set it using the system clock
{
  struct timeval time;
  struct timezone tz;
  gettimeofday(&time, &tz);
  if(mcmcseed==0) mcmcseed = time.tv_usec;
}







