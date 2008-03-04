#include <mcmc.h>


// MCMC routine
//****************************************************************************************************************************************************  
void mcmc(struct runpar *run, int networksize, struct interferometer *ifo[])
//****************************************************************************************************************************************************  
{
  if(MvdSdebug) printf("MCMC\n");
  struct parset state;                        // MCMC/template parameter set
  
  struct mcmcvariables mcmc;                  // MCMC variables
  mcmc.ntemps = run->ntemps;                  // From the input file; copy from run to mcmc struct
  mcmc.npar = npar;                           // Copy from the global variable
  if(partemp==0) mcmc.ntemps=1;
  mcmc.temp = max(temp0,1.0);
  mcmc.logL0 = run->logL0;                    // Copy the log 'null-likelihood' from the ru struct to the mcmc struct
  mcmc.seed = run->mcmcseed;                  // Copy mcmcseed from run to mcmc struct
  mcmc.ran = gsl_rng_alloc(gsl_rng_mt19937);  // GSL random-number seed
  gsl_rng_set(mcmc.ran, mcmc.seed);           // Set seed for this run
  
  
  
  char outfilename[99];
  
  if(nburn0>=nburn) {
    printf("\n\n   *** Warning: nburn0 > nburn, setting nburn0 = nburn*0.9 ***\n\n");
    nburn0 = (int)(0.9*(double)nburn);
  }
  
  
  
  // ********************************************************************************************************************************************************************************
  // ***  MEMORY ALLOCATION  ********************************************************************************************************************************************************
  // ********************************************************************************************************************************************************************************
  
  int i=0,i1=0,j=0,j1=0,j2=0,iteri=0,tempj=0;
  double tmpdbl=0.0,tempratio=1.0;
  int nstart=0;
  double mu[npar],stdev[npar];
  
  int *corrupdate,*acceptelems;
  corrupdate = (int*)calloc(mcmc.ntemps,sizeof(int));
  acceptelems = (int*)calloc(mcmc.ntemps,sizeof(int));
  for(i=0;i<mcmc.ntemps;i++) {
    corrupdate[i] = 0;
    acceptelems[i] = 0;
  }
  
  double *newtemps,*tempampl,*maxdlogL,*sumdlogL,*avgdlogL,*expdlogL;
  mcmc.temps = (double*)calloc(mcmc.ntemps,sizeof(double));
  newtemps = (double*)calloc(mcmc.ntemps,sizeof(double));
  tempampl = (double*)calloc(mcmc.ntemps,sizeof(double));
  mcmc.logL = (double*)calloc(mcmc.ntemps,sizeof(double));
  mcmc.nlogL = (double*)calloc(mcmc.ntemps,sizeof(double));
  mcmc.dlogL = (double*)calloc(mcmc.ntemps,sizeof(double));
  maxdlogL = (double*)calloc(mcmc.ntemps,sizeof(double));
  sumdlogL = (double*)calloc(mcmc.ntemps,sizeof(double));
  avgdlogL = (double*)calloc(mcmc.ntemps,sizeof(double));
  expdlogL = (double*)calloc(mcmc.ntemps,sizeof(double));
  for(i=0;i<mcmc.ntemps;i++) {
    mcmc.temps[i] = 0.0;
    newtemps[i] = 0.0;
    tempampl[i] = 0.0;
    mcmc.logL[i] = 0.0;
    mcmc.nlogL[i] = 0.0;
    mcmc.dlogL[i] = 0.0;
    maxdlogL[i] = -1.e30;
    sumdlogL[i] = 0.0;
    avgdlogL[i] = 0.0;
    expdlogL[i] = 0.0;
  }
  
  int *ihist,*swapTs1,*swapTs2;
  mcmc.corrsig = (double*)calloc(mcmc.ntemps,sizeof(double));
  swapTs1 = (int*)calloc(mcmc.ntemps,sizeof(int));
  swapTs2 = (int*)calloc(mcmc.ntemps,sizeof(int));
  mcmc.acceptprior = (int*)calloc(mcmc.ntemps,sizeof(int));
  ihist = (int*)calloc(mcmc.ntemps,sizeof(int));
  for(i=0;i<mcmc.ntemps;i++) {
    mcmc.corrsig[i] = 1.0;
    swapTs1[i] = 0;   //Number of swaps between parallel chains
    swapTs2[i] = 0;   //Number of swaps between parallel chains
    mcmc.acceptprior[i] = 1;
    ihist[i] = 0;
  }
  
  mcmc.accepted = (int**)calloc(mcmc.ntemps,sizeof(int*));         // Count accepted proposals
  mcmc.swapTss = (int**)calloc(mcmc.ntemps,sizeof(int*));          // Count swaps between chains
  mcmc.param = (double**)calloc(mcmc.ntemps,sizeof(double*));      // The old parameters for all chains
  mcmc.nparam = (double**)calloc(mcmc.ntemps,sizeof(double*));     // The new parameters for all chains
  mcmc.maxparam = (double**)calloc(mcmc.ntemps,sizeof(double*));   // The best parameters for all chains (max logL)
  mcmc.sig = (double**)calloc(mcmc.ntemps,sizeof(double*));        // The standard deviation of the gaussian to draw the jump size from
  mcmc.scale = (double**)calloc(mcmc.ntemps,sizeof(double*));      // The rate of adaptation
  for(i=0;i<mcmc.ntemps;i++) {
    mcmc.accepted[i] = (int*)calloc(npar,sizeof(int));
    mcmc.swapTss[i] = (int*)calloc(mcmc.ntemps,sizeof(int));
    mcmc.param[i] = (double*)calloc(npar,sizeof(double));
    mcmc.nparam[i] = (double*)calloc(npar,sizeof(double));
    mcmc.maxparam[i] = (double*)calloc(npar,sizeof(double));
    mcmc.sig[i] = (double*)calloc(npar,sizeof(double));
    mcmc.scale[i] = (double*)calloc(npar,sizeof(double));
  }
  
  double **tempcovar;
  tempcovar = (double**)calloc(npar,sizeof(double*)); // A temp Cholesky-decomposed matrix
  for(i=0;i<npar;i++) {
    tempcovar[i] = (double*)calloc(npar,sizeof(double));
  }
  
  double ***covar1;
  mcmc.hist    = (double***)calloc(mcmc.ntemps,sizeof(double**)); // Store past iterations to calculate the covariances
  mcmc.covar  = (double***)calloc(mcmc.ntemps,sizeof(double**)); // The Cholesky-decomposed matrix
  covar1  = (double***)calloc(mcmc.ntemps,sizeof(double**)); // The actual covariance matrix
  for(i=0;i<mcmc.ntemps;i++) {
    mcmc.hist[i]    = (double**)calloc(npar,sizeof(double*));
    mcmc.covar[i]  = (double**)calloc(npar,sizeof(double*));
    covar1[i]  = (double**)calloc(npar,sizeof(double*));
    for(j=0;j<npar;j++) {
      mcmc.hist[i][j]    = (double*)calloc(ncorr,sizeof(double));
      mcmc.covar[i][j]  = (double*)calloc(npar,sizeof(double));
      covar1[i][j]  = (double*)calloc(npar,sizeof(double));
    }
  }
  
  
  
  
  
  // ********************************************************************************************************************************************************************************
  // ***  INITIALISE PARALLEL TEMPERING  ********************************************************************************************************************************************
  // ********************************************************************************************************************************************************************************
  
  // *** Set up temperature ladder ***
  if(mcmc.ntemps>1) {  
    tempratio = exp(log(tempmax)/(double)(mcmc.ntemps-1));
    if(prpartempinfo>0) {
      printf("   Temperature ladder:\n     Number of chains:%3d,  Tmax:%7.2lf, Ti/Ti-1:%7.3lf\n",mcmc.ntemps,tempmax,tempratio);
      if(partemp==1) printf("     Using fixed temperatures for the chains\n");
      if(partemp==2) printf("     Using sinusoid temperatures for the chains\n");
      if(partemp==3) printf("     Using a manual temperature ladder with fixed temperatures for the chains\n");
      if(partemp==4) printf("     Using a manual temperature ladder with sinusoid temperatures for the chains\n");
      printf("     Chain     To     Ampl.    Tmin     Tmax\n");
    }
    for(tempi=0;tempi<mcmc.ntemps;tempi++) {
      mcmc.temps[tempi] = pow(10.0,log10(tempmax)/(double)(mcmc.ntemps-1)*(double)tempi);
      if(partemp==3 || partemp==4) mcmc.temps[tempi] = run->temps[tempi];  //Set manual ladder
      
      //if(tempi>0) tempampl[tempi] = (mcmc.temps[tempi] - mcmc.temps[tempi-1])/(tempratio+1.0)*tempratio;  //Temperatures of adjacent chains just touch at extrema (since in antiphase)
      //if(tempi>0) tempampl[tempi] = 1.5*(mcmc.temps[tempi] - mcmc.temps[tempi-1])/(tempratio+1.0)*tempratio;  //Temperatures of adjacent chains overlap somewhat at extrema (since in antiphase)
      if(tempi>0) tempampl[tempi] = min(3.0*(mcmc.temps[tempi] - mcmc.temps[tempi-1])/(tempratio+1.0)*tempratio , fabs(mcmc.temps[tempi]-mcmc.temps[tempi-1]));  //Temperatures of adjacent chains overlap a lot at extrema (since in antiphase), make sure Ti,min>=T-i1,0
      if(mcmc.ntemps>10 && tempi>1) tempampl[tempi] = min(3.0*(mcmc.temps[tempi] - mcmc.temps[tempi-1])/(tempratio+1.0)*tempratio , fabs(mcmc.temps[tempi]-mcmc.temps[tempi-2]));  //Temperatures of adjacent chains overlap a lot at extrema (since in antiphase), make sure Ti,min>=T-i1,0
      //if(tempi>0) tempampl[tempi] = fabs(mcmc.temps[tempi]-mcmc.temps[tempi-1]);  //Temperatures of adjacent chains overlap: Amplitude = (T(i) - T(i-1))  (may be a bit smallish for large ntemps)
      if(tempi==0 || partemp<=1 || partemp==3) tempampl[tempi] = 0.0;
      if(prpartempinfo>0) printf("     %3d  %7.2lf  %7.2lf  %7.2lf  %7.2lf\n",tempi,mcmc.temps[tempi],tempampl[tempi],mcmc.temps[tempi]-tempampl[tempi],mcmc.temps[tempi]+tempampl[tempi]);
    }
    if(prpartempinfo>0) printf("\n\n");
  }
  if(mcmc.ntemps==1) mcmc.temps[0] = 1.0;
  tempi = 0;  //MUST be zero
  
  
  
  
  
  
  // ********************************************************************************************************************************************************************************
  // ***  WRITE 'HEADER' TO SCREEN AND FILE  ****************************************************************************************************************************************
  // ********************************************************************************************************************************************************************************
  
  
  // *** Print run parameters to screen ***
  if(offsetmcmc==0) printf("   Starting MCMC with the true initial parameters\n\n");
  if(offsetmcmc==1) printf("   Starting MCMC with offset initial parameters\n\n");
  printf("%10s  %10s  %6s  %20s  %6s  ","niter","nburn","seed","null likelihood","ndet");
  for(i=0;i<networksize;i++) {
    printf("%16s%4s  ",ifo[i]->name,"SNR");
  }
  printf("\n%10d  %10d  %6d  %20.10lf  %6d  ",iter,nburn,mcmc.seed,mcmc.logL0,networksize);
  for(i=0;i<networksize;i++) {
    printf("%20.10lf  ",ifo[i]->snr);
  }
  printf("\n");
  
  
  // *** Open the output file and write run parameters in the header ***
  for(tempi=0;tempi<mcmc.ntemps;tempi++) {
    if(tempi==0 || saveallchains==1) {
      sprintf(outfilename,"mcmc.output.%6.6d.%2.2d",mcmc.seed,tempi);
      mcmc.fout = fopen(outfilename,"w"); //In current dir, allows for multiple copies to run
      fprintf(mcmc.fout, "%10s  %10s  %6s  %20s  %6s %8s   %6s  %8s  %10s\n","Niter","Nburn","seed","null likelihood","Ndet","Ncorr","Ntemps","Tmax","Tchain");
      
      fprintf(mcmc.fout, "%10d  %10d  %6d  %20.10lf  %6d %8d   %6d%10d%12.1f\n",iter,nburn,mcmc.seed,mcmc.logL0,networksize,ncorr,mcmc.ntemps,(int)tempmax,mcmc.temps[tempi]);
      fprintf(mcmc.fout, "\n%16s  %16s  %10s  %10s  %10s  %10s  %20s  %15s  %12s  %12s  %12s\n",
	      "Detector","SNR","f_low","f_high","before tc","after tc","Sample start (GPS)","Sample length","Sample rate","Sample size","FT size");
      for(i=0;i<networksize;i++) {
	fprintf(mcmc.fout, "%16s  %16.8lf  %10.2lf  %10.2lf  %10.2lf  %10.2lf  %20.8lf  %15.7lf  %12d  %12d  %12d\n",
		ifo[i]->name,ifo[i]->snr,ifo[i]->lowCut,ifo[i]->highCut,ifo[i]->before_tc,ifo[i]->after_tc,
		ifo[i]->FTstart,ifo[i]->deltaFT,ifo[i]->samplerate,ifo[i]->samplesize,ifo[i]->FTsize);
      }
      fprintf(mcmc.fout, "\n%12s %20s  %32s  %32s  %37s  %32s  %32s  %32s  %32s  %32s  %32s  %32s  %32s  %32s\n",
	      "cycle","logL-logLo","Mc","eta","tc","logdL","spin","kappa","RA","sinlati","phase","sinthJ0","phiJ0","alpha");
      fprintf(mcmc.fout, "%33s  %32s  %32s  %37s  %32s  %32s  %32s  %32s  %32s  %32s  %32s  %32s  %32s\n",
	      " ","parameter     sigma accept","parameter     sigma accept","parameter     sigma accept","parameter     sigma accept","parameter     sigma accept","parameter     sigma accept","parameter     sigma accept","parameter     sigma accept","parameter     sigma accept","parameter     sigma accept","parameter     sigma accept","parameter     sigma accept");
      fclose(mcmc.fout);
    }
  }
  tempi = 0;  //MUST be zero
  
  
  
  
  
  // ********************************************************************************************************************************************************************************
  // ***  INITIALISE MARKOV CHAIN  **************************************************************************************************************************************************
  // ********************************************************************************************************************************************************************************
  
  
  // *** Get true (or best-guess) values for signal ***
  settrueparameters(&state);
  state.loctc    = (double*)calloc(networksize,sizeof(double));
  state.localti  = (double*)calloc(networksize,sizeof(double));
  state.locazi   = (double*)calloc(networksize,sizeof(double));
  state.locpolar = (double*)calloc(networksize,sizeof(double));
  
  
  // *** Write true/best-guess values to screen and file ***
  par2arr(&state, mcmc.param);  //Put the variables in their array
  localpar(&state, ifo, networksize);
  mcmc.logL[tempi] = net_loglikelihood(&state, networksize, ifo);  //Calculate the likelihood
  
  mcmc.iteri = -1;
  mcmc.tempi = tempi;
  write_mcmc_output(mcmc);  //Write to output line to screen and/or file
  tempi = 0;  //MUST be zero

  
  int nparfit=0; //Determine the number of parameters that is actually fitted/varied (i.e. not kept fixed at the true values)
  for(i=0;i<npar;i++) {
    if(fitpar[i]==1) nparfit += 1;
  }
  
  
  // *** Initialise covariance matrix (initially diagonal), to do updates in the first block ***
  corrupdate[0] = 0;
  if(corrupd>0) {
    for(j1=0;j1<npar;j1++) {
      mcmc.covar[tempi][j1][j1] = pdfsigs[j1];
    }
    corrupdate[0] = 1; //Use the matrix above and don't change it
    if(corrupd==2) corrupdate[0] = 2; //Use the matrix above and update it every ncorr iterations
  }
  
  
  
  
  
  // ********************************************************************************************************************************************************************************
  // ***  GET OFFSET STARTING VALUES  ***********************************************************************************************************************************************
  // ********************************************************************************************************************************************************************************
  
  
  //Offset starting values (only for the parameters we're fitting)
  nstart = 0;
  par2arr(&state, mcmc.param);  //Put the variables in their array
  if(offsetmcmc==1) {
    printf("\n");
    for(i=0;i<npar;i++) {
      mcmc.nparam[tempi][i] = mcmc.param[tempi][i];  //Temporarily store the true values
    }
    mcmc.logL[tempi] = -1.e30;
    while(mcmc.logL[tempi] < mcmc.logL0 + 1.0) { //Accept only good starting values
      mcmc.acceptprior[tempi] = 1;
      for(i=0;i<npar;i++) {
	if(fitpar[i]==1 && offsetpar[i]==1) {
	  mcmc.param[tempi][i] = mcmc.nparam[tempi][i] + offsetx * (gsl_rng_uniform(mcmc.ran) - 0.5) * pdfsigs[i];
	  //0:Mc, 1:eta, 2:tc, 3:logd, 4:a, 5:kappa, 6:RA, 7:sindec, 8:phi, 9:sintheta_Jo, 10: phi_Jo, 11:alpha
	  if(i==1 && (mcmc.param[tempi][i]<=0.01 || mcmc.param[tempi][i] > 0.25)) mcmc.param[tempi][i] = max(min(gsl_rng_uniform(mcmc.ran)*0.25,1.0),0.01);  //Eta: 0.01<eta<0.25  \__ If it's that far outside, you may as well take a random value
	  if(i==4 && (mcmc.param[tempi][i]<=1.e-5 || mcmc.param[tempi][i] > 1.0)) mcmc.param[tempi][i] = max(min(gsl_rng_uniform(mcmc.ran),1.0),1.e-5);      //Spin: 0<a<1         /   over the range of this parameter
	  if((i==5 || i==7 || i==9) && (mcmc.param[tempi][i] < -2.0 || mcmc.param[tempi][i] > 2.0)) mcmc.param[tempi][i] = gsl_rng_uniform(mcmc.ran)*2.0 - 1.0;
	  if(i==5 || i==7 || i==9) mcmc.param[tempi][i] = fmod(mcmc.param[tempi][i]+1001.0,2.0) - 1.0;
	  if((i==6 || i==8 || i==10 || i==11) && (mcmc.param[tempi][i] < -2.0*pi || mcmc.param[tempi][i] > 4.0*pi)) mcmc.param[tempi][i] = gsl_rng_uniform(mcmc.ran)*tpi;
	  mcmc.acceptprior[tempi] *= prior(&mcmc.param[tempi][i],i);
	}
      }
      if(mcmc.acceptprior[tempi]==1) {                     //Check the value of the likelihood for this draw
	arr2par(mcmc.param, &state);	                      //Get the parameters from their array
	localpar(&state, ifo, networksize);
	mcmc.logL[tempi] = net_loglikelihood(&state, networksize, ifo);  //Calculate the likelihood
      }
      nstart = nstart + 1;
    }
    printf("%10s  %15s  %8s  %8s  %16s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n", "nDraws","logL-logLo","Mc","eta","tc","logdL","spin","kappa","RA","sinlati","phase","sinthJ0","phiJ0","alpha","logT");
    printf("%10d  %15.6lf  %8.5f  %8.5f  %16.6lf  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.3f\n",
	   nstart,mcmc.logL[tempi]-mcmc.logL0,mcmc.param[tempi][0],mcmc.param[tempi][1],mcmc.param[tempi][2],mcmc.param[tempi][3],mcmc.param[tempi][4],mcmc.param[tempi][5],rightAscension(mcmc.param[tempi][6],GMST(mcmc.param[tempi][2])),mcmc.param[tempi][7],mcmc.param[tempi][8],mcmc.param[tempi][9],mcmc.param[tempi][10],mcmc.param[tempi][11],log10(mcmc.temp));
  }
  
  
  // *** Set the NEW array, sigma and scale ***
  for(i=0;i<npar;i++) {
    mcmc.nparam[tempi][i] = mcmc.param[tempi][i];
    mcmc.sig[tempi][i]   = 0.1  * pdfsigs[i];
    if(adapt==1) mcmc.sig[tempi][i] = pdfsigs[i]; //Don't use adaptation
    mcmc.scale[tempi][i] = 10.0 * pdfsigs[i];
    //mcmc.sig[tempi][i]   = 3.0  * pdfsigs[i]; //3-sigma = delta-100%
    //mcmc.scale[tempi][i] = 0.0 * pdfsigs[i]; //No adaptation
    if(i==6 || i==8 || i==10 || i==11) mcmc.sig[tempi][i] = fmod(mcmc.sig[tempi][i]+mtpi,tpi);  //Bring the sigma between 0 and 2pi
  }
  
  
  
  
  // ********************************************************************************************************************************************************************************
  // ***  WRITE STARTING STATE TO SCREEN AND FILE  **********************************************************************************************************************************
  // ********************************************************************************************************************************************************************************
  
  
  arr2par(mcmc.param, &state);                         //Get the parameters from their array
  localpar(&state, ifo, networksize);
  mcmc.logL[tempi] = net_loglikelihood(&state, networksize, ifo);  //Calculate the likelihood

  // *** Write output line to screen and/or file
  printf("\n");
  mcmc.iteri = 0;
  mcmc.tempi = tempi;
  write_mcmc_output(mcmc);  //Write to output line to screen and/or file
  tempi = 0;  //MUST be zero
  
  
  
  
  // ********************************************************************************************************************************************************************************
  // ***  INITIALISE PARALLEL TEMPERING  ********************************************************************************************************************************************
  // ********************************************************************************************************************************************************************************
  
  // *** Put the initial values of the parameters, sigmas etc in the different temperature chains ***
  if(mcmc.ntemps>1) {
    for(tempi=1;tempi<mcmc.ntemps;tempi++) {
      for(j=0;j<npar;j++) {
	mcmc.param[tempi][j] = mcmc.param[0][j];
	mcmc.nparam[tempi][j] = mcmc.nparam[0][j];
	mcmc.sig[tempi][j] = mcmc.sig[0][j];
	mcmc.scale[tempi][j] = mcmc.scale[0][j];
	mcmc.logL[tempi] = mcmc.logL[0];
	mcmc.nlogL[tempi] = mcmc.nlogL[0];

	for(j1=0;j1<npar;j1++) {
	  for(j2=0;j2<=j1;j2++) {
	    mcmc.covar[tempi][j1][j2] = mcmc.covar[0][j1][j2];
	  }
	}
      }
      corrupdate[tempi] = corrupdate[0];
      //corrupdate[tempi] = 0; //Correlated update proposals only for T=1 chain?
      //corrupdate[mcmc.ntemps-1] = 0; //Correlated update proposals not for hottest chain
    }
  }
  
  
  
  // ********************************************************************************************************************************************************************************
  // ***  CREATE MARKOV CHAIN   *****************************************************************************************************************************************************
  // ********************************************************************************************************************************************************************************
      
  
  iteri = 1;
  while(iteri<=iter) {  //loop over Markov-chain states 
    mcmc.iteri = iteri;
      
    for(tempi=0;tempi<mcmc.ntemps;tempi++) {  //loop over temperature ladder
      mcmc.tempi = tempi;
      //printf(" %d  %d  %d\n",tempi,mcmc.ntemps,iteri);
      
      //Set temperature
      if(partemp==1 || partemp==3) { //Chains at fixed T
	mcmc.temp = mcmc.temps[tempi];
      }
      if(partemp==2 || partemp==4) { //Chains with sinusoid T
	if(tempi==0) {
	  mcmc.temp = 1.0;
	}
	else {
	  //mcmc.temp = mcmc.temps[tempi] * (1.0  +  0.5 * pow((-1.0),tempi) * sin(tpi*(double)iteri/((double)ncorr)));  //Sinusoid around the temperature T_i with amplitude 0.5*T_i and period ncorr
	  //mcmc.temp = mcmc.temps[tempi]  +  tempampl[tempi] * pow((-1.0),tempi) * sin(tpi*(double)iteri/((double)ncorr));  //Sinusoid around the temperature T_i with amplitude tempampl and period ncorr
	  //mcmc.temp = mcmc.temps[tempi]  +  tempampl[tempi] * pow((-1.0),tempi) * sin(tpi*(double)iteri/(0.5*(double)ncorr));  //Sinusoid around the temperature T_i with amplitude tempampl and period 1/2 ncorr
	  mcmc.temp = mcmc.temps[tempi]  +  tempampl[tempi] * pow((-1.0),tempi) * sin(tpi*(double)iteri/(5.0*(double)ncorr));  //Sinusoid around the temperature T_i with amplitude tempampl and period 5 * ncorr
	}
	//printf("%4d  %10.3lf\n",tempi,mcmc.temp);
      }
      
      
      
      // ********************************************************************************************************************************************************************************
      // ***  UPDATE MARKOV CHAIN STATE  ************************************************************************************************************************************************
      // ********************************************************************************************************************************************************************************
      
      // *** Uncorrelated updates
      if(corrupdate[tempi]<=0) {
	//if(iteri>nburn && gsl_rng_uniform(mcmc.ran) < run->blockfrac){                            //Block update; only after burnin
	if(gsl_rng_uniform(mcmc.ran) < run->blockfrac){                                             //Block update; always
	  uncorrelated_mcmc_block_update(*ifo, networksize, &state, &mcmc);
	}
	else{                                                                                       //Componentwise update (90%)
	  uncorrelated_mcmc_single_update(*ifo, networksize, &state, &mcmc);
	}
      } //End uncorrelated updates
      
      
      // *** Correlated updates
      if(corrupdate[tempi]>=1) correlated_mcmc_update(*ifo, networksize, &state, &mcmc);
      // ***
      
      
      
      // Update the dlogL = logL - logLo, and remember the parameter values where it has a maximum
      mcmc.dlogL[tempi] = mcmc.logL[tempi]-mcmc.logL0;
      if(mcmc.dlogL[tempi]>maxdlogL[tempi]) {
	maxdlogL[tempi] = mcmc.dlogL[tempi];
	for(i=0;i<npar;i++) {
	  mcmc.maxparam[tempi][i] = mcmc.param[tempi][i];
	}
      }
      
      
      
      
      
      
      
      //0:mc, 1:eta, 2:tc, 3:logdl, 4:spin, 5:kappa, 6: longi (->RA), 7:sinlati, 8:phase, 9:sinthJ0, 10:phiJ0, 11:alpha
      if(mcmc.acceptprior[0]==1) { //Then write output and care about the correlation matrix
	
	
	// ********************************************************************************************************************************************************************************
	// ***  WRITE STATE TO SCREEN AND FILE  *******************************************************************************************************************************************
	// ********************************************************************************************************************************************************************************
	
	write_mcmc_output(mcmc);  //Write to output line to screen and/or file

	
	
	
	
	// ********************************************************************************************************************************************************************************
	// ***  CORRELATION MATRIX  *******************************************************************************************************************************************************
	// ********************************************************************************************************************************************************************************
	
	
	//if(corrupdate[tempi]==2) {  //Calculate correlations only once
	if(corrupdate[tempi]>=2) { //Calculate correlations multiple times
	  
	  // *** Save state to calculate correlations ***
	  if(ihist[tempi]<ncorr) {
	    for(j1=0;j1<npar;j1++){
	      mcmc.hist[tempi][j1][ihist[tempi]] = mcmc.param[tempi][j1];
	    }
	    ihist[tempi] += 1;
	    sumdlogL[tempi] += mcmc.dlogL[tempi];
	  }
	  
	  
	  // *** Calculate covariances/correlations ***
	  if(ihist[tempi]>=ncorr) {
	    
	    //Calculate the mean
	    for(j1=0;j1<npar;j1++){
	      mu[j1]=0.0;
	      for(i1=0;i1<ncorr;i1++){
		mu[j1]+=mcmc.hist[tempi][j1][i1];
	      }
	      mu[j1]/=((double)ncorr);
	    }
	    
	    //Calculate the standard deviation. Only for printing, not used in the code
	    for(j1=0;j1<npar;j1++){
	      stdev[j1]=0.0;
	      for(i1=0;i1<ncorr;i1++){
		stdev[j1] += (mcmc.hist[tempi][j1][i1]-mu[j1])*(mcmc.hist[tempi][j1][i1]-mu[j1]);
	      }
	      stdev[j1] = sqrt(stdev[j1]/(double)(ncorr-1));
	    }
	    
	    //Calculate the covariances and save them in covar1
	    for(j1=0;j1<npar;j1++){
	      for(j2=0;j2<=j1;j2++){
		covar1[tempi][j1][j2]=0.0;
		for(i1=0;i1<ncorr;i1++){
		  covar1[tempi][j1][j2] += (mcmc.hist[tempi][j1][i1] - mu[j1]) * (mcmc.hist[tempi][j2][i1] - mu[j2]);
		}
		covar1[tempi][j1][j2] /= (double)(ncorr-1);
	      }
	    }
	    
	    //Store the covariance matrix in a temporary variable tempcovar, for Cholesky decomposition
	    for(j1=0;j1<npar;j1++){
	      for(j2=0;j2<=j1;j2++){
		tempcovar[j1][j2] = covar1[tempi][j1][j2];
	      }
	    }
	    
	    if(tempi==0 && prmatrixinfo>0) printf("\n\n");
	    
	    //Do Cholesky decomposition
	    chol(npar,tempcovar);
	    
	    //Get conditions to decide whether to accept the new matrix or not
	    acceptelems[tempi] = 0;
	    for(j1=0;j1<npar;j1++) {
	      if(fitpar[j1]==1) {
		if(tempcovar[j1][j1] < mcmc.covar[tempi][j1][j1]) acceptelems[tempi] += 1; //Smaller diagonal element is better; count for how many this is the case
		if(tempcovar[j1][j1]<=0.0 || isnan(tempcovar[j1][j1])!=0 || isinf(tempcovar[j1][j1])!=0) acceptelems[tempi] -= 9999;  //If diagonal element is <0, NaN or Inf
	      }
	    }
	    acceptelems[tempi] = max(acceptelems[tempi],-1); //Now -1 means there is a diagonal element that is 0, NaN or Inf
	  
	    //Update the matrix only if most of the diagonal elements of covar1 are smaller than the previous ones (and non-zero), to avoid non-positive definite matrices and get improved jump sizes
	    if(prmatrixinfo==2 && tempi==0){  //Print matrix information
	      printf("\n  Update for the covariance matrix proposed at iteration:  %10d\n",iteri);
	      printf("\n    Acceptelems: %d\n",acceptelems[tempi]);
	      printf("\n    Covariance matrix:\n");
	      for(j1=0;j1<npar;j1++){
		for(j2=0;j2<=j1;j2++){
		  printf("    %10.3g",covar1[tempi][j1][j2]);
		}
		printf("\n");
	      }
	      printf("\n    Old Cholesky-decomposed matrix:\n");
	      for(j1=0;j1<npar;j1++){
		for(j2=0;j2<=j1;j2++){
		  printf("    %10.3g",mcmc.covar[tempi][j1][j2]);
		}
		printf("\n");
	      }
	      printf("\n    New Cholesky-decomposed matrix:\n");
	      for(j1=0;j1<npar;j1++){
		for(j2=0;j2<=j1;j2++){
		  //for(j2=0;j2<npar;j2++){
		  printf("    %10.3g",tempcovar[j1][j2]);
		}
		printf("\n");
	      }
	    }
	    
	    // Copy the new covariance matrix from tempcovar into mcmc.covar
	    //if((corrupdate[tempi]==2 && acceptelems[tempi]>0) || (double)acceptelems[tempi] >= (double)nparfit*0.75) { //Only accept new matrix if 75% of diagonal elements are better (smaller)
	    if(acceptelems[tempi]>=0) { //Accept new matrix always, except when Infs, NaNs, 0s occur.
	      for(j1=0;j1<npar;j1++){
		for(j2=0;j2<=j1;j2++){
		  mcmc.covar[tempi][j1][j2] = tempcovar[j1][j2]; 
		}
	      }
	      corrupdate[tempi] += 1;
	      if(prmatrixinfo>0 && tempi==0) printf("  Proposed covariance-matrix update at iteration %d accepted.  Acceptelems: %d.  Accepted matrices: %d/%d \n", iteri, acceptelems[tempi], corrupdate[tempi]-2, (int)((double)iteri/(double)ncorr));  // -2 since you start with 2
	    }
	    else{
	      if(prmatrixinfo>0 && tempi==0) printf("  Proposed covariance-matrix update at iteration %d rejected.  Acceptelems: %d.  Accepted matrices: %d/%d \n", iteri, acceptelems[tempi], corrupdate[tempi]-2, (int)((double)iteri/(double)ncorr));  // -2 since you start with 2
	    }
	    
	    
	    
	    
	    
	    
	    
	    // ********************************************************************************************************************************************************************************
	    // ***  PARALLEL TEMPERING  *******************************************************************************************************************************************************
	    // ********************************************************************************************************************************************************************************
	    
	    
	    if(partemp>=1 && prpartempinfo>0) { //Print info on the current temperature chains
	      if(tempi==0) printf("\n\n      Chain  log(T)   dlog(L)  AccEls AccMat     Swap  AccRat  logStdev:    Mc     eta     t_c  log(d)  a_spin   kappa      RA sin(de)   phi_c theta_J   phi_J alpha_c\n");
	      printf("        %3d   %5.3f %9.2f     %3d    %3d   %6.4f  %6.4f         ",
		     tempi,log10(mcmc.temp),sumdlogL[tempi]/(double)ihist[tempi],acceptelems[tempi],corrupdate[tempi]-2,  (double)swapTs1[tempi]/(double)iteri,(double)mcmc.accepted[tempi][0]/(double)iteri);
	      for(j1=0;j1<npar;j1++) {
		tmpdbl = log10(stdev[j1]+1.e-30);
		if(tmpdbl<-9.99) tmpdbl = 0.0;
		printf("  %6.3f",tmpdbl);
	      }
	      printf("\n");
	      
	      
	      if(prpartempinfo==2 && tempi==mcmc.ntemps-1) { //Print swap-rate matrix for parallel tempering
		printf("\n   Chain swap: ");
		for(j1=0;j1<mcmc.ntemps-1;j1++) {
		  printf("  %6d",j1);
		}
		printf("   total\n");
		for(j1=1;j1<mcmc.ntemps;j1++) {
		  printf("            %3d",j1);
		  for(j2=0;j2<mcmc.ntemps-1;j2++) {
		    if(j2<j1) {
		      printf("  %6.4f",(double)mcmc.swapTss[j2][j1]/(double)iteri);
		    } else {
		      printf("        ");
		    }
		  }
		  printf("  %6.4f\n",(double)swapTs2[j1]/(double)iteri);
		}
		printf("          total");
		for(j1=0;j1<mcmc.ntemps-1;j1++) {
		  printf("  %6.4f",(double)swapTs1[j1]/(double)iteri);
		}
		printf("\n");
	      } //if(prpartempinfo==2 && tempi==mcmc.ntemps-1)
	    } //if(partemp>=1)
	    
	    
	    
	    ihist[tempi] = 0;
	    sumdlogL[tempi] = 0.0;
	    
	    if( (prmatrixinfo>0 || prpartempinfo>0) && tempi==mcmc.ntemps-1 ) printf("\n\n");
	  } //if(ihist[tempi]>=ncorr)
	} //if(corrupdate[tempi]>=2)
	// *************************************************************
      
      } //if(mcmc.acceptprior[tempi]==1)
      
      
      
      
    } // for(tempi=0;tempi<mcmc.ntemps;tempi++) {  //loop over temperature ladder
    
    
    
    
    
    
    // ********************************************************************************************************************************************************************************
    // ***  ANNEALING  ****************************************************************************************************************************************************************
    // ********************************************************************************************************************************************************************************
    
    if(partemp==0) {  //Use only when not using parallel tempering (and of course, temp0>1)
      //mcmc.temp = temp0*pow(1.0/((double)(max(iteri,ncorr))),1.0/10.0);
      //mcmc.temp = temp0*pow(((double)(iteri)),-0.1);
      //mcmc.temp = temp0*pow(((double)(iteri)),-0.25);
      //mcmc.temp = temp0 * pow(((double)(iteri)), (-log(temp0)/log((double)nburn)) );  //Temp drops from temp0 to 1.0 in nburn iterations
      //printf("%f\n",-log(temp0)/log((double)nburn));
      //mcmc.temp = min( max( temp0*pow( max((double)iteri-0.1*(double)nburn,1.0) ,-log10(temp0)/log10(0.9*(double)nburn)) , 1.0) , temp0);  //Temp stays at temp0 for 10% of burnin and then drops to 1.0 at the end of the burnin (too drastic for short burnins/high percentages?)
      //mcmc.temp = temp0*pow( max((double)iteri-1000.0,1.0) ,-log10(temp0)/log10((double)nburn-1000.0));  //Temp stays at temp0 for the first 1000 iterations and then drops to 1.0 at the end of the burnin (too drastic for short burnins/high percentages?)
      mcmc.temp = exp(log(temp0) * ((double)(nburn-iteri)/((double)(nburn-nburn0))));
      
      //if(gsl_rng_uniform(mcmc.ran)<1.e-2 && iteri > nburn && mcmc.dlogL[tempi] < 0.95*maxdlogL[tempi]){
      //if(iteri > nburn && mcmc.dlogL[tempi] < 0.85*maxdlogL[tempi]){
      //  temp=pow(temp0,0.5);
      //  mcmc.corrsig[tempi] *= 10.0;
      //}
      mcmc.temp = min( max(mcmc.temp,1.0) , temp0);
    }
      
    
    
    /*
    // *** Adapt temperature ladder ***
    //  This doesn't really seem to work, since getting a particular likelihood can also been done by zooming in on a piece of parameter space and dropping T; most T's go down to 1 this way.
    if(iteri>=ncorr && ihist[0]>=ncorr-1) {
      printf("\n\n");
      
      //Calculate true and expected average logL for each chain
      for(tempi=0;tempi<mcmc.ntemps;tempi++) {
	avgdlogL[tempi] = sumdlogL[tempi]/(double)ncorr;
	expdlogL[tempi] = maxdlogL[0]/(double)(mcmc.ntemps-1)*(double)(mcmc.ntemps-tempi-1);  //The dlogL you expect for this chain if distributed flat in log(L) between L_0 and L_max
	printf("  tempi:  %d,   T: %10.4f,   maxlodL: %10.4f,   avglogL: %10.4f,   explogL: %10.4f\n",tempi,mcmc.temps[tempi],maxdlogL[tempi],avgdlogL[tempi],expdlogL[tempi]);
      }
      
      printf("\n");
      //See between which two average logL's each expected logL falls
      temp1 = 0;
      temp2 = mcmc.ntemps-1;
      for(tempi=1;tempi<mcmc.ntemps-1;tempi++) {
	for(t=0;t<mcmc.ntemps-1;t++) {
	  if(avgdlogL[t]>=expdlogL[tempi]) temp1 = t;
	}
	for(t=mcmc.ntemps;t>0;t--) {
	  if(avgdlogL[t]<expdlogL[tempi]) temp2 = t;
	}
	
	//Interpolate between temp1 and temp2
	tmpdbl = mcmc.temps[temp1] + (mcmc.temps[temp2]-mcmc.temps[temp1])/(avgdlogL[temp2]-avgdlogL[temp1])*(expdlogL[tempi]-avgdlogL[temp1]);
	printf("  t: %d,  t1-2: %d %d,    avglogL:%10.4f,   t2-t1:%10.4f,  avglogL1-avglogL2:%10.4f tmpdbl:%10.4f\n",tempi,temp1,temp2,avgdlogL[tempi], mcmc.temps[temp2]-mcmc.temps[temp1],  avgdlogL[temp1]-avgdlogL[temp2], tmpdbl );
	
	newtemps[tempi] = tmpdbl;
      } //for(tempi=1;tempi<mcmc.ntemps-1;tempi++)

      for(tempi=1;tempi<mcmc.ntemps-1;tempi++) {
	mcmc.temps[tempi] = newtemps[tempi];
      } //for(tempi=1;tempi<mcmc.ntemps-1;tempi++)
      
      printf("\n\n");
    }
    */
    
    
    
    
    
    
    
    
    
    // ********************************************************************************************************************************************************************************
    // ***  PARALLEL TEMPERING  *******************************************************************************************************************************************************
    // ********************************************************************************************************************************************************************************
    
    // *** Swap states between T-chains ***
    if(mcmc.acceptprior[0]==1) {
      if(partemp>=1 && mcmc.ntemps>1) {
	/*
	//Swap parameters and likelihood between two adjacent chains
	for(tempi=1;tempi<mcmc.ntemps;tempi++) {
	  if(exp(max(-30.0,min(0.0,(1.0/mcmc.temps[tempi-1]-1.0/mcmc.temps[tempi])*(mcmc.logL[tempi]-mcmc.logL[tempi-1])))) > gsl_rng_uniform(mcmc.ran)) { //Then swap...
	    for(i=0;i<npar;i++) {
	      tmpdbl = mcmc.param[tempi][i]; //Temp var
	      mcmc.param[tempi][i] = mcmc.param[tempi-1][i];
	      mcmc.param[tempi-1][i] = tmpdbl;
	    }
	    tmpdbl = mcmc.logL[tempi];
	    mcmc.logL[tempi] = mcmc.logL[tempi-1];
	    mcmc.logL[tempi-1] = tmpdbl;
	    swapTs1[tempi-1] += 1;
	  }
	} //tempi
	*/
	
	//Swap parameters and likelihood between any two chains
	for(tempi=0;tempi<mcmc.ntemps-1;tempi++) {
	  for(tempj=tempi+1;tempj<mcmc.ntemps;tempj++) {
	    
	    if(exp(max(-30.0,min(0.0,(1.0/mcmc.temps[tempi]-1.0/mcmc.temps[tempj])*(mcmc.logL[tempj]-mcmc.logL[tempi])))) > gsl_rng_uniform(mcmc.ran)) { //Then swap...
	      for(i=0;i<npar;i++) {
		tmpdbl = mcmc.param[tempj][i]; //Temp var
		mcmc.param[tempj][i] = mcmc.param[tempi][i];
		mcmc.param[tempi][i] = tmpdbl;
	      }
	      tmpdbl = mcmc.logL[tempj];
	      mcmc.logL[tempj] = mcmc.logL[tempi];
	      mcmc.logL[tempi] = tmpdbl;
	      mcmc.swapTss[tempi][tempj] += 1;
	      swapTs1[tempi] += 1;
	      swapTs2[tempj] += 1;
	    }
	  } //tempj
	} //tempi
	
      } //if(partemp>=1 && mcmc.ntemps>1)
      iteri++;
    } //if(mcmc.acceptprior[0]==1)
    
  } // while(iter<=niter) {  //loop over markov chain states 
  
  
  
  
  // ********************************************************************************************************************************************************************************
  // ***  FREE MEMORY  **************************************************************************************************************************************************************
  // ********************************************************************************************************************************************************************************
  
  //fclose(mcmc.fout);
  printf("\n");
  
  gsl_rng_free(mcmc.ran);
  
  free(corrupdate);
  free(acceptelems);
  
  free(mcmc.temps);
  free(newtemps);
  free(tempampl);
  free(mcmc.logL);
  free(mcmc.nlogL);
  free(mcmc.dlogL);
  free(maxdlogL);
  free(sumdlogL);
  free(avgdlogL);
  free(expdlogL);
  
  free(mcmc.corrsig);
  free(mcmc.acceptprior);
  free(ihist);
  free(swapTs1);
  free(swapTs2);
  
  for(i=0;i<mcmc.ntemps;i++) {
    free(mcmc.accepted[i]);
    free(mcmc.swapTss[i]);
    free(mcmc.param[i]);
    free(mcmc.nparam[i]);
    free(mcmc.maxparam[i]);
    free(mcmc.sig[i]);
    free(mcmc.scale[i]);
  }
  free(mcmc.accepted);
  free(mcmc.swapTss);
  free(mcmc.param);
  free(mcmc.nparam);
  free(mcmc.maxparam);
  free(mcmc.sig);
  free(mcmc.scale);
  
  for(i=0;i<npar;i++) {
    free(tempcovar[i]);
  }
  free(tempcovar);
  
  for(i=0;i<mcmc.ntemps;i++) {
    for(j=0;j<npar;j++) {
      free(mcmc.hist[i][j]);
      free(mcmc.covar[i][j]);
      free(covar1[i][j]);
    }
    free(mcmc.hist[i]);
    free(mcmc.covar[i]);
    free(covar1[i]);
  }
  free(mcmc.hist);
  free(mcmc.covar);
  free(covar1);
  
  free(state.loctc);
  free(state.localti);
  free(state.locazi);
  free(state.locpolar); 
  
  
}
// End adaptive MCMC
//****************************************************************************************************************************************************  





//****************************************************************************************************************************************************  
void chol(int n, double **A)
// Performs cholesky decompositon of A and returns result in the same matrix - adapted from PJG Fortran function
// If matrix is not positive definite, return zeroes
{
  int j1=0,j2=0,j3=0,notposdef=0;
  double sum=0.0;
  for(j1=0;j1<n;j1++){
    if(fitpar[j1]==1) {
      sum = A[j1][j1];
      for(j2=0;j2<j1;j2++){
	if(fitpar[j2]==1) {
	  sum -= A[j1][j2]*A[j1][j2];
	}
      }
      if(sum<0.0) {
	notposdef=1;
      }
      else {
	A[j1][j1]=sqrt(sum);
	for(j2=j1+1;j2<n;j2++){
	  if(fitpar[j2]==1){
	    sum = A[j2][j1];
	    for(j3=0;j3<j1;j3++){
	      if(fitpar[j3]==1){
		sum -= A[j2][j3]*A[j1][j3];
	      }
	    }
	  }
	  A[j2][j1] = sum/A[j1][j1];
	}
      }
    }
  }
  if(notposdef==1) {
    //printf("  Chol(): Matrix %d is not positive definite\n",tempi);
    for(j1=0;j1<n;j1++){
      for(j2=0;j2<n;j2++){
	A[j1][j2] = 0.0;
      }
    }
  }
}
//****************************************************************************************************************************************************  





//****************************************************************************************************************************************************  
void par2arr(struct parset *par, double **param)
//Put the mcmc parameters in their array
//0:mc, 1:eta, 2:tc, 3:logdl, 4:spin, 5:kappa, 6: longi (->RA), 7:sinlati, 8:phase, 9:sinthJ0, 10:phiJ0, 11:alpha
{
  param[tempi][0]  =  par->mc      ;
  param[tempi][1]  =  par->eta     ;
  param[tempi][2]  =  par->tc      ;
  param[tempi][3]  =  par->logdl   ;
  param[tempi][4]  =  par->spin    ;
  param[tempi][5]  =  par->kappa   ;
  param[tempi][6]  =  par->longi   ;
  param[tempi][7]  =  par->sinlati ;
  param[tempi][8]  =  par->phase   ;
  param[tempi][9]  =  par->sinthJ0 ;
  param[tempi][10] =  par->phiJ0   ;
  param[tempi][11] =  par->alpha   ;
}
//****************************************************************************************************************************************************  




//****************************************************************************************************************************************************  
void arr2par(double **param, struct parset *par)
//Get the mcmc parameters from their array
//0:mc, 1:eta, 2:tc, 3:logdl, 4:spin, 5:kappa, 6: longi (->RA), 7:sinlati, 8:phase, 9:sinthJ0, 10:phiJ0, 11:alpha
{
  par->mc      =  param[tempi][0]   ;
  par->eta     =  param[tempi][1]   ;
  par->tc      =  param[tempi][2]   ;
  par->logdl   =  param[tempi][3]   ;
  par->spin    =  param[tempi][4]   ;
  par->kappa   =  param[tempi][5]   ;
  par->longi   =  param[tempi][6]   ;
  par->sinlati =  param[tempi][7]   ;
  par->phase   =  param[tempi][8]   ;
  par->sinthJ0 =  param[tempi][9]   ;
  par->phiJ0   =  param[tempi][10]  ;
  par->alpha   =  param[tempi][11]  ;
}
//****************************************************************************************************************************************************  








//****************************************************************************************************************************************************  
int prior(double *par, int p)
//****************************************************************************************************************************************************  
//Contains boundary conditions and prior information for the adaptive MCMC.  Trying to avoid returning 0, to increase jump sizes
//  'Stick to the wall method' should not be used, since it is asymmetric
//  0:mc, 1:eta, 2:tc, 3:logdl, 4:spin, 5:kappa, 6: longi (->RA), 7:sinlati, 8:phase, 9:sinthJ0, 10:phiJ0, 11:alpha
{
  int prior = 1;
  double dt = 0.0;
  double *lb,*ub;
  lb = (double*)calloc(12,sizeof(double));
  ub = (double*)calloc(12,sizeof(double));
  
  //Lower and upper boundaries:
  lb[0] = 1.0; //Mc
  ub[0] = 6.0;
  
  lb[1] = 0.03; //eta
  ub[1] = 0.245;
  
  dt = 0.05; //This should be dt/2  For known signals
  //dt = 0.5; //This should be dt/2  For unknown signals
  lb[2] = prior_tc_mean - dt; //t_c
  ub[2] = prior_tc_mean + dt;
  
  lb[3] = -6.9; //log(d_L)
  ub[3] = 4.6;
  
  lb[4] = 1.e-10; //a_spin
  ub[4] = 0.999999;
  
  lb[5] = -0.999999; //kappa
  ub[5] = 0.999999;
  
  lb[7] = -0.999999; //sin(dec)
  ub[7] = 0.999999;
  
  lb[9] = -0.999999; //sin(theta_J0)
  ub[9] = 0.999999;
  
  //Set prior to 0 if outside range
  /*
    if(p==0) if(*par < 1.0 || *par > 5.0) prior = 0;                // Chirp mass <1 or >5
    if(p==1) if(*par < 0.03 || *par > 0.245) prior = 0;                // Mass ratio, refect if <0.03 or >0.245
    if(p==2) if(*par <= prior_tc_mean -0.05 || *par > prior_tc_mean+0.05) prior = 0;    // If time outside range
    if(p==3) if(*par <= -6.91 || *par > 4.6) prior = 0;                  // If distance <1kpc or >100Mpc
    if(p==4) if(*par <= 0.0 || *par > 1.0) prior = 0;                     // If spin<0 or >1
    if(p==5 || p==7 || p==9)  if(*par < -1.0 || *par > 1.0) prior = 0;    // If sin or cos <-1 or >1
  */
  
  //printf("%2d  %20.6f  %20.6f  %20.6f  ",p,lb[p],ub[p],*par);
  
  //Bounce back from the wall:
  //if((p==0 || p==1 || p==2 || p==3 || p==4 || p==5 || p==7 || p==9) && *par<lb[p])  *par = lb[p] + fabs(*par - lb[p]); //Does only 1 step
  //if((p==0 || p==1 || p==2 || p==3 || p==4 || p==5 || p==7 || p==9) && *par>ub[p])  *par = ub[p] - fabs(*par - ub[p]);
  if(p==0 || p==1 || p==2 || p==3 || p==4 || p==5 || p==7 || p==9) {
    while(*par<lb[p] || *par>ub[p]) {
      while(*par<lb[p])  *par = lb[p] + fabs(*par - lb[p]);
      while(*par>ub[p])  *par = ub[p] - fabs(*par - ub[p]);
    }
  }
  
  if(p==6 || p==8 || p==10 || p==11) *par = fmod(*par+mtpi,tpi);  //Periodic boundary condition to bring the variable between 0 and 2pi
  
  //printf("%20.6f  %d  \n",*par,prior);
  
  free(lb);
  free(ub);
  
  return prior;
}
//****************************************************************************************************************************************************  






//Do a correlated block update
//****************************************************************************************************************************************************  
void correlated_mcmc_update(struct interferometer ifo[], int networksize, struct parset *state, struct mcmcvariables *mcmc)
//****************************************************************************************************************************************************  
{
  int j1=0, j2=0, tempi=mcmc->tempi;
  double work[mcmc->npar], dparam=0.0;
  
  for(j1=0;j1<mcmc->npar;j1++) {
    work[j1] = gsl_ran_gaussian(mcmc->ran,1.0) * mcmc->corrsig[tempi];   //Univariate gaussian random number, with sigma=1, times the adaptable sigma_correlation
  }
  mcmc->acceptprior[tempi] = 1;
  for(j1=0;j1<mcmc->npar;j1++){
    if(fitpar[j1]==1) {
      dparam = 0.0;
      for(j2=0;j2<=j1;j2++){
	dparam += mcmc->covar[tempi][j1][j2]*work[j2];   //Work is now a univariate gaussian random number
      }
      mcmc->nparam[tempi][j1] = mcmc->param[tempi][j1] + dparam;       //Jump from the previous parameter value
      mcmc->sig[tempi][j1] = fabs(dparam);                       //This isn't really sigma, but the proposed jump size
      mcmc->acceptprior[tempi] *= prior(&mcmc->nparam[tempi][j1],j1);
    }
  }
  
  if(mcmc->acceptprior[tempi]==1) {                                    //Then calculate the likelihood
    arr2par(mcmc->nparam, state);	                      //Get the parameters from their array
    localpar(state, &ifo, networksize);
    mcmc->nlogL[tempi] = net_loglikelihood(state, networksize, &ifo); //Calculate the likelihood
    par2arr(state, mcmc->nparam);	                      //Put the variables back in their array
    
    if(exp(max(-30.0,min(0.0,mcmc->nlogL[tempi]-mcmc->logL[tempi]))) > pow(gsl_rng_uniform(mcmc->ran),mcmc->temp) && mcmc->nlogL[tempi] > mcmc->logL0) {  //Accept proposal
      for(j1=0;j1<mcmc->npar;j1++){
	if(fitpar[j1]==1) {
	  mcmc->param[tempi][j1] = mcmc->nparam[tempi][j1];
	  mcmc->accepted[tempi][j1] += 1;
	}
      }
      mcmc->logL[tempi] = mcmc->nlogL[tempi];
      if(adapt==1){ 
	mcmc->corrsig[tempi] *= 10.0;  //Increase sigma
      }
    }
    else {                                                      //Reject proposal because of low likelihood
      if(adapt==1){ 
	mcmc->corrsig[tempi] *= 0.5;  //Decrease sigma
      }
    }
  }
  else {                                                      //Reject proposal because of boundary conditions.  Perhaps one should increase the step size, or at least not decrease it?
    /*
      if(adapt==1){ 
      mcmc->corrsig[tempi] *= 0.8;
      }
    */
  }
}
//****************************************************************************************************************************************************  












//Do an uncorrelated single-parameter update
//****************************************************************************************************************************************************  
void uncorrelated_mcmc_single_update(struct interferometer ifo[], int networksize, struct parset *state, struct mcmcvariables *mcmc)
//****************************************************************************************************************************************************  
{
  int j1=0;
  double gamma=0.0,alphastar=0.25;
  
  for(j1=0;j1<npar;j1++){
    if(fitpar[j1]==1) mcmc->nparam[tempi][j1] = mcmc->param[tempi][j1];
  }
  for(j1=0;j1<npar;j1++){
    if(fitpar[j1]==1) {
      mcmc->nparam[tempi][j1] = mcmc->param[tempi][j1] + gsl_ran_gaussian(mcmc->ran,mcmc->sig[tempi][j1]);
      
      mcmc->acceptprior[tempi] = prior(&mcmc->nparam[tempi][j1],j1);
      if(mcmc->acceptprior[tempi]==1) {
	arr2par(mcmc->nparam, state);                          //Get the parameters from their array
	localpar(state, &ifo, networksize);
	mcmc->nlogL[tempi] = net_loglikelihood(state, networksize, &ifo);   //Calculate the likelihood
	par2arr(state, mcmc->nparam);                          //Put the variables back in their array
	
	if(exp(max(-30.0,min(0.0,mcmc->nlogL[tempi]-mcmc->logL[tempi]))) > pow(gsl_rng_uniform(mcmc->ran),mcmc->temp) && mcmc->nlogL[tempi] > mcmc->logL0) {  //Accept proposal
	  mcmc->param[tempi][j1] = mcmc->nparam[tempi][j1];
	  mcmc->logL[tempi] = mcmc->nlogL[tempi];
	  if(adapt==1){
	    gamma = mcmc->scale[tempi][j1]*pow(1.0/((double)(mcmc->iteri+1)),1.0/6.0);
	    mcmc->sig[tempi][j1] = max(0.0,mcmc->sig[tempi][j1] + gamma*(1.0 - alphastar)); //Accept
	    if(j1==6 || j1==8 || j1==10 || j1==11) mcmc->sig[tempi][j1] = fmod(mcmc->sig[tempi][j1]+mtpi,tpi);  //Bring the sigma between 0 and 2pi
	  }
	  mcmc->accepted[tempi][j1] += 1;
	}
	else{                                                      //Reject proposal
	  mcmc->nparam[tempi][j1] = mcmc->param[tempi][j1];
	  if(adapt==1){
	    gamma = mcmc->scale[tempi][j1]*pow(1.0/((double)(mcmc->iteri+1)),1.0/6.0);
	    mcmc->sig[tempi][j1] = max(0.0,mcmc->sig[tempi][j1] - gamma*alphastar); //Reject
	    if(j1==6 || j1==8 || j1==10 || j1==11) mcmc->sig[tempi][j1] = fmod(mcmc->sig[tempi][j1]+mtpi,tpi);  //Bring the sigma between 0 and 2pi
	    //mcmc->sig[tempi][j1] = max(0.01*mcmc->sig[tempi][j1], mcmc->sig[tempi][j1] - gamma*alphastar);
	  }
	}
      }
      else{  //If new state not within boundaries
	mcmc->nparam[tempi][j1] = mcmc->param[tempi][j1];
	if(adapt==1) {
	  gamma = mcmc->scale[tempi][j1]*pow(1.0/((double)(mcmc->iteri+1)),1.0/6.0);
	  mcmc->sig[tempi][j1] = max(0.0,mcmc->sig[tempi][j1] - gamma*alphastar); //Reject
	  if(j1==6 || j1==8 || j1==10 || j1==11) mcmc->sig[tempi][j1] = fmod(mcmc->sig[tempi][j1]+mtpi,tpi);  //Bring the sigma between 0 and 2pi
	}
      } //if(mcmc->acceptprior[tempi]==1)
    } //if(fitpar[j1]==1)
  } //j1
}
//****************************************************************************************************************************************************  

	




//Do an uncorrelated block update
//****************************************************************************************************************************************************  
void uncorrelated_mcmc_block_update(struct interferometer ifo[], int networksize, struct parset *state, struct mcmcvariables *mcmc)
//****************************************************************************************************************************************************  
{
  int j1=0;
  
  mcmc->acceptprior[tempi] = 1;
  for(j1=0;j1<mcmc->npar;j1++){
    if(fitpar[j1]==1) {
      mcmc->nparam[tempi][j1] = mcmc->param[tempi][j1] + gsl_ran_gaussian(mcmc->ran,mcmc->sig[tempi][j1]);
      mcmc->acceptprior[tempi] *= prior(&mcmc->nparam[tempi][j1],j1);
    }
  }
  
  if(mcmc->acceptprior[tempi]==1) {
    arr2par(mcmc->nparam, state);	                              //Get the parameters from their array
    localpar(state, &ifo, networksize);                               //Calculate local variables
    mcmc->nlogL[tempi] = net_loglikelihood(state, networksize, &ifo);  //Calculate the likelihood
    par2arr(state, mcmc->nparam);	                              //Put the variables back in their array
    
    if(exp(max(-30.0,min(0.0,mcmc->nlogL[tempi]-mcmc->logL[tempi]))) > pow(gsl_rng_uniform(mcmc->ran),mcmc->temp) && mcmc->nlogL[tempi] > mcmc->logL0){  //Accept proposal if L>Lo
      for(j1=0;j1<mcmc->npar;j1++){
	if(fitpar[j1]==1) {
	  mcmc->param[tempi][j1] = mcmc->nparam[tempi][j1];
	  mcmc->accepted[tempi][j1] += 1;
	}
      }
      mcmc->logL[tempi] = mcmc->nlogL[tempi];
    }
  }
}
//****************************************************************************************************************************************************  



	
//Write an MCMC line to screen and/or file
//****************************************************************************************************************************************************  
void write_mcmc_output(struct mcmcvariables mcmc)
//****************************************************************************************************************************************************  
{
  int j1=0, tempi=mcmc.tempi, iteri=mcmc.iteri;
  char outfilename[99];
  
  //printf("%d  %d",tempi,iteri);
  
  // *** Write output to screen ***
  if(tempi==0) { //Only for the T=1 chain
    if((iteri % (50*screenoutput))==0 || iteri<0)  printf("\n%10s  %15s  %8s  %8s  %16s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n",
							  "cycle","logL-logLo","Mc","eta","tc","logdL","spin","kappa","RA","sinlati","phase","sinthJ0","phiJ0","alpha","logT");
    if((iteri % screenoutput)==0 || iteri<0)  printf("%10d  %15.6lf  %8.5f  %8.5f  %16.6lf  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.3f\n",
						     iteri,mcmc.logL[tempi]-mcmc.logL0,mcmc.param[tempi][0],mcmc.param[tempi][1],mcmc.param[tempi][2],mcmc.param[tempi][3],mcmc.param[tempi][4],mcmc.param[tempi][5],rightAscension(mcmc.param[tempi][6],GMST(mcmc.param[tempi][2])),mcmc.param[tempi][7],mcmc.param[tempi][8],mcmc.param[tempi][9],mcmc.param[tempi][10],mcmc.param[tempi][11],log10(mcmc.temp));
  }
  
  
  double *accrat;
  accrat = (double*)calloc(mcmc.npar,sizeof(double));
  for(j1=0;j1<mcmc.npar;j1++) {
    accrat[j1] = 0.0;
  }
  if(iteri > 0) {
    for(j1=0;j1<mcmc.npar;j1++) {
      accrat[j1] = mcmc.accepted[tempi][j1]/(double)iteri;
    }
  }
  
  // *** Write output to file ***
  if(tempi==0 || saveallchains==1) { //For all T-chains if desired, otherwise the T=1 chain only
    if((iteri % skip)==0){
      sprintf(outfilename,"mcmc.output.%6.6d.%2.2d",mcmc.seed,tempi);
      mcmc.fout = fopen(outfilename,"a");
      fprintf(mcmc.fout, "%12d %20.10lf  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %20.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f\n",
	      iteri,mcmc.logL[tempi]-mcmc.logL0,
	      /*
	      mcmc.param[tempi][0],mcmc.sig[tempi][0],mcmc.accepted[tempi][0]/(double)iteri,  mcmc.param[tempi][1],mcmc.sig[tempi][1],mcmc.accepted[tempi][1]/(double)iteri,  
	      mcmc.param[tempi][2],mcmc.sig[tempi][2],mcmc.accepted[tempi][2]/(double)iteri,  mcmc.param[tempi][3],mcmc.sig[tempi][3],mcmc.accepted[tempi][3]/(double)iteri, 
	      mcmc.param[tempi][4],mcmc.sig[tempi][4],mcmc.accepted[tempi][4]/(double)iteri,  mcmc.param[tempi][5],mcmc.sig[tempi][5],mcmc.accepted[tempi][5]/(double)iteri,
	      rightAscension(mcmc.param[tempi][6],GMST(mcmc.param[tempi][2])),mcmc.sig[tempi][6],mcmc.accepted[tempi][6]/(double)iteri,
	      mcmc.param[tempi][7],mcmc.sig[tempi][7],mcmc.accepted[tempi][7]/(double)iteri,  mcmc.param[tempi][8],mcmc.sig[tempi][8],mcmc.accepted[tempi][8]/(double)iteri,
	      mcmc.param[tempi][9],mcmc.sig[tempi][9],mcmc.accepted[tempi][9]/(double)iteri,  
	      mcmc.param[tempi][10],mcmc.sig[tempi][10],mcmc.accepted[tempi][10]/(double)iteri,  mcmc.param[tempi][11],mcmc.sig[tempi][11],mcmc.accepted[tempi][11]/(double)iteri);
	      */
	      mcmc.param[tempi][0],mcmc.sig[tempi][0],accrat[0],  mcmc.param[tempi][1],mcmc.sig[tempi][1],accrat[1],  
	      mcmc.param[tempi][2],mcmc.sig[tempi][2],accrat[2],  mcmc.param[tempi][3],mcmc.sig[tempi][3],accrat[3], 
	      mcmc.param[tempi][4],mcmc.sig[tempi][4],accrat[4],  mcmc.param[tempi][5],mcmc.sig[tempi][5],accrat[5],
	      rightAscension(mcmc.param[tempi][6],GMST(mcmc.param[tempi][2])),mcmc.sig[tempi][6],accrat[6],
	      mcmc.param[tempi][7],mcmc.sig[tempi][7],accrat[7],  mcmc.param[tempi][8],mcmc.sig[tempi][8],accrat[8],
	      mcmc.param[tempi][9],mcmc.sig[tempi][9],accrat[9],  
	      mcmc.param[tempi][10],mcmc.sig[tempi][10],accrat[10],  mcmc.param[tempi][11],mcmc.sig[tempi][11],accrat[11]);
      
      //fflush(mcmc.fout); //Make sure any 'snapshot' you take halfway is complete
      fclose(mcmc.fout);
    }
  } //if(tempi==0)
  
  free(accrat);
}

  
