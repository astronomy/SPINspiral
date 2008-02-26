#include <mcmc.h>


void mcmc(int networksize, struct interferometer *ifo[])
{
  if(MvdSdebug) printf("MCMC\n");
  struct parset state;
  //const gsl_rng_type*T;
  gsl_rng *r; 
  r = gsl_rng_alloc(gsl_rng_mt19937);
  FILE *fout;
  char outfilename[99]; 
  
  if(partemp==0) ntemps=1;
  if(nburn0>=nburn) {
    printf("\n\n   *** Warning: nburn0 > nburn, setting nburn0 = nburn*0.9 ***\n\n");
    nburn0 = (int)(0.9*(double)nburn);
  }
  
  
  
  // ********************************************************************************************************************************************************************************
  // ***  MEMORY ALLOCATION  ********************************************************************************************************************************************************
  // ********************************************************************************************************************************************************************************
  
  //const int npar=12; //Make it a global variable
  int i=0,i1=0,j=0,j1=0,j2=0,iteri=0,tempj=0;
  double scale[npar],gamma=0.0,alphastar=0.25,temp=temp0;
  double work[npar];
  
  double dparam=0.0,tmpdbl=0.0,tempratio=1.0;
  int nstart=0;
  double mu[npar],stdev[npar];
  
  int *corrupdate,*acceptelems;
  corrupdate = (int*)calloc(ntemps,sizeof(int));
  acceptelems = (int*)calloc(ntemps,sizeof(int));
  for(i=0;i<ntemps;i++) {
    corrupdate[i] = 0;
    acceptelems[i] = 0;
  }
  
  double *temps,*newtemps,*tempampl,*logL,*nlogL,*dlogL,*maxdlogL,*sumdlogL,*avgdlogL,*expdlogL;
  temps = (double*)calloc(ntemps,sizeof(double));
  newtemps = (double*)calloc(ntemps,sizeof(double));
  tempampl = (double*)calloc(ntemps,sizeof(double));
  logL = (double*)calloc(ntemps,sizeof(double));
  nlogL = (double*)calloc(ntemps,sizeof(double));
  dlogL = (double*)calloc(ntemps,sizeof(double));
  maxdlogL = (double*)calloc(ntemps,sizeof(double));
  sumdlogL = (double*)calloc(ntemps,sizeof(double));
  avgdlogL = (double*)calloc(ntemps,sizeof(double));
  expdlogL = (double*)calloc(ntemps,sizeof(double));
  for(i=0;i<ntemps;i++) {
    temps[i] = 0.0;
    newtemps[i] = 0.0;
    tempampl[i] = 0.0;
    logL[i] = 0.0;
    nlogL[i] = 0.0;
    dlogL[i] = 0.0;
    maxdlogL[i] = -1.e30;
    sumdlogL[i] = 0.0;
    avgdlogL[i] = 0.0;
    expdlogL[i] = 0.0;
  }
  
  double *corrsig,*corrscale;
  int *acceptprior,*ihist,*swapTs1,*swapTs2;
  corrsig = (double*)calloc(ntemps,sizeof(double));
  corrscale = (double*)calloc(ntemps,sizeof(double));
  swapTs1 = (int*)calloc(ntemps,sizeof(int));
  swapTs2 = (int*)calloc(ntemps,sizeof(int));
  acceptprior = (int*)calloc(ntemps,sizeof(int));
  ihist = (int*)calloc(ntemps,sizeof(int));
  for(i=0;i<ntemps;i++) {
    corrsig[i] = 1.0;
    corrscale[i] = 10.0;
    swapTs1[i] = 0;   //Number of swaps between parallel chains
    swapTs2[i] = 0;   //Number of swaps between parallel chains
    acceptprior[i] = 1;
    ihist[i] = 0;
  }
  
  int **accepted,**swapTss;
  double **param,**nparam,**maxparam,**sig;
  accepted = (int**)calloc(ntemps,sizeof(int*)); // 
  swapTss = (int**)calloc(ntemps,sizeof(int*)); // 
  param = (double**)calloc(ntemps,sizeof(double*)); // The old parameters for all chains
  nparam = (double**)calloc(ntemps,sizeof(double*)); // The new parameters for all chains
  sig = (double**)calloc(ntemps,sizeof(double*)); // The standard deviation of the gaussian to draw the jump size from
  maxparam = (double**)calloc(ntemps,sizeof(double*)); // The best parameters for all chains (max logL)
  for(i=0;i<ntemps;i++) {
    accepted[i] = (int*)calloc(npar,sizeof(int));
    swapTss[i] = (int*)calloc(ntemps,sizeof(int));
    param[i] = (double*)calloc(npar,sizeof(double));
    nparam[i] = (double*)calloc(npar,sizeof(double));
    sig[i] = (double*)calloc(npar,sizeof(double));
    maxparam[i] = (double*)calloc(npar,sizeof(double));
  }
  
  double **covar2a;
  covar2a = (double**)calloc(npar,sizeof(double*)); // A temp Cholesky-decomposed matrix
  for(i=0;i<npar;i++) {
    covar2a[i] = (double*)calloc(npar,sizeof(double));
  }
  
  double ***hist,***covar1,***covar2;
  hist    = (double***)calloc(ntemps,sizeof(double**)); // Store past iterations to calculate the covariances
  covar1  = (double***)calloc(ntemps,sizeof(double**)); // The actual covariance matrix
  covar2  = (double***)calloc(ntemps,sizeof(double**)); // The Cholesky-decomposed matrix
  for(i=0;i<ntemps;i++) {
    hist[i]    = (double**)calloc(npar,sizeof(double*));
    covar1[i]  = (double**)calloc(npar,sizeof(double*));
    covar2[i]  = (double**)calloc(npar,sizeof(double*));
    for(j=0;j<npar;j++) {
      hist[i][j]    = (double*)calloc(ncorr,sizeof(double));
      covar1[i][j]  = (double*)calloc(npar,sizeof(double));
      covar2[i][j]  = (double*)calloc(npar,sizeof(double));
    }
  }
  
  
  
  // ********************************************************************************************************************************************************************************
  // ***  INITIALISE PARALLEL TEMPERING  ********************************************************************************************************************************************
  // ********************************************************************************************************************************************************************************
  
  // *** Set up temperature ladder ***
  if(ntemps>1) {  
    tempratio = exp(log(tempmax)/(double)(ntemps-1));
    if(prpartempinfo>0) {
      printf("   Temperature ladder:\n     Number of chains:%3d,  Tmax:%7.2lf, Ti/Ti-1:%7.3lf\n",ntemps,tempmax,tempratio);
      if(partemp==1) printf("     Using fixed temperatures for the chains\n");
      if(partemp==2) printf("     Using sinusoid temperatures for the chains\n");
      if(partemp==3) printf("     Using a manual temperature ladder with fixed temperatures for the chains\n");
      if(partemp==4) printf("     Using a manual temperature ladder with sinusoid temperatures for the chains\n");
      printf("     Chain     To     Ampl.    Tmin     Tmax\n");
    }
    for(tempi=0;tempi<ntemps;tempi++) {
      temps[tempi] = pow(10.0,log10(tempmax)/(double)(ntemps-1)*(double)tempi);
      if(partemp==3 || partemp==4) {  //Set chains manually; for Tmax = 50 and ntemps=6
	temps[0] = 1.0;
	temps[1] = 2.66;
	temps[2] = 7.07;
	temps[3] = 11.0; //Extra T between 2 and 4, the rest is for 5 chains between 1 and 50
	temps[4] = 18.80;
	temps[5] = 50.0;
      }
      //if(tempi>0) tempampl[tempi] = (temps[tempi] - temps[tempi-1])/(tempratio+1.0)*tempratio;  //Temperatures of adjacent chains just touch at extrema (since in antiphase)
      //if(tempi>0) tempampl[tempi] = 1.5*(temps[tempi] - temps[tempi-1])/(tempratio+1.0)*tempratio;  //Temperatures of adjacent chains overlap somewhat at extrema (since in antiphase)
      if(tempi>0) tempampl[tempi] = min(3.0*(temps[tempi] - temps[tempi-1])/(tempratio+1.0)*tempratio , fabs(temps[tempi]-temps[tempi-1]));  //Temperatures of adjacent chains overlap a lot at extrema (since in antiphase), make sure Ti,min>=T-i1,0
      if(ntemps>10 && tempi>1) tempampl[tempi] = min(3.0*(temps[tempi] - temps[tempi-1])/(tempratio+1.0)*tempratio , fabs(temps[tempi]-temps[tempi-2]));  //Temperatures of adjacent chains overlap a lot at extrema (since in antiphase), make sure Ti,min>=T-i1,0
      //if(tempi>0) tempampl[tempi] = fabs(temps[tempi]-temps[tempi-1]);  //Temperatures of adjacent chains overlap: Amplitude = (T(i) - T(i-1))  (may be a bit smallish for large ntemps)
      if(tempi==0 || partemp<=1 || partemp==3) tempampl[tempi] = 0.0;
      if(prpartempinfo>0) printf("     %3d  %7.2lf  %7.2lf  %7.2lf  %7.2lf\n",tempi,temps[tempi],tempampl[tempi],temps[tempi]-tempampl[tempi],temps[tempi]+tempampl[tempi]);
    }
    if(prpartempinfo>0) printf("\n\n");
  }
  if(ntemps==1) temps[0] = 1.0;
  tempi = 0;  //MUST be zero
  
  
  
  
  
  
  // ********************************************************************************************************************************************************************************
  // ***  WRITE INITIAL STATE TO SCREEN AND FILE  ***********************************************************************************************************************************
  // ********************************************************************************************************************************************************************************
  
  
  // *** Print run parameters to screen ***
  if(offsetmcmc==0) printf("   Starting MCMC with the true initial parameters\n\n");
  if(offsetmcmc==1) printf("   Starting MCMC with offset initial parameters\n\n");
  printf("%10s  %10s  %6s  %20s  %6s  ","niter","nburn","seed","null likelihood","ndet");
  for(i=0;i<networksize;i++) {
    printf("%16s%4s  ",ifo[i]->name,"SNR");
  }
  printf("\n%10d  %10d  %6d  %20.10lf  %6d  ",iter,nburn,mcmcseed,NullLikelihood,networksize);
  for(i=0;i<networksize;i++) {
    printf("%20.10lf  ",ifo[i]->snr);
  }
  printf("\n\n%10s  %15s  %8s  %8s  %16s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n","cycle","logL-logLo","Mc","eta","tc","logdL","spin","kappa","RA","sinlati","phase","sinthJ0","phiJ0","alpha","logT");
  
  
  // *** Open the output file and write run parameters in the header ***
  for(tempi=0;tempi<ntemps;tempi++) {
    if(tempi==0 || saveallchains==1) {
      sprintf(outfilename,"mcmc.output.%6.6d.%2.2d",mcmcseed,tempi);
      fout = fopen(outfilename,"w"); //In current dir, allows for multiple copies to run
      fprintf(fout, "%10s  %10s  %6s  %20s  %6s %8s   %6s  %8s  %10s\n","Niter","Nburn","seed","null likelihood","Ndet","Ncorr","Ntemps","Tmax","Tchain");
      
      fprintf(fout, "%10d  %10d  %6d  %20.10lf  %6d %8d   %6d%10d%12.1f\n",iter,nburn,mcmcseed,NullLikelihood,networksize,ncorr,ntemps,(int)tempmax,temps[tempi]);
      fprintf(fout, "\n%16s  %16s  %10s  %10s  %10s  %10s  %20s  %15s  %12s  %12s  %12s\n",
	      "Detector","SNR","f_low","f_high","before tc","after tc","Sample start (GPS)","Sample length","Sample rate","Sample size","FT size");
      for(i=0;i<networksize;i++) {
	fprintf(fout, "%16s  %16.8lf  %10.2lf  %10.2lf  %10.2lf  %10.2lf  %20.8lf  %15.7lf  %12d  %12d  %12d\n",
		ifo[i]->name,ifo[i]->snr,ifo[i]->lowCut,ifo[i]->highCut,ifo[i]->before_tc,ifo[i]->after_tc,
		ifo[i]->FTstart,ifo[i]->deltaFT,ifo[i]->samplerate,ifo[i]->samplesize,ifo[i]->FTsize);
      }
      fprintf(fout, "\n%12s %20s  %32s  %32s  %37s  %32s  %32s  %32s  %32s  %32s  %32s  %32s  %32s  %32s\n",
	      "cycle","logL-logLo","Mc","eta","tc","logdL","spin","kappa","RA","sinlati","phase","sinthJ0","phiJ0","alpha");
      fprintf(fout, "%33s  %32s  %32s  %37s  %32s  %32s  %32s  %32s  %32s  %32s  %32s  %32s  %32s\n",
	      " ","parameter     sigma accept","parameter     sigma accept","parameter     sigma accept","parameter     sigma accept","parameter     sigma accept","parameter     sigma accept","parameter     sigma accept","parameter     sigma accept","parameter     sigma accept","parameter     sigma accept","parameter     sigma accept","parameter     sigma accept");
      fclose(fout);
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
  par2arr(&state, param);  //Put the variables in their array
  localpar(&state, ifo, networksize);
  logL[tempi] = net_loglikelihood(&state, networksize, ifo);  //Calculate the likelihood
  
  printf("%10d  %15.6lf  %8.5f  %8.5f  %16.6lf  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.3f\n",
	 -1,logL[tempi]-NullLikelihood,param[tempi][0],param[tempi][1],param[tempi][2],param[tempi][3],param[tempi][4],param[tempi][5],rightAscension(param[tempi][6],GMST(param[tempi][2])),param[tempi][7],param[tempi][8],param[tempi][9],param[tempi][10],param[tempi][11],log10(temp));
  
  for(tempi=0;tempi<ntemps;tempi++) {
    if(tempi==0 || saveallchains==1) {
      sprintf(outfilename,"mcmc.output.%6.6d.%2.2d",mcmcseed,tempi);
      fout = fopen(outfilename,"a");
      fprintf(fout, "%12d %20.10lf  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %20.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f\n",
	      -1,logL[0]-NullLikelihood,param[0][0],0.0,0.0,param[0][1],0.0,0.0,param[0][2],0.0,0.0,param[0][3],0.0,0.0,param[0][4],0.0,0.0,param[0][5],0.0,0.0,rightAscension(param[0][6],GMST(param[0][2])),0.0,0.0,param[0][7],0.0,0.0,
	      param[0][8],0.0,0.0,param[0][9],0.0,0.0,param[0][10],0.0,0.0,param[0][11],0.0,0.0);
      fclose(fout);
    }
  }
  tempi = 0;  //MUST be zero

  
  //Set seed for this chain
  gsl_rng_set(r, mcmcseed);  
  
  int nparfit=0; //Determine the number of parameters that is actually fitted/varied (i.e. not kept fixed at the true values)
  for(i=0;i<npar;i++) {
    if(fitpar[i]==1) nparfit += 1;
  }
  
  
  // *** Initialise covariance matrix (initially diagonal), to do updates in the first block ***
  corrupdate[0] = 0;
  if(corrupd>0) {
    for(j1=0;j1<npar;j1++) {
      covar2[tempi][j1][j1] = pdfsigs[j1];
    }
    corrupdate[0] = 1; //Use the matrix above and don't change it
    if(corrupd==2) corrupdate[0] = 2; //Use the matrix above and update it every ncorr iterations
  }
  
  
  
  
  
  // ********************************************************************************************************************************************************************************
  // ***  GET OFFSET STARTING VALUES  ***********************************************************************************************************************************************
  // ********************************************************************************************************************************************************************************
  
  
  //Offset starting values (only for the parameters we're fitting)
  nstart = 0;
  par2arr(&state, param);  //Put the variables in their array
  if(offsetmcmc==1) {
    printf("\n");
    for(i=0;i<npar;i++) {
      nparam[tempi][i] = param[tempi][i];  //Temporarily store the true values
    }
    logL[tempi] = -1.e30;
    while(logL[tempi] < NullLikelihood + 1.0) { //Accept only good starting values
      acceptprior[tempi] = 1;
      for(i=0;i<npar;i++) {
	if(fitpar[i]==1 && offsetpar[i]==1) {
	  param[tempi][i] = nparam[tempi][i] + offsetx * (gsl_rng_uniform(r) - 0.5) * pdfsigs[i];
	  //0:Mc, 1:eta, 2:tc, 3:logd, 4:a, 5:kappa, 6:RA, 7:sindec, 8:phi, 9:sintheta_Jo, 10: phi_Jo, 11:alpha
	  if(i==1 && (param[tempi][i]<=0.01 || param[tempi][i] > 0.25)) param[tempi][i] = max(min(gsl_rng_uniform(r)*0.25,1.0),0.01);  //Eta: 0.01<eta<0.25  \__ If it's that far outside, you may as well take a random value
	  if(i==4 && (param[tempi][i]<=1.e-5 || param[tempi][i] > 1.0)) param[tempi][i] = max(min(gsl_rng_uniform(r),1.0),1.e-5);      //Spin: 0<a<1         /   over the range of this parameter
	  if((i==5 || i==7 || i==9) && (param[tempi][i] < -2.0 || param[tempi][i] > 2.0)) param[tempi][i] = gsl_rng_uniform(r)*2.0 - 1.0;
	  if(i==5 || i==7 || i==9) param[tempi][i] = fmod(param[tempi][i]+1001.0,2.0) - 1.0;
	  if((i==6 || i==8 || i==10 || i==11) && (param[tempi][i] < -2.0*pi || param[tempi][i] > 4.0*pi)) param[tempi][i] = gsl_rng_uniform(r)*tpi;
	  acceptprior[tempi] *= prior(&param[tempi][i],i);
	}
      }
      if(acceptprior[tempi]==1) {                     //Check the value of the likelihood for this draw
	arr2par(param, &state);	                      //Get the parameters from their array
	localpar(&state, ifo, networksize);
	logL[tempi] = net_loglikelihood(&state, networksize, ifo);  //Calculate the likelihood
      }
      nstart = nstart + 1;
    }
    printf("%10s  %15s  %8s  %8s  %16s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n", "nDraws","logL-logLo","Mc","eta","tc","logdL","spin","kappa","RA","sinlati","phase","sinthJ0","phiJ0","alpha","logT");
    printf("%10d  %15.6lf  %8.5f  %8.5f  %16.6lf  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.3f\n",
	   nstart,logL[tempi]-NullLikelihood,param[tempi][0],param[tempi][1],param[tempi][2],param[tempi][3],param[tempi][4],param[tempi][5],rightAscension(param[tempi][6],GMST(param[tempi][2])),param[tempi][7],param[tempi][8],param[tempi][9],param[tempi][10],param[tempi][11],log10(temp));
  }
  
  
  // *** Set the NEW array, sigma and scale ***
  for(i=0;i<npar;i++) {
    nparam[tempi][i] = param[tempi][i];
    sig[tempi][i]   = 0.1  * pdfsigs[i];
    if(adapt==1) sig[tempi][i] = pdfsigs[i]; //Don't use adaptation
    scale[i] = 10.0 * pdfsigs[i];
    //sig[tempi][i]   = 3.0  * pdfsigs[i]; //3-sigma = delta-100%
    //scale[i] = 0.0 * pdfsigs[i]; //No adaptation
    if(i==6 || i==8 || i==10 || i==11) sig[tempi][i] = fmod(sig[tempi][i]+mtpi,tpi);  //Bring the sigma between 0 and 2pi
  }
  
  
  
  
  // ********************************************************************************************************************************************************************************
  // ***  WRITE STARTING STATE TO SCREEN AND FILE  **********************************************************************************************************************************
  // ********************************************************************************************************************************************************************************
  
  
  arr2par(param, &state);                         //Get the parameters from their array
  localpar(&state, ifo, networksize);
  logL[tempi] = net_loglikelihood(&state, networksize, ifo);  //Calculate the likelihood
  
  // *** Write to screen ***
  printf("\n\n%10s  %15s  %8s  %8s  %16s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n","cycle","logL-logLo","Mc","eta","tc","logdL","spin","kappa","RA","sinlati","phase","sinthJ0","phiJ0","alpha","logT");
  printf("%10d  %15.6lf  %8.5f  %8.5f  %16.6lf  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.3f\n",
	 0,logL[tempi]-NullLikelihood,param[tempi][0],param[tempi][1],param[tempi][2],param[tempi][3],param[tempi][4],param[tempi][5],rightAscension(param[tempi][6],GMST(param[tempi][2])),param[tempi][7],param[tempi][8],param[tempi][9],param[tempi][10],param[tempi][11],log10(temp));
  
  // *** Write to file ***
  for(tempi=0;tempi<ntemps;tempi++) {
    if(tempi==0 || saveallchains==1) {
      sprintf(outfilename,"mcmc.output.%6.6d.%2.2d",mcmcseed,tempi);
      //printf("  3:  %d  %s\n",tempi,outfilename);
      fout = fopen(outfilename,"a");
      fprintf(fout, "%12d %20.10lf  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %20.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f\n",
	      0,logL[0]-NullLikelihood,param[0][0],sig[0][0],0.0,param[0][1],sig[0][1],0.0,param[0][2],sig[0][2],0.0,param[0][3],sig[0][3],0.0,param[0][4],sig[0][4],0.0,param[0][5],sig[0][5],0.0,rightAscension(param[0][6],GMST(param[0][2])),sig[0][6],0.0,param[0][7],sig[0][7],0.0,
	      param[0][8],sig[0][8],0.0,param[0][9],sig[0][9],0.0,param[0][10],sig[0][10],0.0,param[0][11],sig[0][11],0.0);
      fclose(fout);
    }
  }
  tempi = 0;  //MUST be zero
  
  
  
  
  // ********************************************************************************************************************************************************************************
  // ***  INITIALISE PARALLEL TEMPERING  ********************************************************************************************************************************************
  // ********************************************************************************************************************************************************************************
  
  // *** Put the initial values of the parameters, sigmas etc in the different temperature chains ***
  if(ntemps>1) {
    for(tempi=1;tempi<ntemps;tempi++) {
      for(j=0;j<npar;j++) {
	param[tempi][j] = param[0][j];
	nparam[tempi][j] = nparam[0][j];
	sig[tempi][j] = sig[0][j];
	logL[tempi] = logL[0];
	nlogL[tempi] = nlogL[0];

	for(j1=0;j1<npar;j1++) {
	  for(j2=0;j2<=j1;j2++) {
	    covar2[tempi][j1][j2] = covar2[0][j1][j2];
	  }
	}
	
      }
      corrupdate[tempi] = corrupdate[0];
      //corrupdate[tempi] = 0; //Correlated update proposals only for T=1 chain?
      //corrupdate[ntemps-1] = 0; //Correlated update proposals not for hottest chain
    }
  }
  
  
  
  // ********************************************************************************************************************************************************************************
  // ***  CREATE MARKOV CHAIN   *****************************************************************************************************************************************************
  // ********************************************************************************************************************************************************************************
      
  
  iteri = 1;
  while(iteri<=iter) {  //loop over Markov-chain states 
      
    for(tempi=0;tempi<ntemps;tempi++) {  //loop over temperature ladder
      //printf(" %d  %d  %d\n",tempi,ntemps,iteri);
      
      if(partemp==1 || partemp==3) { //Chains at fixed T
	temp = temps[tempi];
      }
      if(partemp==2 || partemp==4) { //Chains with sinusoid T
	if(tempi==0) {
	  temp = 1.0;
	}
	else {
	  //temp = temps[tempi] * (1.0  +  0.5 * pow((-1.0),tempi) * sin(tpi*(double)iteri/((double)ncorr)));  //Sinusoid around the temperature T_i with amplitude 0.5*T_i and period ncorr
	  //temp = temps[tempi]  +  tempampl[tempi] * pow((-1.0),tempi) * sin(tpi*(double)iteri/((double)ncorr));  //Sinusoid around the temperature T_i with amplitude tempampl and period ncorr
	  //temp = temps[tempi]  +  tempampl[tempi] * pow((-1.0),tempi) * sin(tpi*(double)iteri/(0.5*(double)ncorr));  //Sinusoid around the temperature T_i with amplitude tempampl and period 1/2 ncorr
	  temp = temps[tempi]  +  tempampl[tempi] * pow((-1.0),tempi) * sin(tpi*(double)iteri/(5.0*(double)ncorr));  //Sinusoid around the temperature T_i with amplitude tempampl and period 5 * ncorr
	}
	//printf("%4d  %10.3lf\n",tempi,temp);
      }
      
      
      
      
      // ********************************************************************************************************************************************************************************
      // ***  UPDATE MARKOV CHAIN STATE  ************************************************************************************************************************************************
      // ********************************************************************************************************************************************************************************
      
      
      // ****************************************************************************************    Uncorrelated updates
      // ****************************************************************************************    
      if(corrupdate[tempi]<=0) {
      
	// *************************************************************    Blockwise update (10%)
	if(iteri>nburn && gsl_rng_uniform(r)<blockfrac){
	  acceptprior[tempi] = 1;
	  for(j1=0;j1<npar;j1++){
	    if(fitpar[j1]==1) {
	      nparam[tempi][j1] = param[tempi][j1] + gsl_ran_gaussian(r,sig[tempi][j1]);
	      acceptprior[tempi] *= prior(&nparam[tempi][j1],j1);
	    }
	  }
	
	  if(acceptprior[tempi]==1) {
	    arr2par(nparam, &state);	                      //Get the parameters from their array
	    localpar(&state, ifo, networksize);
	    nlogL[tempi] = net_loglikelihood(&state, networksize, ifo); //Calculate the likelihood
	    par2arr(&state, nparam);	                      //Put the variables back in their array
	  
	    if(exp(max(-30.0,min(0.0,nlogL[tempi]-logL[tempi]))) > pow(gsl_rng_uniform(r),temp) && nlogL[tempi] > NullLikelihood){  //Accept proposal if L>Lo
	      for(j1=0;j1<npar;j1++){
		if(fitpar[j1]==1) {
		  param[tempi][j1] = nparam[tempi][j1];
		  accepted[tempi][j1] += 1;
		}
	      }
	      logL[tempi] = nlogL[tempi];
	    }
	  
	  }
	} 
	// *************************************************************
	
	
	// *************************************************************    Componentwise update (90%)
	else{
	  for(j1=0;j1<npar;j1++){
	    if(fitpar[j1]==1) nparam[tempi][j1] = param[tempi][j1];
	  }
	  for(j1=0;j1<npar;j1++){
	    if(fitpar[j1]==1) {
	      nparam[tempi][j1] = param[tempi][j1] + gsl_ran_gaussian(r,sig[tempi][j1]);
	      
	      acceptprior[tempi] = prior(&nparam[tempi][j1],j1);
	      if(acceptprior[tempi]==1) {
		arr2par(nparam, &state);                          //Get the parameters from their array
		localpar(&state, ifo, networksize);
		nlogL[tempi] = net_loglikelihood(&state, networksize, ifo);   //Calculate the likelihood
		par2arr(&state, nparam);                          //Put the variables back in their array
	      
		if(exp(max(-30.0,min(0.0,nlogL[tempi]-logL[tempi]))) > pow(gsl_rng_uniform(r),temp) && nlogL[tempi] > NullLikelihood) {  //Accept proposal
		  param[tempi][j1] = nparam[tempi][j1];
		  logL[tempi] = nlogL[tempi];
		  if(adapt==1){
		    gamma = scale[j1]*pow(1.0/((double)(iteri+1)),1.0/6.0);
		    sig[tempi][j1] = max(0.0,sig[tempi][j1] + gamma*(1.0 - alphastar)); //Accept
		    if(j1==6 || j1==8 || j1==10 || j1==11) sig[tempi][j1] = fmod(sig[tempi][j1]+mtpi,tpi);  //Bring the sigma between 0 and 2pi
		  }
		  accepted[tempi][j1] += 1;
		}
		else{                                                      //Reject proposal
		  nparam[tempi][j1] = param[tempi][j1];
		  if(adapt==1){
		    gamma = scale[j1]*pow(1.0/((double)(iteri+1)),1.0/6.0);
		    sig[tempi][j1] = max(0.0,sig[tempi][j1] - gamma*alphastar); //Reject
		    if(j1==6 || j1==8 || j1==10 || j1==11) sig[tempi][j1] = fmod(sig[tempi][j1]+mtpi,tpi);  //Bring the sigma between 0 and 2pi
		    //sig[tempi][j1] = max(0.01*sig[tempi][j1], sig[tempi][j1] - gamma*alphastar);
		  }
		}
	      }
	      else{  //If new state not within boundaries
		nparam[tempi][j1] = param[tempi][j1];
		if(adapt==1) {
		  gamma = scale[j1]*pow(1.0/((double)(iteri+1)),1.0/6.0);
		  sig[tempi][j1] = max(0.0,sig[tempi][j1] - gamma*alphastar); //Reject
		  if(j1==6 || j1==8 || j1==10 || j1==11) sig[tempi][j1] = fmod(sig[tempi][j1]+mtpi,tpi);  //Bring the sigma between 0 and 2pi
		}
	      } //if(acceptprior[tempi]==1)
	    } //if(fitpar[j1]==1)
	  } //j1
	} //End componentwise update
      } //if(corrupdate[tempi]<=0): end uncorrelated updates
      // ****************************************************************************************    End uncorrelated updates
      // ****************************************************************************************    
      
      
      
      // ****************************************************************************************    Correlated block updates
      // ****************************************************************************************    
      if(corrupdate[tempi]>=1) {
      
	for(j1=0;j1<npar;j1++) {
	  work[j1] = gsl_ran_gaussian(r,1.0) * corrsig[tempi];   //Univariate gaussian random number, with sigma=1, times the adaptable sigma_correlation
	}
	acceptprior[tempi] = 1;
	for(j1=0;j1<npar;j1++){
	  if(fitpar[j1]==1) {
	    dparam = 0.0;
	    for(j2=0;j2<=j1;j2++){
	      dparam += covar2[tempi][j1][j2]*work[j2];   //Work is now a univariate gaussian random number
	    }
	    nparam[tempi][j1] = param[tempi][j1] + dparam;       //Jump from the previous parameter value
	    sig[tempi][j1] = fabs(dparam);                       //This isn't really sigma, but the proposed jump size
	    acceptprior[tempi] *= prior(&nparam[tempi][j1],j1);
	  }
	}
      
	corrscale[tempi] = corrsig[tempi];
	if(acceptprior[tempi]==1) {                                    //Then calculate the likelihood
	  arr2par(nparam, &state);	                      //Get the parameters from their array
	  localpar(&state, ifo, networksize);
	  nlogL[tempi] = net_loglikelihood(&state, networksize, ifo); //Calculate the likelihood
	  par2arr(&state, nparam);	                      //Put the variables back in their array
	  
	  if(exp(max(-30.0,min(0.0,nlogL[tempi]-logL[tempi]))) > pow(gsl_rng_uniform(r),temp) && nlogL[tempi] > NullLikelihood) {  //Accept proposal
	    for(j1=0;j1<npar;j1++){
	      if(fitpar[j1]==1) {
		param[tempi][j1] = nparam[tempi][j1];
		accepted[tempi][j1] += 1;
	      }
	    }
	    logL[tempi] = nlogL[tempi];
	    if(adapt==1){ 
	      corrsig[tempi] *= 10.0;  //Increase sigma
	    }
	  }
	  else {                                                      //Reject proposal because of low likelihood
	    if(adapt==1){ 
	      corrsig[tempi] *= 0.5;  //Decrease sigma
	    }
	  }
	}
	else {                                                      //Reject proposal because of boundary conditions.  Perhaps one should increase the step size, or at least not decrease it?
	  /*
	  if(adapt==1){ 
	    corrsig[tempi] *= 0.8;
	  }
	  */
	}

      } //if(corrupdate[tempi]>=1): End of correlated updates
      // ****************************************************************************************    End correlated block updates
      // ****************************************************************************************    
      
      
      
      // Update the dlogL = logL - logLo, and remember the parameter values where it has a maximum
      dlogL[tempi] = logL[tempi]-NullLikelihood;
      if(dlogL[tempi]>maxdlogL[tempi]) {
	maxdlogL[tempi] = dlogL[tempi];
	for(i=0;i<npar;i++) {
	  maxparam[tempi][i] = param[tempi][i];
	}
      }
      
      
      
      
      
      
      
      // ********************************************************************************************************************************************************************************
      // ***  WRITE STATE TO SCREEN AND FILE  *******************************************************************************************************************************************
      // ********************************************************************************************************************************************************************************
      
      
      //0:mc, 1:eta, 2:tc, 3:logdl, 4:spin, 5:kappa, 6: longi (->RA), 7:sinlati, 8:phase, 9:sinthJ0, 10:phiJ0, 11:alpha
      if(acceptprior[0]==1) { //Then write output and care about the correlation matrix
	
	// *** Write output to screen ***
	if(tempi==0) { //Only for the T=1 chain
	  if((iteri % (50*screenoutput))==0)  printf("\n%10s  %15s  %8s  %8s  %16s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n",
						     "cycle","logL - logLo","Mc","eta","tc","logdL","spin","kappa","RA","sinlati","phase","sinthJ0","phiJ0","alpha","logT");
	  if((iteri % screenoutput)==0)  printf("%10d  %15.6lf  %8.5f  %8.5f  %16.6lf  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.3f\n",
						iteri,dlogL[tempi],param[tempi][0],param[tempi][1],param[tempi][2],param[tempi][3],param[tempi][4],param[tempi][5],rightAscension(param[tempi][6],GMST(param[tempi][2])),param[tempi][7],param[tempi][8],param[tempi][9],param[tempi][10],param[tempi][11],log10(temp));
	}
	
	// *** Write output to file ***
	if(tempi==0 || saveallchains==1) { //For all T-chains if desired, otherwise the T=1 chain only
	  if((iteri % skip)==0){
	    sprintf(outfilename,"mcmc.output.%6.6d.%2.2d",mcmcseed,tempi);
	    fout = fopen(outfilename,"a");
	    fprintf(fout, "%12d %20.10lf  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %20.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f  %15.10lf %9.6f %6.4f\n",
		    iteri,logL[tempi]-NullLikelihood,
		    param[tempi][0],sig[tempi][0],accepted[tempi][0]/(double)iteri,param[tempi][1],sig[tempi][1],accepted[tempi][1]/(double)iteri,param[tempi][2],sig[tempi][2],accepted[tempi][2]/(double)iteri,
		    param[tempi][3],sig[tempi][3],accepted[tempi][3]/(double)iteri,param[tempi][4],sig[tempi][4],accepted[tempi][4]/(double)iteri,param[tempi][5],sig[tempi][5],accepted[tempi][5]/(double)iteri,
		    rightAscension(param[tempi][6],GMST(param[tempi][2])),sig[tempi][6],accepted[tempi][6]/(double)iteri,
		    param[tempi][7],sig[tempi][7],accepted[tempi][7]/(double)iteri,param[tempi][8],sig[tempi][8],accepted[tempi][8]/(double)iteri,param[tempi][9],sig[tempi][9],accepted[tempi][9]/(double)iteri,
		    param[tempi][10],sig[tempi][10],accepted[tempi][10]/(double)iteri,param[tempi][11],sig[tempi][11],accepted[tempi][11]/(double)iteri);
	    
	    //fflush(fout); //Make sure any 'snapshot' you take halfway is complete
	    fclose(fout);
	  }
	} //if(tempi==0)

	
	
	
	
	
      // ********************************************************************************************************************************************************************************
      // ***  CORRELATION MATRIX  *******************************************************************************************************************************************************
      // ********************************************************************************************************************************************************************************
      
      
	//if(corrupdate[tempi]==2) {  //Calculate correlations only once
	if(corrupdate[tempi]>=2) { //Calculate correlations multiple times
	
	  // *** Save state to calculate correlations ***
	  if(ihist[tempi]<ncorr) {
	    for(j1=0;j1<npar;j1++){
	      hist[tempi][j1][ihist[tempi]] = param[tempi][j1];
	    }
	    ihist[tempi] += 1;
	    sumdlogL[tempi] += dlogL[tempi];
	  }
	  
	  
	  // *** Calculate covariances/correlations ***
	  if(ihist[tempi]>=ncorr) {
	    
	    //Calculate the mean
	    for(j1=0;j1<npar;j1++){
	      mu[j1]=0.0;
	      for(i1=0;i1<ncorr;i1++){
		mu[j1]+=hist[tempi][j1][i1];
	      }
	      mu[j1]/=((double)ncorr);
	    }
	    
	    //Calculate the standard deviation. Only for printing, not used in the code
	    for(j1=0;j1<npar;j1++){
	      stdev[j1]=0.0;
	      for(i1=0;i1<ncorr;i1++){
		stdev[j1] += (hist[tempi][j1][i1]-mu[j1])*(hist[tempi][j1][i1]-mu[j1]);
	      }
	      stdev[j1] = sqrt(stdev[j1]/(double)(ncorr-1));
	    }
	    
	    //Calculate the covariances
	    for(j1=0;j1<npar;j1++){
	      for(j2=0;j2<=j1;j2++){
		covar1[tempi][j1][j2]=0.0;
		for(i1=0;i1<ncorr;i1++){
		  covar1[tempi][j1][j2] += (hist[tempi][j1][i1] - mu[j1]) * (hist[tempi][j2][i1] - mu[j2]);
		}
		covar1[tempi][j1][j2] /= (double)(ncorr-1);
	      }
	    }
	    
	    //Store the covariance matrix in a different variable, for Cholesky decomposition
	    for(j1=0;j1<npar;j1++){
	      for(j2=0;j2<=j1;j2++){
		covar2a[j1][j2] = covar1[tempi][j1][j2];
	      }
	    }
	    
	    if(tempi==0 && prmatrixinfo>0) printf("\n\n");
	    
	    //Do Cholesky decomposition
	    chol(npar,covar2a);
	    
	    //Get conditions to decide whether to accept the new matrix or not
	    acceptelems[tempi] = 0;
	    for(j1=0;j1<npar;j1++) {
	      if(fitpar[j1]==1) {
		if(covar2a[j1][j1] < covar2[tempi][j1][j1]) acceptelems[tempi] += 1; //Smaller diagonal element is better; count for how many this is the case
		if(covar2a[j1][j1]<=0.0 || isnan(covar2a[j1][j1])!=0 || isinf(covar2a[j1][j1])!=0) acceptelems[tempi] -= 9999;  //If diagonal element is <0, NaN or Inf
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
		  printf("    %10.3g",covar2[tempi][j1][j2]);
		}
		printf("\n");
	      }
	      printf("\n    New Cholesky-decomposed matrix:\n");
	      for(j1=0;j1<npar;j1++){
		for(j2=0;j2<=j1;j2++){
		  //for(j2=0;j2<npar;j2++){
		  printf("    %10.3g",covar2a[j1][j2]);
		}
		printf("\n");
	      }
	    }
	    
	    //if((corrupdate[tempi]==2 && acceptelems[tempi]>0) || (double)acceptelems[tempi] >= (double)nparfit*0.75) { //Only accept new matrix if 75% of diagonal elements are better (smaller)
	    if(acceptelems[tempi]>=0) { //Accept new matrix always, except when Infs, NaNs, 0s occur.
	      for(j1=0;j1<npar;j1++){
		for(j2=0;j2<=j1;j2++){
		  covar2[tempi][j1][j2] = covar2a[j1][j2]; 
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
		     tempi,log10(temp),sumdlogL[tempi]/(double)ihist[tempi],acceptelems[tempi],corrupdate[tempi]-2,  (double)swapTs1[tempi]/(double)iteri,(double)accepted[tempi][0]/(double)iteri);
	      for(j1=0;j1<npar;j1++) {
		tmpdbl = log10(stdev[j1]+1.e-30);
		if(tmpdbl<-9.99) tmpdbl = 0.0;
		printf("  %6.3f",tmpdbl);
	      }
	      printf("\n");
	      
	      
	      if(prpartempinfo==2 && tempi==ntemps-1) { //Print swap-rate matrix for parallel tempering
		printf("\n   Chain swap: ");
		for(j1=0;j1<ntemps-1;j1++) {
		  printf("  %6d",j1);
		}
		printf("   total\n");
		for(j1=1;j1<ntemps;j1++) {
		  printf("            %3d",j1);
		  for(j2=0;j2<ntemps-1;j2++) {
		    if(j2<j1) {
		      printf("  %6.4f",(double)swapTss[j2][j1]/(double)iteri);
		    } else {
		      printf("        ");
		    }
		  }
		  printf("  %6.4f\n",(double)swapTs2[j1]/(double)iteri);
		}
		printf("          total");
		for(j1=0;j1<ntemps-1;j1++) {
		  printf("  %6.4f",(double)swapTs1[j1]/(double)iteri);
		}
		printf("\n");
	      } //if(prpartempinfo==2 && tempi==ntemps-1)
	    } //if(partemp>=1)
	    
	    ihist[tempi] = 0;
	    sumdlogL[tempi] = 0.0;
	    
	    if( (prmatrixinfo>0 || prpartempinfo>0) && tempi==ntemps-1 ) printf("\n\n");
	  } //if(ihist[tempi]>=ncorr)
	} //if(corrupdate[tempi]>=2)
	// *************************************************************
      
      } //if(acceptprior[tempi]==1)
      
      
      
      
    } // for(tempi=0;tempi<ntemps;tempi++) {  //loop over temperature ladder
    
    
    
    
    
    
    // ********************************************************************************************************************************************************************************
    // ***  ANNEALING  ****************************************************************************************************************************************************************
    // ********************************************************************************************************************************************************************************
    
    if(partemp==0) {  //Use only when not using parallel tempering (and of course, temp0>1)
      //temp = temp0*pow(1.0/((double)(max(iteri,ncorr))),1.0/10.0);
      //temp = temp0*pow(((double)(iteri)),-0.1);
      //temp = temp0*pow(((double)(iteri)),-0.25);
      //temp = temp0 * pow(((double)(iteri)), (-log(temp0)/log((double)nburn)) );  //Temp drops from temp0 to 1.0 in nburn iterations
      //printf("%f\n",-log(temp0)/log((double)nburn));
      //temp = min( max( temp0*pow( max((double)iteri-0.1*(double)nburn,1.0) ,-log10(temp0)/log10(0.9*(double)nburn)) , 1.0) , temp0);  //Temp stays at temp0 for 10% of burnin and then drops to 1.0 at the end of the burnin (too drastic for short burnins/high percentages?)
      //temp = temp0*pow( max((double)iteri-1000.0,1.0) ,-log10(temp0)/log10((double)nburn-1000.0));  //Temp stays at temp0 for the first 1000 iterations and then drops to 1.0 at the end of the burnin (too drastic for short burnins/high percentages?)
      temp = exp(log(temp0) * ((double)(nburn-iteri)/((double)(nburn-nburn0))));
      
      //if(gsl_rng_uniform(r)<1.e-2 && iteri > nburn && dlogL[tempi] < 0.95*maxdlogL[tempi]){
      //if(iteri > nburn && dlogL[tempi] < 0.85*maxdlogL[tempi]){
      //  temp=pow(temp0,0.5);
      //  corrsig[tempi] *= 10.0;
      //}
      temp = min( max(temp,1.0) , temp0);
    }
      
    
    
    /*
    // *** Adapt temperature ladder ***
    //  This doesn't really seem to work, since getting a particular likelihood can also been done by zooming in on a piece of parameter space and dropping T; most T's go down to 1 this way.
    if(iteri>=ncorr && ihist[0]>=ncorr-1) {
      printf("\n\n");
      
      //Calculate true and expected average logL for each chain
      for(tempi=0;tempi<ntemps;tempi++) {
	avgdlogL[tempi] = sumdlogL[tempi]/(double)ncorr;
	expdlogL[tempi] = maxdlogL[0]/(double)(ntemps-1)*(double)(ntemps-tempi-1);  //The dlogL you expect for this chain if distributed flat in log(L) between L_0 and L_max
	printf("  tempi:  %d,   T: %10.4f,   maxlodL: %10.4f,   avglogL: %10.4f,   explogL: %10.4f\n",tempi,temps[tempi],maxdlogL[tempi],avgdlogL[tempi],expdlogL[tempi]);
      }
      
      printf("\n");
      //See between which two average logL's each expected logL falls
      temp1 = 0;
      temp2 = ntemps-1;
      for(tempi=1;tempi<ntemps-1;tempi++) {
	for(t=0;t<ntemps-1;t++) {
	  if(avgdlogL[t]>=expdlogL[tempi]) temp1 = t;
	}
	for(t=ntemps;t>0;t--) {
	  if(avgdlogL[t]<expdlogL[tempi]) temp2 = t;
	}
	
	//Interpolate between temp1 and temp2
	tmpdbl = temps[temp1] + (temps[temp2]-temps[temp1])/(avgdlogL[temp2]-avgdlogL[temp1])*(expdlogL[tempi]-avgdlogL[temp1]);
	printf("  t: %d,  t1-2: %d %d,    avglogL:%10.4f,   t2-t1:%10.4f,  avglogL1-avglogL2:%10.4f tmpdbl:%10.4f\n",tempi,temp1,temp2,avgdlogL[tempi], temps[temp2]-temps[temp1],  avgdlogL[temp1]-avgdlogL[temp2], tmpdbl );
	
	newtemps[tempi] = tmpdbl;
      } //for(tempi=1;tempi<ntemps-1;tempi++)

      for(tempi=1;tempi<ntemps-1;tempi++) {
	temps[tempi] = newtemps[tempi];
      } //for(tempi=1;tempi<ntemps-1;tempi++)
      
      printf("\n\n");
    }
    */
    
    
    
    
    
    
    
    
    
    // ********************************************************************************************************************************************************************************
    // ***  PARALLEL TEMPERING  *******************************************************************************************************************************************************
    // ********************************************************************************************************************************************************************************
    
    // *** Swap states between T-chains ***
    if(acceptprior[0]==1) {
      if(partemp>=1 && ntemps>1) {
	/*
	//Swap parameters and likelihood between two adjacent chains
	for(tempi=1;tempi<ntemps;tempi++) {
	  if(exp(max(-30.0,min(0.0,(1.0/temps[tempi-1]-1.0/temps[tempi])*(logL[tempi]-logL[tempi-1])))) > gsl_rng_uniform(r)) { //Then swap...
	    for(i=0;i<npar;i++) {
	      tmpdbl = param[tempi][i]; //Temp var
	      param[tempi][i] = param[tempi-1][i];
	      param[tempi-1][i] = tmpdbl;
	    }
	    tmpdbl = logL[tempi];
	    logL[tempi] = logL[tempi-1];
	    logL[tempi-1] = tmpdbl;
	    swapTs1[tempi-1] += 1;
	  }
	} //tempi
	*/
	
	//Swap parameters and likelihood between any two chains
	for(tempi=0;tempi<ntemps-1;tempi++) {
	  for(tempj=tempi+1;tempj<ntemps;tempj++) {
	    
	    if(exp(max(-30.0,min(0.0,(1.0/temps[tempi]-1.0/temps[tempj])*(logL[tempj]-logL[tempi])))) > gsl_rng_uniform(r)) { //Then swap...
	      for(i=0;i<npar;i++) {
		tmpdbl = param[tempj][i]; //Temp var
		param[tempj][i] = param[tempi][i];
		param[tempi][i] = tmpdbl;
	      }
	      tmpdbl = logL[tempj];
	      logL[tempj] = logL[tempi];
	      logL[tempi] = tmpdbl;
	      swapTss[tempi][tempj] += 1;
	      swapTs1[tempi] += 1;
	      swapTs2[tempj] += 1;
	    }
	  } //tempj
	} //tempi
	
      } //if(partemp>=1 && ntemps>1)
      iteri++;
    } //if(acceptprior[0]==1)
    
  } // while(iter<=niter) {  //loop over markov chain states 
  
  
  
  
  // ********************************************************************************************************************************************************************************
  // ***  FREE MEMORY  **************************************************************************************************************************************************************
  // ********************************************************************************************************************************************************************************
  
  //fclose(fout);
  printf("\n");
  
  gsl_rng_free(r);
  
  free(corrupdate);
  free(acceptelems);
  
  free(temps);
  free(newtemps);
  free(tempampl);
  free(logL);
  free(nlogL);
  free(dlogL);
  free(maxdlogL);
  free(sumdlogL);
  free(avgdlogL);
  free(expdlogL);
  
  free(corrsig);
  free(corrscale);
  free(acceptprior);
  free(ihist);
  free(swapTs1);
  free(swapTs2);
  
  for(i=0;i<ntemps;i++) {
    free(accepted[i]);
    free(swapTss[i]);
    free(param[i]);
    free(nparam[i]);
    free(sig[i]);
    free(maxparam[i]);
  }
  free(accepted);
  free(swapTss);
  free(param);
  free(nparam);
  free(sig);
  free(maxparam);
  
  for(i=0;i<npar;i++) {
    free(covar2a[i]);
  }
  free(covar2a);
  
  for(i=0;i<ntemps;i++) {
    for(j=0;j<npar;j++) {
      free(hist[i][j]);
      free(covar1[i][j]);
      free(covar2[i][j]);
    }
    free(hist[i]);
    free(covar1[i]);
    free(covar2[i]);
  }
  free(hist);
  free(covar1);
  free(covar2);
  
  free(state.loctc);
  free(state.localti);
  free(state.locazi);
  free(state.locpolar); 
  
  
}
// End adaptive MCMC





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



int prior(double *par, int p)
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
  
  dt = 0.05; //This should be dt/2
  lb[2] = prior_tc_mean-dt; //t_c
  ub[2] = prior_tc_mean+dt;
  
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





