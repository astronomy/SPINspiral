// mcmc_signal.c:
// Routines to calculate likelihood, SNR, match, etc.


#include <mcmc.h>




double match(struct parset *par, struct interferometer *ifo[], int i, int networksize)
//Calculate the match between two waveforms
{
  if(MvdSdebug) printf("    Match %d\n",i);
  double match = 0.0;
  //double rate = (double) ifo[i]->samplerate; //Double instead of int, not needed here
  int    j;
  fftw_complex *FTout1,*FTout2;                  /* FT output (type here identical to `(double) complex')    */
  FTout1 = fftw_malloc(sizeof(fftw_complex) * (ifo[i]->FTsize));
  FTout2 = fftw_malloc(sizeof(fftw_complex) * (ifo[i]->FTsize));
  struct parset truepar;
  double m1m2=0.0,m1m1=0.0,m2m2=0.0;
  
  /* fill `ifo[i]->FTin' with time-domain template: */
  template(par, ifo, i);  
  /* window template:                               */
  for (j=0; j<ifo[i]->samplesize; ++j) ifo[i]->FTin[j] *= ifo[i]->FTwindow[j];
  /* execute Fourier transform of signal template:  */
  fftw_execute(ifo[i]->FTplan);
  for (j=ifo[i]->lowIndex; j<=ifo[i]->highIndex; ++j) FTout1[j] = ifo[i]->FTout[j];
  
  //Get the true parameters
  settrueparameters(&truepar);
  truepar.loctc    = (double*)calloc(networksize,sizeof(double));
  truepar.localti  = (double*)calloc(networksize,sizeof(double));
  truepar.locazi   = (double*)calloc(networksize,sizeof(double));
  truepar.locpolar = (double*)calloc(networksize,sizeof(double));
  localpar(&truepar, ifo, networksize);
  template(&truepar, ifo, i);  
  
  for (j=0; j<ifo[i]->samplesize; ++j) ifo[i]->FTin[j] *= ifo[i]->FTwindow[j];
  /* execute Fourier transform of signal template:  */
  fftw_execute(ifo[i]->FTplan);
  for (j=ifo[i]->lowIndex; j<=ifo[i]->highIndex; ++j) FTout2[j] = ifo[i]->FTout[j];
  
  
  /* sum over the Fourier frequencies within operational range (some 40-1500 Hz or similar): */
  for (j=ifo[i]->lowIndex; j<=ifo[i]->highIndex; ++j){
    m1m2 += (double)( (conj(FTout1[j])*FTout2[j] + FTout1[j]*conj(FTout2[j])) / (exp(ifo[i]->noisePSD[j-ifo[i]->lowIndex])) ); // /rate is cancelled out
    m1m1 += (double)( (conj(FTout1[j])*FTout1[j] + FTout1[j]*conj(FTout1[j])) / (exp(ifo[i]->noisePSD[j-ifo[i]->lowIndex])) );
    m2m2 += (double)( (conj(FTout2[j])*FTout2[j] + FTout2[j]*conj(FTout2[j])) / (exp(ifo[i]->noisePSD[j-ifo[i]->lowIndex])) );
  }
  fftw_free(FTout1);
  fftw_free(FTout2);
  //  printf("%g %g %g\n",m1m2,m1m1,m2m2);
  match = m1m2/sqrt(m1m1*m2m2);
  return match;
}


double net_loglikelihood(struct parset *par, int networksize, struct interferometer *ifo[])
//Calculate the loglikelihood for a network of IFOs
{
  if(MvdSdebug) printf("  Net_loglikelihood\n");
  double result = 0.0;
  int i;
  for (i=0; i<networksize; ++i){
    result += ifo_loglikelihood(par, ifo, i);
  }
  return result;
}



double ifo_loglikelihood(struct parset *par, struct interferometer *ifo[], int i)
//Calculate the loglikelihood for a (1) given IFO
{
  if(MvdSdebug) printf("    Ifo_loglikelihood %d\n",i);
  double absdiff, result = 0.0;
  double rate = (double) ifo[i]->samplerate; /* double instead of int*/
  int    j;

  /* fill `ifo[i]->FTin' with time-domain template: */
  template(par, ifo, i);
  
  /* window template:                               */
  for (j=0; j<ifo[i]->samplesize; ++j) ifo[i]->FTin[j] *= ifo[i]->FTwindow[j];
  /* execute Fourier transform of signal template:  */
  fftw_execute(ifo[i]->FTplan);
  /* sum over the Fourier frequencies within operational range (some 40-1500 Hz or similar): */
  for (j=ifo[i]->lowIndex; j<=ifo[i]->highIndex; ++j){
    absdiff = cabs(ifo[i]->raw_dataTrafo[j] - (ifo[i]->FTout[j] / rate));
    /* squared (absolute) difference divided by noise PSD(f): */
    result += exp(2.0*log(absdiff) - ifo[i]->noisePSD[j-ifo[i]->lowIndex]);
  }
  result *= -2.0/ifo[i]->deltaFT; /* divide by (FT'd) data segment length */
  return result;
}


double signaltonoiseratio(struct parset *par, struct interferometer *ifo[], int i)
/* SNR of signal corresponding to parameter set, w.r.t. i-th interferometer's noise. */
/* (see SNR definition in Christensen/Meyer/Libson(2004), page 323)                  */
{
  if(MvdSdebug) printf("Signaltonoiseratio %d\n",i);
  double signalabs, result = 0.0;
  double rate = (double) ifo[i]->samplerate; /* double instead of int*/
  double lograte = log(rate);
  int    j;
  /* fill `ifo[i]->FTin' with time-domain template: */
  template(par, ifo, i);
  
  /* window template:                               */
  for (j=0; j<ifo[i]->samplesize; ++j) ifo[i]->FTin[j] *= ifo[i]->FTwindow[j];
  /* execute Fourier transform of signal template:  */
  fftw_execute(ifo[i]->FTplan);
  /* sum over the Fourier frequencies within operational range (some 40-1500 Hz or similar): */
  for (j=ifo[i]->lowIndex; j<=ifo[i]->highIndex; ++j){
    signalabs = cabs(ifo[i]->FTout[j]);
    /* squared absolute signal divided by noise PSD(f): */
    result += exp(2.0*(log(signalabs)-lograte) - ifo[i]->noisePSD[j-ifo[i]->lowIndex]);
  }
  result /= ifo[i]->deltaFT; /* divide by (FT'd) data segment length */
  result = 2.0*sqrt(result);
  return result;
}


void writesignaltodisc(struct parset *par, struct interferometer *ifo[], int i)
// Write data (signal, noise, and PSDs) to disc
{
  if(MvdSdebug) printf("Writesignaltodisc %d\n",i);
  double signalabs, result = 0.0;
  double rate = (double) ifo[i]->samplerate; /* double instead of int*/
  double lograte = log(rate);
  int    j;
  
  //  printf(" %g %g %g\n",rate,(double)ifo[i]->FTsize,ifo[i]->deltaFT);
  
  /* fill `ifo[i]->FTin' with time-domain template: */
  template(par, ifo, i);
  
  //Write signal in time domain
  char filename[1000]="";
  sprintf(filename, "%s-signal.dat", ifo[i]->name);  //Write in current dir
  FILE *dump = fopen(filename,"w");
  fprintf(dump,"%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n","m1","m2","mc","eta","tc","dl","lat","lon","phase","spin","kappa","thJ0","phJ0","alpha");
  fprintf(dump,"%12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	  par->m1,par->m2,par->mc,par->eta,par->tc,exp(par->logdl),asin(par->sinlati)*r2d,par->longi*r2d,par->phase,par->spin,par->kappa,par->sinthJ0,par->phiJ0,par->alpha);
  fprintf(dump,"t Ht\n");
  for (j=0; j<ifo[i]->samplesize; ++j)
    fprintf(dump, "%9.9f %.6e\n", ifo[i]->FTstart+(((double)j)/rate), 
	    ifo[i]->FTin[j]);
  fclose(dump); if(intscrout) printf(" : (signal written to file)\n");
  
  
  double f=0.0;
  
  //Write noise PSD
  if(1==2) {
    sprintf(filename, "%s-noisePSD.dat", ifo[i]->name);  //Write in current dir
    FILE *dump1 = fopen(filename,"w");
    
    fprintf(dump1,"%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n","m1","m2","mc","eta","tc","dl","lat","lon","phase","spin","kappa","thJ0","phJ0","alpha");
    fprintf(dump1,"%12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	    par->m1,par->m2,par->mc,par->eta,par->tc,exp(par->logdl),asin(par->sinlati)*r2d,par->longi*r2d,par->phase,par->spin,par->kappa,par->sinthJ0,par->phiJ0,par->alpha);
    fprintf(dump1,"f Nf\n");
    /* loop over the Fourier frequencies within operational range (some 40-1500 Hz or similar): */
    f=0.0;
    for (j=0; j<ifo[i]->indexRange; ++j){
      f = (((double)(j+ifo[i]->lowIndex))/((double)ifo[i]->FTsize*2.0)) * rate;
      fprintf(dump1, "%9.9f %.6e\n",log10(f), log10(2.0*sqrt(exp(ifo[i]->noisePSD[j])/ifo[i]->deltaFT))  );
    }
    fclose(dump1); if(intscrout) printf(" : (noise PSD written to file)\n");
  }
  
  /* window template:                               */
  for (j=0; j<ifo[i]->samplesize; ++j) ifo[i]->FTin[j] *= ifo[i]->FTwindow[j];
  /* execute Fourier transform of signal template:  */
  fftw_execute(ifo[i]->FTplan);
  /* sum over the Fourier frequencies within operational range (some 40-1500 Hz or similar): */
  for (j=ifo[i]->lowIndex; j<=ifo[i]->highIndex; ++j){
    signalabs = cabs(ifo[i]->FTout[j]);
    /* squared absolute signal divided by noise PSD(f): */
    result += exp(2.0*(log(signalabs)-lograte) - ifo[i]->noisePSD[j-ifo[i]->lowIndex]);
  }
  result /= ifo[i]->deltaFT; /* divide by (FT'd) data segment length */
  result = 2.0*sqrt(result);
  
  
  //Write signal PSD
  sprintf(filename, "%s-signalPSD.dat", ifo[i]->name);  //Write in current dir
  FILE *dump2 = fopen(filename,"w");
  fprintf(dump2,"%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n","m1","m2","mc","eta","tc","dl","lat","lon","phase","spin","kappa","thJ0","phJ0","alpha");
  fprintf(dump2,"%12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	  par->m1,par->m2,par->mc,par->eta,par->tc,exp(par->logdl),asin(par->sinlati)*r2d,par->longi*r2d,par->phase,par->spin,par->kappa,par->sinthJ0,par->phiJ0,par->alpha);
  fprintf(dump2,"f Nf\n");
  /* loop over the Fourier frequencies within operational range (some 40-1500 Hz or similar): */
  for (j=ifo[i]->lowIndex; j<=ifo[i]->highIndex; ++j){
    f = (((double)j)/((double)ifo[i]->FTsize*2.0)) * rate;
    fprintf(dump2, "%9.9f %.6e\n",log10(f), log10(2.0*sqrt(exp(2.0*(log(cabs(ifo[i]->FTout[j]))-lograte))/ifo[i]->deltaFT))  );
  }
  fclose(dump2); if(intscrout) printf(" : (signal PSD written to file)\n");
  return;
}

