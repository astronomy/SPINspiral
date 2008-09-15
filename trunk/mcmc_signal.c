// mcmc_signal.c:
// Routines to calculate likelihood, SNR, match, etc.


#include <mcmc.h>



double net_loglikelihood(struct parset *par, int networksize, struct interferometer *ifo[])
//Calculate the loglikelihood for a network of IFOs
{
  //if(MvdSdebug) printf("  Net_loglikelihood\n");
  double result = 0.0;
  int i;
  for (i=0; i<networksize; ++i){
    result += ifo_loglikelihood(par, ifo, i);
  }
  return result;
}



double ifo_loglikelihood(struct parset *par, struct interferometer *ifo[], 
	int ifonr)
//Calculate the loglikelihood for a (1) given IFO
{
  int j=0;

  // Fill `ifo[ifonr]->FTin' with time-domain template:
  template(par, ifo, ifonr);

  // Window template:
  for(j=0; j<ifo[ifonr]->samplesize; ++j) 
	ifo[ifonr]->FTin[j] *= ifo[ifonr]->FTwindow[j];

  // Execute Fourier transform of signal template:
  fftw_execute(ifo[ifonr]->FTplan);

  // Compute the overlap between waveform and data:
  double overlaphd = vecoverlap(ifo[ifonr]->raw_dataTrafo, 
	ifo[ifonr]->FTout, ifo[ifonr]->noisePSD,
        ifo[ifonr]->lowIndex, ifo[ifonr]->highIndex, ifo[ifonr]->deltaFT);
  //correct FFT for sampling rate of waveform
  overlaphd/=((double)ifo[ifonr]->samplerate);  

  // Compute the overlap between waveform and itself:
  double overlaphh = vecoverlap(ifo[ifonr]->FTout,
        ifo[ifonr]->FTout, ifo[ifonr]->noisePSD,
        ifo[ifonr]->lowIndex, ifo[ifonr]->highIndex, ifo[ifonr]->deltaFT);
  //correct FFT for sampling rate of waveform
  overlaphh/=((double)ifo[ifonr]->samplerate);
  overlaphh/=((double)ifo[ifonr]->samplerate);

  return (overlaphd-0.5*overlaphh);

/*
  //Alternative: about 8% slower because of extra copies of FFT output
  //   and division of each FFT point by ((double)ifo[ifonr]->samplerate)
  fftw_complex *FFTwaveform =
        fftw_malloc(sizeof(fftw_complex) * (ifo[ifonr]->FTsize));
  signalFFT(FFTwaveform, par, ifo, ifonr);                             
  // Compute the overlap between waveform and data:
  double overlaphd = vecoverlap(ifo[ifonr]->raw_dataTrafo,
        FFTwaveform, ifo[ifonr]->noisePSD,
        ifo[ifonr]->lowIndex, ifo[ifonr]->highIndex, ifo[ifonr]->deltaFT);
  // Compute the overlap between waveform and itself:
  double overlaphh = vecoverlap(
       FFTwaveform, FFTwaveform, ifo[ifonr]->noisePSD,
        ifo[ifonr]->lowIndex, ifo[ifonr]->highIndex, ifo[ifonr]->deltaFT);
  return (overlaphd-0.5*overlaphh);
*/
  
/*
  // Clean, but computes waveform thrice for a slowdown by ~x3
  return (overlapwithdata(par, ifo, ifonr) 
	- 0.5*paroverlap(par, par, ifo, ifonr));
*/ 
}



double signaltonoiseratio(struct parset *par, struct interferometer *ifo[], int ifonr)
// SNR of signal corresponding to parameter set, w.r.t. i-th interferometer's noise.
// (see SNR definition in Christensen/Meyer/Libson (2004), p.323)
{
  int j=0;

  // Fill `ifo[ifonr]->FTin' with time-domain template:
  template(par, ifo, ifonr);

  // Window template:
  for(j=0; j<ifo[ifonr]->samplesize; ++j)
        ifo[ifonr]->FTin[j] *= ifo[ifonr]->FTwindow[j];

  // Execute Fourier transform of signal template:
  fftw_execute(ifo[ifonr]->FTplan);

  // Compute the overlap between waveform and itself:
  double overlaphh = vecoverlap(ifo[ifonr]->FTout,
        ifo[ifonr]->FTout, ifo[ifonr]->noisePSD,
        ifo[ifonr]->lowIndex, ifo[ifonr]->highIndex, ifo[ifonr]->deltaFT);
  //correct FFT for sampling rate of waveform
  overlaphh/=((double)ifo[ifonr]->samplerate);
  overlaphh/=((double)ifo[ifonr]->samplerate);

  return sqrt(overlaphh);

/*
  //Clean, but recomputes waveform multiple times
  return sqrt(paroverlap(par, par, ifo, ifonr));
*/
}



void writesignaltodisc(struct parset *par, struct interferometer *ifo[], int i)
// Write data (signal, noise, and PSDs) to disc
{
  if(MvdSdebug) printf("Writesignaltodisc %d\n",i);
  double rate = (double)ifo[i]->samplerate; // Double instead of int
  double lograte = log(rate);
  double f=0.0;
  int j=0;
  
  //printf(" %g %g %g\n",rate,(double)ifo[i]->FTsize,ifo[i]->deltaFT);
  
  
  // Fill `ifo[i]->FTin' with time-domain template:
  template(par, ifo, i);
  
  
  // Write signal in time domain:
  char filename[1000]="";
  sprintf(filename, "%s-signal.dat", ifo[i]->name);  // Write in current dir
  FILE *dump = fopen(filename,"w");
  fprintf(dump,"%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n","m1","m2","mc","eta","tc","dl","lat","lon","phase","spin","kappa","thJ0","phJ0","alpha");
  fprintf(dump,"%12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	  par->m1,par->m2,par->mc,par->eta,par->tc,exp(par->logdl),asin(par->sinlati)*r2d,par->longi*r2d,par->phase,par->spin,par->kappa,par->sinthJ0,par->phiJ0,par->alpha);
  
  fprintf(dump,"       GPS time (s)         H(t)\n");
  for(j=0; j<ifo[i]->samplesize; ++j)  fprintf(dump, "%9.9f %.6e\n", ifo[i]->FTstart+(((double)j)/rate),ifo[i]->FTin[j]);
  fclose(dump);
  if(intscrout) printf(" : (signal written to file)\n");
  
  
  
  
  
  //Write noise ASD  (Square root of the estimated noise PSD (i.e., no injected signal))
  sprintf(filename, "%s-noiseASD.dat", ifo[i]->name);  // Write in current dir
  FILE *dump1 = fopen(filename,"w");
  
  fprintf(dump1,"%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n","m1","m2","mc","eta","tc","dl","lat","lon","phase","spin","kappa","thJ0","phJ0","alpha");
  fprintf(dump1,"%12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	  par->m1,par->m2,par->mc,par->eta,par->tc,exp(par->logdl),asin(par->sinlati)*r2d,par->longi*r2d,par->phase,par->spin,par->kappa,par->sinthJ0,par->phiJ0,par->alpha);
  fprintf(dump1,"       f (Hz)          H(f)\n");
  
  // Loop over the Fourier frequencies within operational range (some 40-1500 Hz or similar):
  f=0.0;
  double fact1a = rate/(2.0 * (double)ifo[i]->FTsize);
  double fact1b = sqrt(2.0)*2.0/sqrt(ifo[i]->deltaFT);  //Extra factor of sqrt(2) to get the numbers right with the outside world
  for(j=0; j<ifo[i]->indexRange; ++j){
    f = fact1a * (double)(j+ifo[i]->lowIndex);
    fprintf(dump1, "%13.6f %13.6e %13.6e\n",f, fact1b * sqrt(ifo[i]->noisePSD[j]), 0.0 ); // Add 0 to create output compatible with FFT
  }
  fclose(dump1);
  if(intscrout) printf(" : (noise ASD written to file)\n");
  
  
  
  
  // Write signal FFT to disc (i.e., amplitude spectrum of signal w/o noise):
  sprintf(filename, "%s-signalFFT.dat", ifo[i]->name);  // Write in current dir
  FILE *dump2 = fopen(filename,"w");
  fprintf(dump2,"%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n","m1","m2","mc","eta","tc","dl","lat","lon","phase","spin","kappa","thJ0","phJ0","alpha");
  fprintf(dump2,"%12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	  par->m1,par->m2,par->mc,par->eta,par->tc,exp(par->logdl),asin(par->sinlati)*r2d,par->longi*r2d,par->phase,par->spin,par->kappa,par->sinthJ0,par->phiJ0,par->alpha);
  fprintf(dump2,"       f (Hz)    real(H(f))    imag(H(f))\n");
  
  double fact2a = rate / (2.0*(double)ifo[i]->FTsize);
  double fact2b = sqrt(2.0)*2.0/(rate*sqrt(ifo[i]->deltaFT));  //Extra factor of sqrt(2) to get the numbers right with the outside world
  // Loop over the Fourier frequencies within operational range (some 40-1500 Hz or similar):
  for(j=ifo[i]->lowIndex; j<=ifo[i]->highIndex; ++j){
    f = fact2a * (double)j;
    double complex tempvar = fact2b * ifo[i]->FTout[j];
    //fprintf(dump2, "%13.6e %13.6e\n",log10(f), log10(sqrt(0.5)*cabs(tempvar))  ); //Produces the old output
    fprintf(dump2, "%13.6e %13.6e %13.6e\n",f, creal(tempvar), cimag(tempvar) );  //Save the real and imaginary parts of the signal FFT
  }
  fclose(dump2); 
  if(intscrout) printf(" : (signal FFT written to file)\n");
  return;
}



double match(struct parset *par, struct interferometer *ifo[], int i, int networksize)
//Calculate the match between two waveforms
{
  if(MvdSdebug) printf("    Match %d\n",i);
  double match = 0.0;
  int j=0;
  fftw_complex *FTout1,*FTout2;                  // FT output (type here identical to `(double) complex')
  FTout1 = fftw_malloc(sizeof(fftw_complex) * (ifo[i]->FTsize));
  FTout2 = fftw_malloc(sizeof(fftw_complex) * (ifo[i]->FTsize));
  struct parset truepar;
  double m1m2=0.0,m1m1=0.0,m2m2=0.0;
  
  // Fill `ifo[i]->FTin' with time-domain template:
  template(par, ifo, i); 
  
  // Window template:
  for(j=0; j<ifo[i]->samplesize; ++j) ifo[i]->FTin[j] *= ifo[i]->FTwindow[j];
  
  // Execute Fourier transform of signal template:
  fftw_execute(ifo[i]->FTplan);
  for(j=ifo[i]->lowIndex; j<=ifo[i]->highIndex; ++j) FTout1[j] = ifo[i]->FTout[j];
  
  //Get the true parameters and the corresponding waveform template:
  gettrueparameters(&truepar);
  truepar.loctc    = (double*)calloc(networksize,sizeof(double));
  truepar.localti  = (double*)calloc(networksize,sizeof(double));
  truepar.locazi   = (double*)calloc(networksize,sizeof(double));
  truepar.locpolar = (double*)calloc(networksize,sizeof(double));
  localpar(&truepar, ifo, networksize);
  template(&truepar, ifo, i);
  
  
  // Window template:
  for(j=0; j<ifo[i]->samplesize; ++j)  ifo[i]->FTin[j] *= ifo[i]->FTwindow[j];
  
  // Execute Fourier transform of signal template:
  fftw_execute(ifo[i]->FTplan);
  for(j=ifo[i]->lowIndex; j<=ifo[i]->highIndex; ++j) FTout2[j] = ifo[i]->FTout[j];
  
  // Sum over the Fourier frequencies within operational range (some 40-1500 Hz or similar):
  for (j=ifo[i]->lowIndex; j<=ifo[i]->highIndex; ++j){
    m1m2 += (double)( (conj(FTout1[j])*FTout2[j] + FTout1[j]*conj(FTout2[j])) / (ifo[i]->noisePSD[j-ifo[i]->lowIndex]) ); // 1/rate cancels out
    m1m1 += (double)( (conj(FTout1[j])*FTout1[j] + FTout1[j]*conj(FTout1[j])) / (ifo[i]->noisePSD[j-ifo[i]->lowIndex]) );
    m2m2 += (double)( (conj(FTout2[j])*FTout2[j] + FTout2[j]*conj(FTout2[j])) / (ifo[i]->noisePSD[j-ifo[i]->lowIndex]) );
  }
  fftw_free(FTout1);
  fftw_free(FTout2);
  //  printf("%g %g %g\n",m1m2,m1m1,m2m2);
  match = m1m2/sqrt(m1m1*m2m2);
  return match;
}




double parmatch(struct parset * par1,struct parset * par2, struct interferometer *ifo[], int networksize)
// Compute match between waveforms with parameter sets par1 and par2
{
  double match=0.0,ovrlp11=0.0,ovrlp12=0.0,ovrlp22=0.0;
  int ifonr=0;
  
  for(ifonr=0; ifonr<networksize; ifonr++){
    ovrlp11 = paroverlap(par1,par1,ifo,ifonr);
    ovrlp22 = paroverlap(par2,par2,ifo,ifonr);
    ovrlp12 = paroverlap(par1,par2,ifo,ifonr);
    match += ovrlp12/sqrt(ovrlp11*ovrlp22);
  }
  match /= (double)networksize;
  return match;
}


double overlapwithdata(struct parset *par, struct interferometer *ifo[], int ifonr)
//compute frequency domain overlap of waveform of given parameters with raw data
{
  int j=0;

  fftw_complex *FFTwaveform = 
	fftw_malloc(sizeof(fftw_complex) * (ifo[ifonr]->FTsize));
  signalFFT(FFTwaveform, par, ifo, ifonr);

  double overlap=
     vecoverlap(ifo[ifonr]->raw_dataTrafo, FFTwaveform, ifo[ifonr]->noisePSD, 
	ifo[ifonr]->lowIndex, ifo[ifonr]->highIndex, ifo[ifonr]->deltaFT);

  free(FFTwaveform);
  return overlap;
}



double paroverlap(struct parset * par1, struct parset * par2, struct interferometer *ifo[], int ifonr)
//Compute the overlap in the frequency domain between two waveforms with parameter sets par1 and par2
{
  double overlap = 0.0;
  int j=0;
  fftw_complex *FFT1 = fftw_malloc(sizeof(fftw_complex) * (ifo[ifonr]->FTsize));
  fftw_complex *FFT2 = fftw_malloc(sizeof(fftw_complex) * (ifo[ifonr]->FTsize));
  
  // Get waveforms, FFT them and store them in FFT1,2
  signalFFT(FFT1, par1, ifo, ifonr);
  signalFFT(FFT2, par2, ifo, ifonr);
  
  // Compute the overlap between the vectors FFT1,2, between index i1 and i2:
  overlap = vecoverlap(FFT1, FFT2, ifo[ifonr]->noisePSD, 
	ifo[ifonr]->lowIndex, ifo[ifonr]->highIndex, ifo[ifonr]->deltaFT);
  
  fftw_free(FFT1);
  fftw_free(FFT2);
  
  return overlap;
}



double vecoverlap(fftw_complex *vec1, fftw_complex *vec2, double * noise, int j1, int j2, double deltaFT)
//Compute the overlap of vectors vec1 and vec2, between indices j1 and j2
{
  int j=0;
  double overlap = 0.0;
  
  // Sum over the Fourier frequencies within operational range (some 40-1500 Hz or similar):
  for (j=j1; j<=j2; ++j){
    overlap += 4.0*creal( vec1[j]*conj(vec2[j]) / noise[j-j1] ) / deltaFT;   
  }
  return overlap;
}



void signalFFT(fftw_complex * FFTout, struct parset *par, 
	struct interferometer *ifo[], int ifonr)
//Compute the FFT of a waveform with given parameter set
{
  int j=0;
  
  if(FFTout==NULL)
  {
	printf("Memory should be allocated for FFTout vector");
	printf(" before call to signalFFT()\n");
	exit(1);
  }

  // Fill `ifo[i]->FTin' with time-domain template:
  template(par, ifo, ifonr);

  // Window template:
  for(j=0; j<ifo[ifonr]->samplesize; ++j) 
	ifo[ifonr]->FTin[j] *= ifo[ifonr]->FTwindow[j];

  // Execute Fourier transform of signal template:
  fftw_execute(ifo[ifonr]->FTplan);

  for(j=0; j<ifo[ifonr]->FTsize; j++)
	FFTout[j]=ifo[ifonr]->FTout[j]/((double)ifo[ifonr]->samplerate);

/*
  //Allocate the parset vectors, compute the local parameters and the time-domain template:
  localpar(&par, ifo, networksize);
  template(&par, ifo, ifonr); 
  
  // Window the time-domain template:
  for(j=0; j<ifo[ifonr]->samplesize; ++j) ifo[ifonr]->FTin[j] *= ifo[ifonr]->FTwindow[j];
  
  // Execute Fourier transform of signal template:
  fftw_execute(ifo[ifonr]->FTplan);
  for(j=ifo[ifonr]->lowIndex; j<=ifo[ifonr]->highIndex; ++j) FFT[j] = ifo[ifonr]->FTout[j];
*/
}




void computeFishermatrixIFO(struct parset *par, int npar, struct interferometer *ifo[], int networksize, int ifonr, double **matrix)
// Compute  Fisher matrix for parameter set par for a given IFO
{
/*  int ip=0,jp=0,j=0,j1=0,j2=0,nFT=0;
  struct parset par1;
  allocparset(&par1, networksize);
  double pars[npar];
  
  nFT = ifo[ifonr]->FTsize;
  double dx[npar], noise[nFT];
  double _Complex FFT0[nFT], FFT1[nFT], dFFTs[npar][nFT];
  
  j1 = ifo[ifonr]->lowIndex;
  j2 = ifo[ifonr]->highIndex;
  
  // Compute the FFTed signal for the default parameter set FFT0
  signalFFT(FFT0, par, ifo, ifonr);
  
  for(ip=0;ip<npar;ip++) {
    // Change parameter ip with dx
    dx[ip] = 1.e-5;       // Should this be the same for each parameter?
    par2arr(par,pars);    // Stick the default parameter set par into the array pars
    pars[ip] += dx[ip];   // Change parameter ip with dx
    arr2par(pars,&par1);  // Put the changed parameter set into struct par1
    
    // Compute the FFTed signal for this parameter set FFT1
    signalFFT(par1, ifo, networksize, ifonr, FFT1);
    
    // Compute the partial derivative to parameter ip
    for(j=j1;j<=j2;j++) {
      dFFTs[ip][j] = (FFT1[j]-FFT0[j])/dx[ip];
    }
  }
  
  
  // Compute the actual Fisher matrix (diagonal + lower triangle)
  for(ip=0;ip<npar;ip++) {
    for(jp=0;jp<=ip;jp++) {
      matrix[ip][jp] = vecoverlap(dFFTs[ip], dFFTs[jp], 
	ifo[ifonr]->noisePSD, j1, j2, ifo[ifonr]->deltaFT);
    }
  }
  
  // Copy the lower to the upper triangle to get a complete matrix
  for(ip=0;ip<npar;ip++) {
    for(jp=ip;jp<npar;jp++) {
      matrix[ip][jp] = matrix[jp][ip];
    }
  }
  
  freeparset(&par1);
*/
}


void computeFishermatrix(struct parset *par, int npar, struct interferometer *ifo[], int networksize, double **matrix)
// Compute the Fisher matrix for a network of IFOs, using computeFishermatrixIFO to compute the elements per IFO
{
  int ip=0,jp=0,ifonr=0;
  double **dmatrix  = (double**)calloc(npar,sizeof(double*));
  for(ip=0;ip<npar;ip++) dmatrix[ip]  = (double*)calloc(npar,sizeof(double));
  
  for(ifonr=0;ifonr<networksize;ifonr++) {
    computeFishermatrixIFO(par,npar,ifo,networksize,ifonr,dmatrix);
    
    for(ip=0;ip<npar;ip++) {
      for(jp=0;jp<npar;jp++) {
	matrix[ip][jp] += dmatrix[ip][jp];
      }
    }
  }
}




void printparset(struct parset par)
// Print the parameter set par to screen
{
  printf("\n");
  printf("  Mc:  %20.10lf,   eta:     %20.10lf,   tc:      %20.10lf,   logdl:   %20.10lf\n",par.mc,par.eta,par.tc,par.logdl);
  printf("  a:   %20.10lf,   kappa:   %20.10lf,   longi:   %20.10lf,   sinlati: %20.10lf\n",par.spin,par.kappa,par.longi,par.sinlati);
  printf("  phi: %20.10lf,   sinthJ0: %20.10lf,   phithJ0: %20.10lf,   alpha:   %20.10lf\n",par.phase,par.sinthJ0,par.phiJ0,par.alpha);
  printf("\n");
}
