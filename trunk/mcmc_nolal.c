// File with dummy routine to compile and run the MCMC code without LAL


#include <mcmc.h>


//Use the LAL 3.5/2.5 PN spinning waveform, with 1 spinning object (12 parameters)
void templateLAL12(struct parset *par, struct interferometer *ifo[], int ifonr)
{
  //Get rid of warnings at compile time (x set but never used)
  double x = 0.0;
  x = par->par[0];
  x = ifo[0]->lati; 
  x = (double)ifonr;
  x = x;
  
  printf("\n\n");
  printf("\n    *** ERROR ***");
  printf("\n    You have compiled the MCMC code without LAL support, but are trying to run it with the LAL waveform.");
  printf("\n    Please change mcmcWaveform or injectionWaveform to use a non-LAL waveform  *OR*  compile with LAL support.");
  printf("\n\n\n");
  exit(1);
}

//Use the LAL 3.5/2.5 PN spinning waveform, with 2 spinning objects (15 parameters)
void templateLAL15(struct parset *par, struct interferometer *ifo[], int ifonr)
{
  //Get rid of warnings at compile time (x set but never used)
  double x = 0.0;
  x = par->par[0];
  x = ifo[0]->lati; 
  x = (double)ifonr;
  x = x;
  
  printf("\n\n");
  printf("\n    *** ERROR ***");
  printf("\n    You have compiled the MCMC code without LAL support, but are trying to run it with the LAL waveform.");
  printf("\n    Please change mcmcWaveform or injectionWaveform to use a non-LAL waveform  *OR*  compile with LAL support.");
  printf("\n\n\n");
  exit(1);
}

