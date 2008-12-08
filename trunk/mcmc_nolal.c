// File with dummy routine to compile and run the MCMC code without LAL


#include <mcmc.h>


//Use the LAL 3.5/2.5 PN spinning waveform, with 2 spinning objects (15 parameters)
void template15(struct parset *par, struct interferometer *ifo[], int ifonr)
{
  double x = 0.0;
  x = par->m1;
  x = ifo[0]->lati; 
  x = (double)ifonr;
  
  //Get rid of warnings at compile time (x set but never used)
  double y = x;
  x = y;
  
  printf("\n\n");
  printf("\n    *** ERROR ***");
  printf("\n    You have compiled the MCMC code without LAL support, but are trying to run it with the LAL waveform.");
  printf("\n    Please set  waveformversion = 1  in the main input file (e.g. mcmc.input),  *OR*  compile with LAL support.");
  printf("\n\n\n");
  exit(1);
}

