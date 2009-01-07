#ifndef mcmc_lal_h
#define mcmc_lal_h


#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/GenerateInspiral.h>

#include <mcmc.h>


// Declare function prototypes:
//void LALHpHc(CoherentGW *waveform, double *hplus, double *hcross, int *l, int length, struct parset *par, struct interferometer *ifo, int ifonr);
void LALHpHc12(CoherentGW *waveform, int *l, struct parset *par, struct interferometer *ifo);
void LALHpHc15(CoherentGW *waveform, int *l, struct parset *par, struct interferometer *ifo);
//double LALFpFc(CoherentGW *waveform, double *wave, int *l, int length, struct parset *par, int ifonr);
double LALFpFc(CoherentGW *waveform, double *wave, int length, struct parset *par, int ifonr);
void LALfreedom(CoherentGW *waveform);


#endif
