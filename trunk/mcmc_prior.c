#include <mcmc.h>



// *** Routines to handle Priors ***

int prior(double *par, int p, int waveformVersion, struct mcmcvariables mcmc)
// Call the prior routine, the variable waveformVersion determines which one
{
  if(waveformVersion==1) {
    return prior1(par, p, mcmc);  // Apostolatos 12-parameter template
  } 
  else if(waveformVersion==2) {
    return prior1(par, p, mcmc);  // LAL 12-parameter template
  }
  else if(waveformVersion==3) {
    return prior2(par, p, mcmc);  // LAL 15-parameter template
  }
  else return 0;
}



double uncorrelated_mcmc_single_update_angle_prior(double sigma, int p, int waveformVersion)
// Call the prior routine, the variable waveformVersion determines which one
{
  if(waveformVersion==1) { // Apostolatos 12-parameter template
    if(p==6 || p==8 || p==10 || p==11) return min(tpi,sigma);  //Bring the sigma between 0 and 2pi;
	else return sigma; 
  } 
  else if(waveformVersion==2) { // LAL 12-parameter template
    if(p==6 || p==8 || p==10 || p==11) return min(tpi,sigma);  //Bring the sigma between 0 and 2pi;
    else return sigma; 
  }
  else if(waveformVersion==3) { // LAL 15-parameter template
    if(p==4 || p==7 || p==8 || p==11 || p==14) return min(tpi,sigma); //Bring the sigma between 0 and 2pi;
	else return sigma; 
  }
else return 0.0;
}




//****************************************************************************************************************************************************  
int prior1(double *par, int p, struct mcmcvariables mcmc) //Apostolatos 12-parameter priors
//****************************************************************************************************************************************************  
//Contains boundary conditions and prior information for the adaptive MCMC.  Trying to avoid returning 0, to increase jump sizes
//  'Stick to the wall method' should not be used, since it is asymmetric
//  0:mc, 1:eta, 2:tc, 3:logdl, 4:spin, 5:kappa, 6: longi (->RA), 7:sindec, 8:phase, 9:sinthJ0, 10:phiJ0, 11:alpha
{
  int prior = 1;
  double dt = 0.0;
  double *lb,*ub;
  lb = (double*)calloc(mcmc.nMCMCpar,sizeof(double));
  ub = (double*)calloc(mcmc.nMCMCpar,sizeof(double));
  
  //Lower and upper boundaries:
  //Mc:
  lb[0] = 1.0;
  ub[0] = 6.0;
  
  
  //eta:
  lb[1] = 0.03;
  ub[1] = 0.25;
  
  //t_c:
  dt = 0.05; //This is dt/2  For known signals
  //dt = 0.5; //This is dt/2  For unknown signals
  lb[2] = prior_tc_mean - dt; //t_c
  ub[2] = prior_tc_mean + dt;
  
  lb[3] = -6.9; //ln(d_L); ln(-6.9) = 0.001Mpc = 1kpc
  ub[3] = 4.6;  //ln(4.6) = 100Mpc
  
  lb[4] = 1.e-10; //a_spin
  ub[4] = 0.999999;
  
  lb[5] = -0.999999; //kappa
  ub[5] = 0.999999;
  
  lb[7] = -0.999999; //sin(dec)
  ub[7] = 0.999999;
  
  lb[9] = -0.999999; //sin(theta_J0)
  ub[9] = 0.999999;
  
  
  if(p==6 || p==8 || p==10 || p==11) {                                    // Periodic boundary condition to bring the variable between 0 and 2pi
    *par = fmod(*par+mtpi,tpi);
  } else {                                                                // Bounce back from the wall
    if(*par<lb[p] || *par>ub[p]) {                                        // Do only one bounce
      if(*par<lb[p]) {
	*par = lb[p] + fabs(*par - lb[p]);
      } else {
	*par = ub[p] - fabs(*par - ub[p]);
      }
      if(*par<lb[p] || *par>ub[p]) prior = 0;                             // If, after bouncing once, still outside the range, reject
    }
  }
  
  free(lb);
  free(ub);
  
  return prior;
}
//End prior
//****************************************************************************************************************************************************  


//****************************************************************************************************************************************************  
int prior2(double *par, int p, struct mcmcvariables mcmc) //LAL 15-parameter priors
//****************************************************************************************************************************************************  
//Contains boundary conditions and prior information for the adaptive MCMC.  Trying to avoid returning 0, to increase jump sizes
//  'Stick to the wall method' should not be used, since it is asymmetric
//  0:mc, 1:eta, 2:tc, 3:logdl, 4:spin, 5:kappa, 6: longi (->RA), 7:sindec, 8:phase, 9:sinthJ0, 10:phiJ0, 11:alpha
{
  int prior = 1;
  double dt = 0.0;
  double *lb,*ub;
  lb = (double*)calloc(mcmc.nMCMCpar,sizeof(double));
  ub = (double*)calloc(mcmc.nMCMCpar,sizeof(double));
  
  //Lower and upper boundaries:
  //Mc:
  lb[0] = 1.0;
  ub[0] = 6.0;

  lb[1] = 0.03;//eta
  ub[1] = 0.25;
  
  //t_c:
  dt = 0.05; //This is dt/2  For known signals
  //dt = 0.5; //This is dt/2  For unknown signals
  lb[2] = prior_tc_mean - dt; //t_c
  ub[2] = prior_tc_mean + dt;
  
  lb[3] = -6.9; //ln(d_L); ln(-6.9) = 0.001Mpc = 1kpc
  ub[3] = 4.6;  //ln(4.6) = 100Mpc
  
  lb[4] = 0.0; //RA
  ub[4] = tpi;
  
  lb[5] = -0.999999; //sin(dec)
  ub[5] = 0.999999;
  
  lb[6] = -0.999999; //cos(i)
  ub[6] = 0.999999;
  
  lb[7] = 0.0; //phi_c
  ub[7] = tpi;
  
  lb[8] = 0.0; //psi
  ub[8] = tpi;
  
  lb[9] = 1.e-10; //a_spin1
  ub[9] = 0.999999;
  
  lb[10] = -0.999999; //cos(theta_1)
  ub[10] = 0.999999;

  lb[11] = 0.0; //phi_1
  ub[11] = tpi;
  
  lb[12] = 1.e-10; //a_spin2
  ub[12] = 0.999999;
  
  lb[13] = -0.999999; //cos(theta_2)
  ub[13] = 0.999999;

  lb[14] = 0.0; //phi_2
  ub[14] = tpi;
   
  //Set prior to 0 if outside range: Seems very inefficient for correlated update proposals
  /*
    if(p==0) if(*par < 1.0 || *par > 5.0) prior = 0;                // Chirp mass <1 or >5
    if(p==1) if(*par < 0.03 || *par > 0.245) prior = 0;                // Mass ratio, refect if <0.03 or >0.245
    if(p==2) if(*par <= prior_tc_mean -0.05 || *par > prior_tc_mean+0.05) prior = 0;    // If time outside range
    if(p==3) if(*par <= -6.91 || *par > 4.6) prior = 0;                  // If distance <1kpc or >100Mpc
    if(p==4) if(*par <= 0.0 || *par > 1.0) prior = 0;                     // If spin<0 or >1
    if(p==5 || p==7 || p==9)  if(*par < -1.0 || *par > 1.0) prior = 0;    // If sin or cos <-1 or >1
  */
  
  //printf("%2d  %20.6f  %20.6f  %20.6f  ",p,lb[p],ub[p],*par);
  
  
 if(p==4 || p==7 || p==8 || p==11 || p==14) {                                    // Periodic boundary condition to bring the variable between 0 and 2pi
    *par = fmod(*par+mtpi,tpi);
  } else {																
  

																		    // Bounce back from the wall
    //while(*par<lb[p] || *par>ub[p]) {                                     // Do as many bounces as necessary to get between the walls
    //  if(*par<lb[p])  *par = lb[p] + fabs(*par - lb[p]);
    //  if(*par>ub[p])  *par = ub[p] - fabs(*par - ub[p]);
    //}
           if(*par<lb[p] || *par>ub[p]) {                                        // Do only one bounce
      if(*par<lb[p]) {
	*par = lb[p] + fabs(*par - lb[p]);
      } else {
	*par = ub[p] - fabs(*par - ub[p]);
      }
      if(*par<lb[p] || *par>ub[p]) prior = 0;                             // If, after bouncing once, still outside the range, reject
    }
	 
  }
  
  //printf("%20.6f  %d  \n",*par,prior);
  
  free(lb);
  free(ub);
  
  return prior;
}
//End prior
//****************************************************************************************************************************************************  


