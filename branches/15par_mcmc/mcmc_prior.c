#include <mcmc.h>



// *** Routines to handle Priors ***

int prior(double *par, int p)
// Call the prior routine, the global variable waveformversion determines which one
{
  if(waveformversion==1) {
    prior1(par, p);  // Apostolatos 12-parameter template
  } 
  else if(waveformversion==2) {
    prior1(par, p);  // LAL 12-parameter template
  }
  else if(waveformversion==3) {
    prior2(par, p);  // LAL 15-parameter template
  }

}



void setrandomtrueparameters(struct runpar *run)
// Call the prior for the random injection, the global variable waveformversion determines which one
{
  if(waveformversion==1) {
    setrandomtrueparameters1(run);  // Apostolatos 12-parameter template
  } 
  else if(waveformversion==2) {
    setrandomtrueparameters1(run);  // LAL 12-parameter template
  }
  else if(waveformversion==3) {
    setrandomtrueparameters2(run);  // LAL 15-parameter template
  }

}



//****************************************************************************************************************************************************  
int prior1(double *par, int p) //Apostolatos 12-parameter priors
//****************************************************************************************************************************************************  
//Contains boundary conditions and prior information for the adaptive MCMC.  Trying to avoid returning 0, to increase jump sizes
//  'Stick to the wall method' should not be used, since it is asymmetric
//  0:mc, 1:eta, 2:tc, 3:logdl, 4:spin, 5:kappa, 6: longi (->RA), 7:sindec, 8:phase, 9:sinthJ0, 10:phiJ0, 11:alpha
{
  int prior = 1;
  double dt = 0.0;
  double *lb,*ub;
  lb = (double*)calloc(npar,sizeof(double));
  ub = (double*)calloc(npar,sizeof(double));
  
  //Lower and upper boundaries:
  //Mc:
  lb[0] = 1.0;
  ub[0] = 6.0;
  //lb[0] = 0.609; //Mc for 1.4+1.4Mo: /2 - *2:  1.4+1.4->Mc=1.22; use range Mc/2 - Mc*2: 0.609-2.44
  //ub[0] = 2.44;
  //lb[0] = 1.30; //Mc for 3+3Mo: /2 - *2
  //ub[0] = 5.22;
  //lb[0] = 1.50; //Mc for 10+1.4Mo: /2 - *2
  //ub[0] = 5.99;
  //lb[0] = 4.35; //Mc for 10+10Mo: /2 - *2
  //ub[0] = 17.4;

  
  //eta:
  //lb[1] = 0.001; //Chains get stuck at low \eta and very high L!
  lb[1] = 0.03;
  //ub[1] = 0.245;
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
  
  
  if(p==6 || p==8 || p==10 || p==11) {                                    // Periodic boundary condition to bring the variable between 0 and 2pi
    *par = fmod(*par+mtpi,tpi);
  } else {                                                                // Bounce back from the wall
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


void setrandomtrueparameters1(struct runpar *run)  //Get random values for the 'true' parameters for the 12-parameter spinning template. Contain priors for the injection, not the MCMC. 
// *** This changes the injected signal!!! ***
{
  int i=0;
  gsl_rng *ran;
  double rannr = 0.0;
  ran = gsl_rng_alloc(gsl_rng_mt19937);  // GSL random-number seed
  if(1==2) {  //Select a random seed, *** ONLY FOR TESTING ***
    printf("\n  *** SELECTING RANDOM SEED ***  This should only be done while testing!!! setrandomtrueparameters() \n\n");
    run->ranparseed = 0;
    setseed(&run->ranparseed);
    //printf("  Seed: %d\n", run->ranparseed);
  }
  gsl_rng_set(ran, run->ranparseed);     // Set seed
  
  //Lower and upper boundaries:
  double *lb,*ub,db,dt;
  lb = (double*)calloc(npar,sizeof(double));
  ub = (double*)calloc(npar,sizeof(double));
  lb[0] = 2.0;       //M1 (Mo)
  ub[0] = 4.0;
  lb[1] = 0.05;       //M2 (Mo)
  ub[1] = 0.15;
  dt = 0.5; //This is dt/2
  lb[2] = prior_tc_mean - dt; //t_c
  ub[2] = prior_tc_mean + dt;
  lb[3] = 2.3;      //d_L (Mpc)  !Linear!
  ub[3] = 3.4;
  lb[4] = 1.e-10;    //a_spin (0-1)
  ub[4] = 0.999999;
  lb[5] = -0.999999; //kappa
  ub[5] = 0.999999;
  lb[6] = 0.0;       //RA (h)
  ub[6] = tpi;
  lb[7] = -0.999999; //sin(dec)
  ub[7] = 0.999999;
  lb[8] = 0.0;       //phi_c (deg)
  ub[8] = tpi;
  lb[9] = -0.999999; //sin(theta_J0)
  ub[9] = 0.999999;
  lb[10] = 0.0;      //phi_Jo (deg)
  ub[10] = tpi;
  lb[11] = 0.0;      //alpha_c (deg)
  ub[11] = tpi;
  
  for(i=0;i<npar;i++) {
    db = ub[i]-lb[i];
    rannr = gsl_rng_uniform(ran);                                                        //This assures you always draw the same number of random variables
    if(run->setranpar[i]==1) truepar[i] = rannr*db + lb[i];
   // if(i==5 && run->setranpar[i]==1) truepar[i] = acos(rannr*2.0 - 1.0)*r2d;             //kappa -> th_SL
   // if((i==7 || i==9)  && run->setranpar[i]==1) truepar[i] = asin(rannr*2.0 - 1.0)*r2d;  //sin(dec)->dec, sin(th_J0)->th_J0
    //printf("  %d  %lf  %lf  %lf  %lf\n",i,lb[i],ub[i],db,truepar[i]);
  }
  
  free(lb);
  free(ub);
  gsl_rng_free(ran);
}




//****************************************************************************************************************************************************  
int prior2(double *par, int p) //LAL 15-parameter priors
//****************************************************************************************************************************************************  
//Contains boundary conditions and prior information for the adaptive MCMC.  Trying to avoid returning 0, to increase jump sizes
//  'Stick to the wall method' should not be used, since it is asymmetric
//  0:mc, 1:eta, 2:tc, 3:logdl, 4:spin, 5:kappa, 6: longi (->RA), 7:sindec, 8:phase, 9:sinthJ0, 10:phiJ0, 11:alpha
{
  int prior = 1;
  double dt = 0.0;
  double *lb,*ub;
  lb = (double*)calloc(npar,sizeof(double));
  ub = (double*)calloc(npar,sizeof(double));
  
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
  
  lb[10] = 0.0; //theta_1
  ub[10] = tpi;

  lb[11] = 0.0; //phi_1
  ub[11] = pi;
  
  lb[12] = 1.e-10; //a_spin2
  ub[12] = 0.999999;
  
  lb[13] = 0.0; //theta_2
  ub[13] = tpi;

  lb[14] = 0.0; //phi_2
  ub[14] = pi;
   
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
  
  
 if(p==4 || p==7 || p==8 || p==10 || p==13) {                                    // Periodic boundary condition to bring the variable between 0 and 2pi
    *par = fmod(*par+mtpi,tpi);
  } else {                                                                // Bounce back from the wall
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


void setrandomtrueparameters2(struct runpar *run)  //Get random values for the 'true' parameters for the 15-parameter spinning template. Contain priors for the injection, not the MCMC. 
// *** This changes the injected signal!!! ***
{
  int i=0;
  gsl_rng *ran;
  double rannr = 0.0;
  ran = gsl_rng_alloc(gsl_rng_mt19937);  // GSL random-number seed
  if(1==2) {  //Select a random seed, *** ONLY FOR TESTING ***
    printf("\n  *** SELECTING RANDOM SEED ***  This should only be done while testing!!! setrandomtrueparameters() \n\n");
    run->ranparseed = 0;
    setseed(&run->ranparseed);
    //printf("  Seed: %d\n", run->ranparseed);
  }
  gsl_rng_set(ran, run->ranparseed);     // Set seed
  
  //Lower and upper boundaries:
  double *lb,*ub,db,dt;
  lb = (double*)calloc(npar,sizeof(double));
  ub = (double*)calloc(npar,sizeof(double));

  //Mc:
  lb[0] = 2.0;
  ub[0] = 4.0;

  lb[1] = 0.05;//eta
  ub[1] = 0.15;
  
  //t_c:
  dt = 0.5; //This is dt/2
  lb[2] = prior_tc_mean - dt;
  ub[2] = prior_tc_mean + dt;
  
  lb[3] = 2.3; // 10Mpc
  ub[3] = 3.4; // 30Mpc
  
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
  
  lb[10] = 0.0; //theta_1
  ub[10] = tpi;

  lb[11] = 0.0; //phi_1
  ub[11] = pi;
  
  lb[12] = 1.e-10; //a_spin2
  ub[12] = 0.999999;
  
  lb[13] = 0.0; //theta_2
  ub[13] = tpi;

  lb[14] = 0.0; //phi_2
  ub[14] = pi;
 
  
  for(i=0;i<npar;i++) {
    db = ub[i]-lb[i];
    rannr = gsl_rng_uniform(ran);                                                        //This assures you always draw the same number of random variables
    if(run->setranpar[i]==1) truepar[i] = rannr*db + lb[i];
	//printf("  %d  %lf  %lf  %lf  %lf\n",i,lb[i],ub[i],db,truepar[i]);
  }
  
  free(lb);
  free(ub);
  gsl_rng_free(ran);
}


