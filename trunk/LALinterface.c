


#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/GenerateInspiral.h>
//#include <lal/LALConstants.h>
#include <lal/DetectorSite.h>

#include <lal/Date.h>
#include <lal/TimeDelay.h>
#include <lal/SkyCoordinates.h>
//#include <lal/DetectorSite.h>

//////////////////////////////////////////
#include <mcmc.h>
//////////////////////////////////////////

//NRCSID(LALSTPNWaveformTestC, "$Id: LALSTPNWaveformTest.c,v 1.1 2004/05/05 20:06:23 thomas Exp");



//void LALHpHc(CoherentGW *waveform, double *hplus, double *hcross, int *l, int length, struct parset *par, struct interferometer *ifo, int ifonr) { //MvdS: ifonr not used //Vivien: just for debugging in a printf
void LALHpHc(CoherentGW *waveform, double *hplus, double *hcross, int *l, int length, struct parset *par, struct interferometer *ifo) {
  // Compute h_+ and h_x
  
  static LALStatus    mystatus;
  
  SimInspiralTable    injParams;
  PPNParamStruc       ppnParams;
  
  //const char        *filename = "wave2.dat";
  //FILE        *outputfile;
  INT4        i;
  int			lengthLAL;
  //REAL8       dt;
  REAL8       a1, a2, phi, shift; //, phi0;
  //REAL8		M1,M2,t_c,d_L,a_spin,th_SL,H_A,dec,phic,th_Jo,phi_Jo,alpha;
  //REAL8		epsilon,gamma,kappa,omega,l,j,s,r,M,f_lower;
  
  pi=M_PI;
  
  
  
  double x;
  double m1=0.0,m2=0.0,M=0.0,mu=0.0;
  double cvec1[3],cvec2[3],cvec3[3];
  double tvec1[3],tvec2[3],tvec4[3],tvec6[3],tvec7[3];
  //int /*i,*/ terminate=0;
  double n_L[3];
  double spin=0.0,samplerate=0.0,inversesamplerate=0.0;//,altitude=0.0,azimuth=0.0,localtc=0.0;
  //double Ltheta0=0.0, Lphi0=0.0, Lphi02=0.0;
  for(i=0;i<3;i++) {
    cvec1[i] = 0.0;
    cvec2[i] = 0.0;
    cvec3[i] = 0.0;
    tvec1[i] = 0.0;
    tvec2[i] = 0.0;
    //tvec3[i] = 0.0;
    tvec4[i] = 0.0;
    //tvec5[i] = 0.0;
    tvec6[i] = 0.0;
    tvec7[i] = 0.0;
    n_L[i] = 0.0;
  }
  
  
  double f_lower=ifo->lowCut;
  
  samplerate = (double)ifo->samplerate;
  inversesamplerate = 1.0/samplerate;
  
  //double altitude   = par->localti[ifonr];
  //double azimuth    = par->locazi[ifonr];
  
  //struct parset par;
  
  //struct parset *par;
  
  //	 par = malloc(sizeof(struct parset));
  
  
  
  memset( &mystatus, 0, sizeof(LALStatus) );
  // memset( &thewaveform, 0, sizeof(CoherentGW) );
  memset( waveform, 0, sizeof(CoherentGW));
  memset( &injParams, 0, sizeof(SimInspiralTable) );
  memset( &ppnParams, 0, sizeof(PPNParamStruc) );
  
  /* --- parameters from Marc's code.... --- */
  
  /*	par->m1       = 10.0;                    // M1 (10.0)
	par->m2       = 1.4;                    // M2  (1.4)
	
	par->m        = par->m1+par->m2;
	
	par->mu       = par->m1*par->m2/par->m;
	//fprintf(stderr, "par->mu = %f\n", par->mu);
	par->eta      = par->mu/par->m;				// mass ratio                
	par->mc       = par->m*pow(par->eta,0.6);      // chirp mass. in Mo         
	par->tc       = 0.0;                    // coalescence time
	par->logdl    = log(13.0);               // log-distance (Mpc) (17.5)             
	
	par->spin     = 0.00001;                    // magnitude of total spin   (0.1)
	par->kappa    = cos(55.0*d2r);           // L^.S^, cos of angle between L^ & S^  (0.819152)
	par->longi    = 91.0*d2r;                // longitude (~"hour angle"?, saved as RA)     (120)
	par->sinlati  = sin(40.0*d2r);           // sin latitude (sin(delta))  (40)     
	
	par->phase    = 11.4592*d2r;                // orbital phase   (phi_c)   (0.2)
	par->sinthJ0  = sin(15.0*d2r);           // sin Theta_J0 ~ polar, 0=NP    (15)
	par->phiJ0    = 125.0*d2r;               // Phi_J0 ~ azimuthal            (125)
	par->alpha    = 51.5662*d2r;               // Alpha_c                       (0.9 rad = 51.566202deg)
	
	par->loctc    = NULL;
	par->localti  = NULL;
	par->locazi   = NULL;
	par->locpolar = NULL;*/
  
  
  
  double n_z[3] = {0.0,0.0,1.0};                                                                                        //North in global coordinates
  double normalvec[3];//,rightvec[3],leftvec[3];                                                                                                  
  for(i=0;i<3;i++) normalvec[i] = ifo->normalvec[i];                                                             //Detector position normal vector = local zenith vector z'
  //  for(i=0;i<3;i++) rightvec[i] = ifo->rightvec[i];
  //  for(i=0;i<3;i++) leftvec[i] = ifo->leftvec[i];
  
  // for(i=0;i<3;i++) normalvec[i]=n_z[i];
  
  
  //double D_L = exp(par->logdl)*Mpcs;                                                                                    //Source luminosity distance, in seconds
  double coslati = sqrt(1.0-par->sinlati*par->sinlati);
  // double n_N[3] = {sin(par->longi)*coslati,cos(par->longi)*coslati,par->sinlati};                                       //n_N: Position unit vector = N^
  //  double n_N[3] = {cos(par->longi)*par->sinlati,sin(par->longi)*par->sinlati,coslati};
  double n_N[3] = {cos(par->longi)*coslati,sin(par->longi)*coslati,par->sinlati};
  
  double sthJ0   = par->sinthJ0;                                                                                        //n_J0: 'total' AM unit vector, J0^  (almost equal to the real J, see Eq.15)
  double cthJ0   = sqrt(1. - sthJ0*sthJ0);
  // double n_J0[3] = {sin(par->phiJ0)*cthJ0,cos(par->phiJ0)*cthJ0,sthJ0};  //???  ----- Swap sin(phi) & cos(phi)?  Or keep (by definition equal to n_N?
  // double n_J0[3] = {cos(par->phiJ0)*sthJ0,sin(par->phiJ0)*sthJ0,cthJ0};
  
  double n_J0[3] = { cos(par->phiJ0)*cthJ0 , sin(par->phiJ0)*cthJ0 , sthJ0 };
  
  //Get masses from Mch and eta
  double root = sqrt(0.25-par->eta);
  double fraction = (0.5-root) / (0.5+root);
  double inversefraction = 1.0/fraction;
  double Mc = par->mc*M0;                                                                                               //Chirp mass in seconds
  x = exp(0.6*log(fraction));
  m1 = Mc * (pow(1.0+fraction,0.2) / x);
  m2 = Mc * (pow(1.0+inversefraction,0.2) * x);
  M = m1+m2;                                                                                                            // Eq.16a
  mu = m1*m2/M;                                                                                                         // Eq.16b
  spin = par->spin*m1*m1;
  
  
  //double beta = 1.0/12.0*(113.0*(m1*m1)/(M*M) + 75.0*par->eta)*par->kappa*spin/(m1*m1);                                 // Eq.20, for S2=0 or m1=m2,S1=S2:  kappa*spin/(m1*m1) = L^.S/m1^2
  
  //double cst1 = 743.0/336.0 + 11.0/4.0*par->eta;
  //double cst2 = (4.0*pi-beta);
  double cst5 = spin*sqrt(1.0-par->kappa*par->kappa);
  //printf("%10.20f\n", cst5);
  
  
  //Constant vector 1 for the construction of Eq.41e
  //facvec(n_J0,-cthJ0,tvec1);      
  facvec(n_J0,-sthJ0,tvec1);      
  //tvec1 = -J0^*cos(theta_J0)
  addvec(n_z,tvec1,tvec2);             
  //tvec2 = n_z - J0^*cos(theta_J0)
  //facvec(tvec2,1.0/sthJ0,cvec1);                                                                                         //cvec1 = (n_z - J0^*cos(theta_J0))/sin(theta_J0)
  facvec(tvec2,1.0/cthJ0,cvec1);
  
  //Constant vector 2 for the construction of Eq.41e
  crossproduct(n_J0,n_z,tvec1);      
  //tvec1 = J0^ x z^
  //facvec(tvec1,1.0/sthJ0,cvec2);                                                                                         //cvec2 = (J0^ x z^) / sin(theta_J0)
  facvec(tvec1,1.0/cthJ0,cvec2);
  
  
  //Constant vector 3 for the construction of Eq.12
  facvec(n_N,-dotproduct(normalvec,n_N),tvec1);                                                                          //tvec1 = -N^(z^'.N^)
  addvec(normalvec,tvec1,cvec3);                                                                                         //cvec3 = z^' - N^(z^'.N^)
  
  //Construct Eq.8ab
  //double cosalti   = cos(altitude);
  //double sin2azi   = sin(2.0*azimuth);
  //double cos2azi   = cos(2.0*azimuth);
  //double cst6  = 0.5*(1.0+cosalti*cosalti)*cos2azi;
  //double cst7  = cosalti*sin2azi;
  
  /***************************************************************************FIXED FREQUENCIES******************************************************************/
  //double omega_low  = pi*f_lower;//ifo[ifonr]->lowCut;   //30 or 40 Hz, translated from f_gw to omega_orb
  //double omega_high = pi*ifo->highCut;  //1600 Hz, translated from f_gw to omega_orb
  //double omega_high = min(pi*ifo[ifonr]->highCut, exp(-1.5*log(cutoff_a) - log(M)) );  //1600 Hz, translated from f_gw to omega_orb, or a/M = cutoff_a, whichever is smaller
  
  
  /*  if(printmuch) {
      double Tcoal = 5.0*pow(8.0*omega_low,-8.0*c3rd)*pow(Mc,-5.0*c3rd) * (1.0 + 4.0*c3rd*cst1*pow(omega_low*M,2.0*c3rd) - 1.6*cst2*(omega_low*M));   //Time between f_low and coalescence
      double t0 = localtc - Tcoal;
      double deltat = (double)length*inversesamplerate;
      printf("Times:  Tcoal: %g,  t0: %g,  localtc: %g,  tstart: %g,  length: %d,  dt: %g,  dt-ltc: %g\n",Tcoal,t0,localtc,tstart,length,deltat,deltat-localtc);
      }*/
  
  //double oldomega = -1.e30;
  double alpha=0.0;//,phi_gw=0.0;
  
  //To display the number of wave and precession cycles:
  //double phi1 = 0.0;
  //double phi2 = 0.0;
  //double alpha1 = 0.0;
  //double alpha2 = 0.0;
  //int i1=0,i2=0;
  
  //double t=0.0,tau=0.0,tau18=0.0,tau28=0.0,tau38=0.0,tau58=0.0,tau_18=0.0,tau_28=0.0,tau_38=0.0,tau_58=0.0,tau_68=0.0;
  double omega_orb=0.0,l_L=0.0,Y=0.0,Gsq=0.0,G=0.0,slamL=0.0,clamL=0.0,LdotN=0.0;
  //double locpolar=0.0,sin2polar=0.0,cos2polar=0.0,Fplus=0.0,Fcross=0.0;
  double cst4=0.0,x1=0.0,x2=0.0,x3=0.0;
  // double taperx[length],omegas[length];
  // for(i=0;i<length;i++) {
  //   taperx[i] = 0.0;
  //   omegas[i] = 0.0;
  // }
  
  /*-- fill 'output' with time-domain template: --*/
  //  for (i=0; i<length; ++i){
  /// determine time left until coalescence, "(t_c-t)" in (4.17)/(11): 
  /*    t = localtc - ((double)i)*inversesamplerate;  // (time to t_c) = "(t_c-t)" in (4.17)
	if(t<0.0) { 
	terminate = 1;
	}
	else {
	tau    = par->eta/(5.0*M)*t;  //t = localtc-t already
	tau18  = exp(0.125*log(tau));
	tau28  = tau18*tau18;
	tau38  = tau28*tau18;
	tau58  = tau28*tau38;
	tau_18 = 1.0/tau18;
	tau_28 = tau_18*tau_18;
	tau_38 = tau_18*tau_28;
	tau_58 = tau_28*tau_38;
	tau_68 = tau_38*tau_38;
	
	omega_orb = 1.0/(8.0*M) * (tau_38 + 0.125*cst1*tau_58 - 0.075*cst2*tau_68);   // Orbital frequency
	omegas[i] = omega_orb;
	}
	
	if ((omega_orb>=omega_low) && (terminate==0)) {  // After source comes into window, before tc and before frequency reaches its maximum  Careful, if t<0, omega = nan!!!
	
	if (omega_orb < oldomega || omega_orb >= omega_high){  // Frequency starts decreasing, or frequency higher than highCut --> terminate signal
  */      //if (omega_orb < oldomega || omega_orb >= omega_high || taperx[i]>0.09){  // Frequency starts decreasing, or frequency higher than highCut, or v_orb>0.3c --> terminate signal
  //     ifo[ifonr]->FTin[i] = 0.0; /****************************************COMMENTED*********************************************/
  //terminate = 1;
  /*        if(omega_orb < oldomega) terminate += 2;
	    if(omega_orb >= omega_high) terminate += 4;
	    //if(taperx[i]>0.09) terminate += 8;
	    }
	    
	    else {              // Frequency still increasing --> keep on computing...
	    if(i1==0) i1=i;  //Save initial i for tapering the beginning of the signal
	    i2 = i;          //Save final i for tapering the end of the signal
	    oldomega = omega_orb;
	    taperx[i] = exp(2.0*c3rd*log(M*omega_orb));                                                                      // x := (M*w)^(2/3)  =  v_orb^2
  */      
  //Compute orbital A.M.
  omega_orb=pi*f_lower;
  
  l_L = m1*m2*exp(-c3rd*log(omega_orb*M));
  
  //GW and orbital phase
  //  phi_gw = par->phase - 2.0/par->eta * (tau58 + 0.625*c3rd*cst1*tau38 - 0.1875*cst2*tau28);                       // GW phase
  //if(fabs(phi1)<1.e-30) phi1 = phi_gw;  //Save initial phi
  //phi2 = phi_gw;                       //Save final phi
  
  Y = spin/l_L;                                                                                                    //Y = |S|/|L|, Eq.43
  Gsq = 1.0 + 2.0*par->kappa*Y + Y*Y;                                                                                      //G^2, Eq.46
  G   = sqrt(Gsq);
  
  cst4 = l_L+par->kappa*spin;
  x = mu*M;
  x1 = x*x*x;
  x = G*l_L;
  x2 = x*x*x;
  x3 = spin*spin*spin;
  alpha = par->alpha - 5.0/(96.0*x1) * (1.0+0.75*m2/m1) * 
    (2.0*x2 - 3.0*par->kappa*spin*cst4*G*l_L - 3.0*par->kappa*x3*(1.0-par->kappa*par->kappa) * asinh(cst4/cst5));                                            //Eq.47
  //if(fabs(alpha1)<1.e-30) alpha1 = alpha;  //Save initial alpha
  //alpha2 = alpha;                         //Save final alpha
  
  slamL = cst5/(l_L*G);                                                                                            //sin(lambda_L), Eq.48a
  clamL = cst4/(l_L*G);                                                                                            //cos(lambda_L), Eq.48b
  
  //Construct Eq.41e
  facvec(n_J0,clamL,tvec1);                                                                                        //tvec1 = J0^*cos(lambda_L)
  facvec(cvec1,slamL*cos(alpha),tvec4);                                                                            //tvec4 = (n_z - J0^*cos(theta_J0))*sin(lambda_L)*cos(alpha)/sin(theta_J0)
  facvec(cvec2,slamL*sin(alpha),tvec6);                                                                            //tvec6 = (J0^ x z^) * sin(lambds_L)*sin(alpha)/sin(theta_J0)
  addvec(tvec1,tvec4,tvec7);                                                                                       //Construct Eq.41e
  addvec(tvec7,tvec6,n_L);
  //Eq.41e: n_L=L^
  //	printf("n_N = %10.20f\n", sqrt(n_N[0]*n_N[0]+n_N[1]*n_N[1]+n_N[2]*n_N[2]));
  //	printf("n_J0 = %10.20f\n", sqrt(n_J0[0]*n_J0[0]+n_J0[1]*n_J0[1]+n_J0[2]*n_J0[2]));
  //	printf("n_L = %10.20f\n", sqrt(n_L[0]*n_L[0]+n_L[1]*n_L[1]+n_L[2]*n_L[2]));
  
  //normalise(n_L);
  
  LdotN  = dotproduct(n_L,n_N);
  
  //printf("%10.20f\n", LdotN);
  //	printf("%10.20f\t%10.20f\t%10.20f\n", tvec4[0], tvec4[1], tvec4[2]);
  //printf("%10.20f\t%10.20f\t%10.20f\n", n_N[0], n_N[1], n_N[2]);
  //printf("%10.20f\n", sqrt(n_N[0]*n_N[0]+n_N[1]*n_N[1]+n_N[2]*n_N[2]));
  //printf("%10.20f\t%10.20f\t%10.20f\n", n_L[0], n_L[1], n_L[2]);
  //printf("%10.20f\n", sqrt(n_L[0]*n_L[0]+n_L[1]*n_L[1]+n_L[2]*n_L[2]));
  //printf("%10.20f\n", sqrt(n_J0[0]*n_J0[0]+n_J0[1]*n_J0[1]+n_J0[2]*n_J0[2]));
  
  //crossproduct(n_L,normalvec,tvec1);
  //locpolar  = atan2(dotproduct(n_L,cvec3),dotproduct(n_N,tvec1));
  
  
  
  //      sin2polar = sin(2.0*locpolar);
  //     cos2polar = cos(2.0*locpolar);
  //	cos2polar = sqrt(1.0-sin2polar*sin2polar);                                                                      //Since 2*locpolar should be between -pi and pi (?)
  //     Fplus     =  cst6*cos2polar + cst7*sin2polar;                                                                   //Eq.8a
  //     Fcross    = -cst6*sin2polar + cst7*cos2polar;                                                                   //Eq.8b
  
  
  
  //double polarph = atan((2*LdotN*Fcross)/((1+LdotN)*(1+LdotN)*Fplus));
  
  
  //M=M1+M2;
  //omega=PI*f_lower;
  double r,e;
  
  r=pow(M/(omega_orb*omega_orb),1.0/3.0);
  //s=a_spin*M1*M1;
  //l=M1*M2*sqrt(r)/sqrt(M);
  //gamma=s/l;
  //kappa=cos(th_SL);
  
  e=(16.0/5.0)*sqrt((M/r)*(M/r)*(M/r))/((1+(3.0/4.0)*m2/m1)*(1.0+2.0*par->kappa*Y+Y*Y));
  //printf("%10.20f\n", e);
  //j=sqrt(l*l+s*s+2*l*s*kappa);
  
  double n_S[3];
  
  //sqrt(spin*spin+l_L*l_L+2*spin*l_L*par->kappa)*n_J0[0]-l_L*n_L[0];
  
  n_S[0]=(-(-n_J0[0]-n_J0[1]*n_L[2]*e+n_L[1]*n_J0[2]*e-n_L[0]*n_L[0]*n_J0[0]*e*e-n_L[0]*n_L[1]*n_J0[1]*e*e-n_L[0]*n_L[2]*n_J0[2]*e*e)/(1+n_L[0]*e*e*e+n_L[1]*n_L[1]*e*e+n_L[2]*n_L[2]*e*e));//-n_L[0];
  
  
  n_S[1]=(-(-n_J0[1]+n_J0[0]*n_L[2]*e-n_L[0]*n_J0[2]*e-n_L[0]*n_J0[0]*n_L[1]*e*e-n_L[1]*n_L[1]*n_J0[1]*e*e-n_L[1]*n_L[2]*n_J0[2]*e*e)/(1+n_L[0]*e*e*e+n_L[1]*n_L[1]*e*e+n_L[2]*n_L[2]*e*e));//-n_L[1];
  
  
  n_S[2]=(-(-n_J0[2]-n_J0[0]*n_L[1]*e+n_L[0]*n_J0[1]*e-n_L[0]*n_J0[0]*n_L[2]*e*e-n_L[1]*n_J0[1]*n_L[2]*e*e-n_L[2]*n_L[2]*n_J0[2]*e*e)/(1+n_L[0]*e*e*e+n_L[1]*n_L[1]*e*e+n_L[2]*n_L[2]*e*e));//-n_L[2];
  
  n_S[0]=n_S[0]*G*l_L-l_L*n_L[0];
  n_S[1]=n_S[1]*G*l_L-l_L*n_L[1];
  n_S[2]=n_S[2]*G*l_L-l_L*n_L[2];
  
  /*
    n_S[0]=n_J0[0]*G*l_L-l_L*n_L[0];
    n_S[1]=n_J0[1]*G*l_L-l_L*n_L[1];
    n_S[2]=n_J0[2]*G*l_L-l_L*n_L[2];
  */
  
  //printf("%10.20f\n", spin);
  //printf("computed a = %10.20f\n", sqrt(n_S[0]*n_S[0]+n_S[1]*n_S[1]+n_S[2]*n_S[2])/(m1*m1));
  //printf("%10.20f\n", sqrt(n_L[0]*n_L[0]+n_L[1]*n_L[1]+n_L[2]*n_L[2]));
  //printf("%10.20f\n", sqrt(n_J0[0]*n_J0[0]+n_J0[1]*n_J0[1]+n_J0[2]*n_J0[2]));
  
  /*
    printf("%10.20f\n", sqrt(rightvec[0]*rightvec[0]+rightvec[1]*rightvec[1]+rightvec[2]*rightvec[2]));
    printf("%10.20f\t%10.20f\t%10.20f\n", rightvec[0], rightvec[1], rightvec[2]);
    printf("%10.20f\n", sqrt(leftvec[0]*leftvec[0]+leftvec[1]*leftvec[1]+leftvec[2]*leftvec[2]));
    printf("%10.20f\t%10.20f\t%10.20f\n", leftvec[0], leftvec[1], leftvec[2]);
    printf("%10.20f\n", sqrt(normalvec[0]*normalvec[0]+normalvec[1]*normalvec[1]+normalvec[2]*normalvec[2]));
    printf("%10.20f\t%10.20f\t%10.20f\n", normalvec[0], normalvec[1], normalvec[2]);
    
    double normtemp[3];
    crossproduct(rightvec,leftvec,normtemp);		
    printf("%10.20f\t%10.20f\t%10.20f\n", normtemp[0], normtemp[1], normtemp[2]);
    
  */	
  
  
  double xloc[3],yloc[3],zloc[3];                        // coordinates in the global frame (in which N is defined) of the local vectors (e.g. z=N)                                                                          
  for(i=0;i<3;i++) zloc[i] = n_N[i];                                                             
  for(i=0;i<3;i++) yloc [i] = 0.0;
  for(i=0;i<3;i++) xloc[i] = 0.0;
  
  //printf("zloc %10.20f\t%10.20f\t%10.20f\n", zloc[0], zloc[1], zloc[2]);
  
  crossproduct(zloc,n_L,yloc);
  normalise(yloc);
  crossproduct(yloc,zloc,xloc);
  
  //	printf("zloc      %10.20f\t%10.20f\t%10.20f\n", zloc[0], zloc[1], zloc[2]);
  
  //	double zlocverif[3];
  //	crossproduct(xloc,yloc,zlocverif);
  //	printf("zlocverif %10.20f\t%10.20f\t%10.20f\n", zlocverif[0], zlocverif[1], zlocverif[2]);
  
  //	double n_Lloc[3];
  double n_Sloc[3];
  normalise(n_S);
  
  n_Sloc[0] = dotproduct(n_S,xloc);
  n_Sloc[1] = dotproduct(n_S,yloc);
  n_Sloc[2] = dotproduct(n_S,zloc);
  
  //	n_Lloc[0] = dotproduct(n_L,xloc);
  //	n_Lloc[1] = dotproduct(n_L,yloc);
  //	n_Lloc[2] = dotproduct(n_L,zloc);
  
  //	for(i=0;i<3;i++) n_L[i] = n_Lloc[i];
  for(i=0;i<3;i++) n_S[i] = n_Sloc[i];
  
  
  //	printf("L %10.20f\t%10.20f\t%10.20f\n", n_L[0], n_L[1], n_L[2]);
  
  
  
  //	printf("kloc = %10.20f\n", dotproduct(n_Lloc,n_Sloc));
  //	printf("k = %10.20f\n", dotproduct(n_L,n_S));//n_S[0]*n_L[0]+n_S[1]*n_L[1]+n_S[2]*n_L[2]);
  //	printf("k = %10.20f\n", par->kappa);
  
  n_S[0]=par->spin*n_S[0];
  n_S[1]=par->spin*n_S[1];
  n_S[2]=par->spin*n_S[2];
  
  //printf("par->spin = %f\n", par->spin);
  /*
    Ltheta0 = acos(n_L[2]);
    
    if(n_L[0]<0) {	
    Lphi0 = pi-asin(n_L[1]/sin(Ltheta0));
    }
    else {
    Lphi0 = asin(n_L[1]/sin(Ltheta0));
    }
    
    if(n_L[1]<0) {	
    Lphi02 = -acos(n_L[0]/sin(Ltheta0));
    }
    else {
    Lphi02 = acos(n_L[0]/sin(Ltheta0));
    }
  */
  
  //printf("%10.20f\t%10.20f\t%10.20f\n", n_L[0], n_L[1], n_L[2]);
  //printf("%10.20f\t%10.20f\t%10.20f\n", cos(Lphi0)*sin(Ltheta0),sin(Lphi0)*sin(Ltheta0),cos(Ltheta0));
  
  /* --- first we fill the SimInspiral structure --- */
  
  injParams.mass1 = (float)par->m1;//M1;
  
  injParams.mass2 = (float)par->m2;//M2;
  
  
  /* MV-20060224: I believe this is not used in the SpinTaylor code! */
  injParams.f_final = 1600.0;
  
  injParams.f_lower = (float)f_lower;
  
  
  LALSnprintf(injParams.waveform,LIGOMETA_WAVEFORM_MAX*sizeof(CHAR),"SpinTaylorthreePointFivePN");//"SpinTaylortwoPN");
  
  /* this is given in Mpc */    
  injParams.distance = (float)exp(par->logdl);//d_L;
  //	
  /* this should not be zero*/
  //injParams.theta0 = Ltheta0;
  //printf("theta0 = %10.20f\n", injParams.theta0);
  //injParams.phi0 = Lphi0;
  //printf("phi0 = %10.20f\n", injParams.phi0);
  // printf("phi0 ?= %10.20f\n", Lphi02);
  
  injParams.inclination = (float)acos(LdotN);//*****************************************************************************************bornes de l'inclinaison a verifier**************
  //	
  //injParams.polarization = locpolar;
  //printf("polarization = %10.20f\n", injParams.polarization);
  
  injParams.spin1x = (float)n_S[0];
  injParams.spin1y = (float)n_S[1];
  injParams.spin1z = (float)n_S[2];
  //	
  
  //
  
  injParams.spin2x = 0.0;
  injParams.spin2y = 0.0;
  injParams.spin2z = 0.0;
  
  // 4 parameters used after the computation of h+ hx ********************//
  injParams.coa_phase = (float)par->phase;
  injParams.longitude = (float)par->longi;
  injParams.latitude = (float)asin(par->sinlati);
  injParams.polarization = (float)par->alpha;    
  
  ppnParams.deltaT = inversesamplerate;//1.0 / 4096.0;
  
  /*
    printf("f_lower = %10.20f\n", injParams.f_lower);
    printf("f_final = %10.20f\n", injParams.f_final);
    
    printf("mass1 = %10.20f\n", injParams.mass1);
    printf("mass2 = %10.20f\n", injParams.mass2);
    printf("distance = %10.20f\n", injParams.distance);
    printf("inclination = %10.20f\n", injParams.inclination);
    printf("spin1x = %10.20f\n", injParams.spin1x);
    printf("spin1y = %10.20f\n", injParams.spin1y);
    printf("spin1z = %10.20f\n", injParams.spin1z);
    printf("a = %10.20f\n", sqrt(injParams.spin1x*injParams.spin1x+injParams.spin1y*injParams.spin1y+injParams.spin1z*injParams.spin1z));
    printf("ifo = %d\n", ifonr);*/
  
  //  fprintf(stderr, "Lower cut-off frequency used will be %fHz\n", injParams.f_lower);
  
  /* --- now we can call the injection function --- */
  // LALGenerateInspiral( &mystatus, &thewaveform, &injParams, &ppnParams );
  LALGenerateInspiral( &mystatus, waveform, &injParams, &ppnParams );
  if ( mystatus.statusCode )
    {
      fprintf( stderr, "LALSTPNWaveformTest: error generating waveform\n" );
      exit( 1 );
    }
  //printf("test10\n");
  /* --- and finally save in a file --- */
  
  // outputfile = fopen(filename,"w");
  
  //lengthLAL  = thewaveform.phi->data->length;
  lengthLAL  = waveform->phi->data->length;
  //
  //double **wave;
  /*wave = malloc( 3 * sizeof(**wave) );
    for ( i = 0 ; i < 3 ; ++i ) {
    (*wave)[i] = malloc( (lengthLAL+2) * sizeof((*wave)[i]) );
    }*/
  
  //printf("lengthLAL = %d\n",lengthLAL);
  
  
  
  //  dt      = thewaveform.phi->deltaT;
  //	dt      = waveform->phi->deltaT;
  //  phi0    = thewaveform.phi->data->data[0];
  
  
  *l = lengthLAL;
  //		printf("length = %d\n",length);
  
  //wave[1][0] = 0.0;
  //wave[2][0] = 0.0;
  
  /* for(i = 0; i < lengthLAL && i < length; i++) {
     a1  = thewaveform.a->data->data[2*i];
     a2  = thewaveform.a->data->data[2*i+1];
     phi     = thewaveform.phi->data->data[i];// - phi0;
     shift   = thewaveform.shift->data->data[i];*/
  
  for(i = 0; i < lengthLAL && i < length; i++) {
    a1  = waveform->a->data->data[2*i];
    a2  = waveform->a->data->data[2*i+1];
    phi     = waveform->phi->data->data[i];// - phi0;
    shift   = waveform->shift->data->data[i];
    
    /*   fprintf(outputfile,"%e %e %e\n",
	 i*dt,
	 a1*cos(shift)*cos(phi) - a2*sin(shift)*sin(phi),
	 a1*sin(shift)*cos(phi) + a2*cos(shift)*sin(phi));*/
    
    
    
    //wave[0][i+1] = i*dt;
    hplus[i] = a1*cos(shift)*cos(phi) - a2*sin(shift)*sin(phi);
    hcross[i] = a1*sin(shift)*cos(phi) + a2*cos(shift)*sin(phi);
    
  }
  
  
  
  
  
  //	printf("tc = %10.10f\n", ppnParams.tc);
  //	printf("tf = %10.10f\n", waveform->phi->data->length*waveform->phi->deltaT);
  //printf("tf = %d\t*%f\t = %10.10f\n", waveform->a->data->length, waveform->a->deltaT, waveform->a->data->length*waveform->a->deltaT);
  //printf("dt = %f\n", inversesamplerate);
  
  //  fclose(outputfile);
  //  fprintf(stdout,"waveform saved in a file\n" );
  
  //printf("%10.20f\n", wave[0][lengthLAL-1]);
  //printf("%10.20f\n", dt);
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //	free(par);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // return wave;
  
  //  CoherentGW *waveform;
  
  //  waveform = &thewaveform;
  
  //return waveform;
  //LALfreedom(&thewaveform);
  
  //return &thewaveform;
  
  /*  LALSDestroyVectorSequence(&mystatus, &( waveform->a->data ));
      LALSDestroyVector(&mystatus, &( waveform->f->data ));
      LALDDestroyVector(&mystatus, &( waveform->phi->data ));
      LALSDestroyVector(&mystatus, &( waveform->shift->data ));
      
      LALFree( waveform->a );
      LALFree( waveform->f ); 
      LALFree( waveform->phi) ;
      LALFree( waveform->shift );*/
  
}











// Compute the detector response for a given detector (ifonr) and h_+,h_x
//double LALFpFc(CoherentGW *waveform, double *wave, int *l, int length, struct parset *par, int ifonr) {  //MvdS: l not used //Vivien: very true.
double LALFpFc(CoherentGW *waveform, double *wave, int length, struct parset *par, int ifonr) {
  
  static LALStatus stat;     // status structure
  
  memset( &stat, 0, sizeof(LALStatus) );
  
  // CHAR *sourcefile = NULL;   // name of sourcefile 
  // CHAR *respfile = NULL;     // name of respfile 
  // CHAR *infile = NULL;       // name of infile 
  // CHAR *outfile = NULL;      // name of outfile 
  //INT4 seed = 0;             // random number seed 
  //INT4 sec = SEC;            // ouput epoch.gpsSeconds 
  //INT4 nsec = NSEC;          // ouput epoch.gpsNanoSeconds 
  //INT4 npt = NPT;            // number of output points
  //REAL8 dt = DT;             // output sampling interval 
  //REAL4 sigma = SIGMA;       // noise amplitude 
  int i;
  
  
  DetectorResponse detector;   // the detector in question 
  
  LALDetector site;
  
  if(ifonr==0) site = lalCachedDetectors[LALDetectorIndexLHODIFF]; 
  if(ifonr==1) site = lalCachedDetectors[LALDetectorIndexLLODIFF];
  if(ifonr==2) site = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  
  //site = lalCachedDetectors[LALDetectorIndexLLODIFF];
  
  detector.site = &site;
  detector.transfer = NULL;
  detector.ephemerides = NULL;
  
  
  // REAL4 m1, m2, dist, inc, phic; // unconverted parameters 
  
  REAL4TimeSeries signal;        // GW signal 
  // REAL8 time;                    // length of GW signal 
  // CHAR timeCode;                 // code for signal time alignment 
  // CHAR message[MSGLEN];          // warning/info messages 
  
  signal.epoch.gpsSeconds = (INT4)par->tc;
  signal.epoch.gpsNanoSeconds = (INT4) 100*(int)(1000000.0*(par->tc - signal.epoch.gpsSeconds));
  
  waveform->f->epoch = waveform->phi->epoch = waveform->a->epoch = signal.epoch; 
  
  signal.deltaT = waveform->phi->deltaT;
  signal.f0 = 0.0;
  signal.data = NULL;
  //  time = ( time + 2.0 )/signal.deltaT;
  
  //	 printf("signal.epoch.gpsSeconds(1) = %d\n",signal.epoch.gpsSeconds);
  //	 printf("signal.epoch.gpsNanoSeconds(1) = %d\n",signal.epoch.gpsNanoSeconds);
  
  LALSCreateVector( &stat, &( signal.data ), (UINT4)waveform->phi->data->length );
  
  
  LALSimulateCoherentGW( &stat, &signal, waveform, &detector );
  
  for ( i = 0; i < signal.data->length && i < length; i++ ){
    wave[i] = signal.data->data[i]; 
    //   printf("%d\t%e\n", i, wave[i]);
  }
  
  /*********TIME DELAY***********/
  
  //LIGOTimeGPS      gps;
  //SkyPosition      source;
  REAL8            delay;
  //REAL8            difference;
  DetTimeAndASource     det1_and_source;
  /* TwoDetsTimeAndASource dets_and_source; */
  LALPlaceAndGPS        det1_and_gps;
  /* LALPlaceAndGPS        det2_and_gps; */
  
  // source = waveform->position;	  
	  
  // gps = signal.epoch;
  
  det1_and_gps.p_detector = detector.site;
  det1_and_gps.p_gps      = &(signal.epoch);
  
  det1_and_source.p_det_and_time = &det1_and_gps;
  det1_and_source.p_source       = &(waveform->position);
  
  LALTimeDelayFromEarthCenter(&stat, &delay, &det1_and_source);
  
  //	 printf("LALdelay1 = %10.10f\n", -delay);
  //	 printf("ifo = %d\n", ifonr);
  //	 printf("signal.epoch.gpsSeconds(2) = %d\n",signal.epoch.gpsSeconds);
  //	 printf("signal.epoch.gpsNanoSeconds(2) = %d\n",signal.epoch.gpsNanoSeconds);
  
  LALSDestroyVector( &stat, &( signal.data ) );
  
  return -delay;
  
}






void LALfreedom(CoherentGW *waveform) {
  // Free LAL stuff  
  static LALStatus stat;     /* status structure */
  
  memset( &stat, 0, sizeof(LALStatus) );
  
  LALSDestroyVectorSequence(&stat, &( waveform->a->data ));
  LALSDestroyVector(&stat, &( waveform->f->data ));
  LALDDestroyVector(&stat, &( waveform->phi->data ));
  LALSDestroyVector(&stat, &( waveform->shift->data ));
  
  LALFree( waveform->a );
  LALFree( waveform->f ); 
  LALFree( waveform->phi) ;
  LALFree( waveform->shift );
  
}

