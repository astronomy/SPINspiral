/* 
   
   SPINspiral:                parameter estimation on binary inspirals detected by LIGO, including spins of the binary members
   mcmc_parameters.c:         routines to read/write input files, set constants and set true and null parameters
   
   
   Copyright 2007, 2008, 2009 Marc van der Sluys, Vivien Raymond, Christian Roever, Ilya Mandel
   
   
   This file is part of SPINspiral.
   
   SPINspiral is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   SPINspiral is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with SPINspiral.  If not, see <http://www.gnu.org/licenses/>.
   
*/


#include <getopt.h>

#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXMLRead.h>

#include <mcmc.h>





// ****************************************************************************************************************************************************  
void readCommandLineOptions(int argc, char* argv[], struct runPar *run)
{
  int c;
  if(argc > 1) printf("   Parsing %i command-line arguments:\n",argc-1);
  
  
  
  //Set up struct with long (--) options:
  static struct option long_options[] =
    {
      {"injXMLfile",      required_argument, 0,               0},
      {"injXMLnr",        required_argument, 0,               0},
      {0, 0, 0, 0}
    };
  
  
  int option_index = 0;
  while( (c = getopt_long(argc, argv, "i:",long_options, &option_index)) != -1) {
    switch(c) {
      
      
      //Treat (untranslated) long options:
    case 0:
      if(strcmp(long_options[option_index].name,"injXMLfile")==0) {
	run->injXMLfilename=(char*)malloc(strlen(optarg)+1);
	strcpy(run->injXMLfilename,optarg);
	printf("     Reading injection parameters from the XML file %s\n",run->injXMLfilename);
      }
      if(strcmp(long_options[option_index].name,"injXMLnr")==0) {
	run->injXMLnr = atoi(optarg);
	printf("     Using injection %d from the injection XML file\n",run->injXMLnr);
      }
      
      break; //For case 0: long options
      
      
      
      //Treat the short and translated long options:
    case 'i':
      strcpy(run->mainFilename,optarg);
      printf("     Using main input file %s\n",run->mainFilename);
      break;
      
    default:
      //fprintf(stderr,"   Unrecognised option: %d\n",c);  // This notice is already produced by getopt_long()
      fprintf(stderr,USAGE); 
      exit(1);
      break;
      
    } // switch()
  } // while()
  
  
  
  
  // Print any remaining command line arguments (the ones that are not options, i.e. lack a starting '-'):
  if(optind < argc) {
    printf("   SPINspiral - unused command-line arguments: ");
    while(optind < argc) printf ("%s ", argv[optind++]);
    printf("\n");
  }
  
} // End void readCommandLineOptions(argc,argv)
// ****************************************************************************************************************************************************  







// Read the main input file.
// All parameters that are read here should be(come) members of the runvar struct and lose their global status
// ****************************************************************************************************************************************************  
void readMainInputfile(struct runPar *run)
{
  int i;
  char bla[500];
  FILE *fin;
  
  if((fin = fopen(run->mainFilename,"r")) == NULL) {
    printf("   Error reading main input file: %s, aborting.\n\n\n",run->mainFilename);
    exit(1);
  } else {
    printf("   Using main input file: %s.\n",run->mainFilename);
  }
  
  
  //Use and l for floats: %lf, %lg, etc, rather than %f, %g
  
  for(i=1;i<=3;i++) fgets(bla,500,fin);  //Read first 3 lines
  
  //Operation and output:
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin);  sscanf(bla,"%d",&doSNR);
  fgets(bla,500,fin);  sscanf(bla,"%d",&doMCMC);
  fgets(bla,500,fin);  sscanf(bla,"%d",&doMatch);
  fgets(bla,500,fin);  sscanf(bla,"%d",&intscrout);
  fgets(bla,500,fin);  sscanf(bla,"%d",&writeSignal);
  fgets(bla,500,fin);  sscanf(bla,"%d",&printMuch);
  
  
  //Secondary input files:
  fgets(bla,500,fin); fgets(bla,500,fin); fgets(bla,500,fin); //Read the empty and comment lines
  fgets(bla,500,fin); sscanf(bla,"%s",run->mcmcFilename);
  fgets(bla,500,fin); sscanf(bla,"%s",run->dataFilename);
  fgets(bla,500,fin); sscanf(bla,"%s",run->injectionFilename);
  fgets(bla,500,fin); sscanf(bla,"%s",run->parameterFilename);
  fgets(bla,500,fin); sscanf(bla,"%s",run->systemFilename);
  
  fclose(fin);
}  //End of readMainInputfile
// ****************************************************************************************************************************************************  










// Read the MCMC input file.
// All parameters that are read in here should be(come) members of the runvar struct and lose their global status
// ****************************************************************************************************************************************************  
void readMCMCinputfile(struct runPar *run)
{
  int i;
  double tmpdbl;
  char bla[500];
  FILE *fin;
  
  if((fin = fopen(run->mcmcFilename,"r")) == NULL) {
    printf("   Error reading MCMC input file: %s, aborting.\n\n\n",run->mcmcFilename);
    exit(1);
  } else {
    printf("   Using MCMC input file: %s.\n",run->mcmcFilename);
  }
  
  
  //Use and l for floats: %lf, %lg, etc, rather than %f, %g
  
  for(i=1;i<=3;i++) { //Read first 3 lines
    fgets(bla,500,fin);
  }
  
  //Basic settings
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->mcmcWaveform);
  //run->injectionWaveform = run->mcmcWaveform;  //For now
  
  
  if(run->mcmcWaveform==1) {
    printf("   Using Apostolatos, 1.5PN, 12-parameter waveform as the MCMC template.\n");
    run->nMCMCpar=12;
  } else if(run->mcmcWaveform==2) {
    printf("   Using LAL, 3.5PN, 12-parameter waveform as the MCMC template.\n");
    run->nMCMCpar=12;
  } else if(run->mcmcWaveform==3) {
    printf("   Using LAL, 3.5PN, 15-parameter waveform as the MCMC template.\n");
    run->nMCMCpar=15;
  } else {
    printf("   Unknown waveform chosen as MCMC template: %d.   Available waveforms are:\n",run->mcmcWaveform);
    printf("     1: Apostolatos, simple precession, 12 parameters\n");
    printf("     2: LAL, single spin, 12 parameters\n");
    printf("     3: LAL, double spin, 15 parameters\n");
    printf("   Please set mcmcWaveform in %s to one of these values.\n\n",run->mcmcFilename);
    exit(1);
  }
  //run->nInjectPar = run->nMCMCpar;  //For now
  
  
  fgets(bla,500,fin);  sscanf(bla,"%lg",&tmpdbl);    
  iter = (int)tmpdbl;
  fgets(bla,500,fin);  sscanf(bla,"%d",&thinOutput);
  fgets(bla,500,fin);  sscanf(bla,"%d",&thinScreenOutput);
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->MCMCseed);
  fgets(bla,500,fin);  sscanf(bla,"%d",&adapt);
  fgets(bla,500,fin);  sscanf(bla,"%lf",&run->blockfrac);
  
  

  //Correlated update proposals:
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin);  sscanf(bla,"%d",&corrupd);
  fgets(bla,500,fin);  sscanf(bla,"%lf",&run->corrfrac);
  fgets(bla,500,fin);  sscanf(bla,"%lg",&tmpdbl);
  ncorr = (int)tmpdbl;
  fgets(bla,500,fin);  sscanf(bla,"%lf",&run->mataccfr);
  fgets(bla,500,fin);  sscanf(bla,"%d",&prmatrixinfo);
  
  
  //Annealing:
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin);  sscanf(bla,"%lf",&temp0);
  fgets(bla,500,fin);  sscanf(bla,"%lg",&tmpdbl);
  nburn = (int)tmpdbl;
  fgets(bla,500,fin);  sscanf(bla,"%lg",&tmpdbl);
  nburn0 = (int)tmpdbl;
  
  //Parallel tempering:
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin);  sscanf(bla,"%d",&partemp);
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->ntemps);
  fgets(bla,500,fin);  sscanf(bla,"%lf",&tempmax);
  fgets(bla,500,fin);  sscanf(bla,"%d",&savehotchains);
  fgets(bla,500,fin);  sscanf(bla,"%d",&prpartempinfo);
  
  //Manual temperature ladder for parallel tempering:
  fgets(bla,500,fin); fgets(bla,500,fin); //Read the empty and comment line
  for(i=0;i<run->ntemps;i++) fscanf(fin,"%lf",&run->temps[i]);  //Read the array directly, because sscanf cannot be in a loop...
  
  fclose(fin);
}
// ****************************************************************************************************************************************************  
//End of void readMCMCinputfile(struct runPar *run)












// Read the data input file.
// ****************************************************************************************************************************************************  
void readDataInputfile(struct runPar *run, struct interferometer ifo[])
{
  int i=0,j=0;
  double lati,longi,leftArm,rightArm;
  char bla[500], subdir[500];
  FILE *fin;
  
  if((fin = fopen(run->dataFilename,"r")) == NULL) {
    printf("   Error reading data file: %s, aborting.\n\n\n",run->dataFilename);
    exit(1);
  }
  
  
  //Use and l for floats: %lf, %lg, etc, rather than %f, %g
  
  for(j=1;j<=3;j++) fgets(bla,500,fin);  //Read first 3 lines
  fgets(run->datasetName,80,fin);  fgets(bla,500,fin);  //Read name of the data set used, and then the rest of the line
  
  //Detector network:
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->networksize);
  for(i=0;i<run->networksize;i++) fscanf(fin,"%d",&run->selectifos[i]);  //Read the array directly, because sscanf cannot be in a loop...
  fgets(bla,500,fin);  //Read the rest of the line
  
  
  //Data handling:
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin); sscanf(bla,"%d",&downsamplefactor);
  fgets(bla,500,fin); sscanf(bla,"%lf",&run->dataBeforeTc);
  fgets(bla,500,fin); sscanf(bla,"%lf",&run->dataAfterTc);
  fgets(bla,500,fin); sscanf(bla,"%lf",&run->lowFrequencyCut);
  fgets(bla,500,fin); sscanf(bla,"%lf",&run->highFrequencyCut);
  fgets(bla,500,fin); sscanf(bla,"%lf",&run->tukeyWin);
  
  
  //Read input for PSD estimation:
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment lines
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->PSDsegmentNumber);
  fgets(bla,500,fin);  sscanf(bla,"%lf",&run->PSDsegmentLength);
  
  
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment lines
  for(i=0;i<run->networksize;i++){
    fgets(bla,500,fin); fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment lines
    
    fgets(bla,500,fin);  sscanf(bla,"%s",ifo[i].name);
    fgets(bla,500,fin);  sscanf(bla,"%lf",&lati);
    fgets(bla,500,fin);  sscanf(bla,"%lf",&longi);
    fgets(bla,500,fin);  sscanf(bla,"%lf",&rightArm);
    fgets(bla,500,fin);  sscanf(bla,"%lf",&leftArm);
    
    ifo[i].lati      = lati     /180.0*pi;
    ifo[i].longi     = longi    /180.0*pi;
    ifo[i].rightArm  = rightArm /180.0*pi;
    ifo[i].leftArm   = leftArm  /180.0*pi;
    
    fgets(bla,500,fin);  //Read the empty line
    
    fgets(bla,500,fin);  sscanf(bla,"%s",ifo[i].ch1name);
    //fgets(bla,500,fin);  sscanf(bla,"%s",&ifo[i].ch1filepath);
    fgets(bla,500,fin);  sscanf(bla,"%s",subdir);
    sprintf(ifo[i].ch1filepath,"%s%s%s",datadir,"/",subdir);
    fgets(bla,500,fin);  sscanf(bla,"%s",ifo[i].ch1fileprefix);
    fgets(bla,500,fin);  sscanf(bla,"%s",ifo[i].ch1filesuffix);
    fgets(bla,500,fin);  sscanf(bla,"%d",&ifo[i].ch1filesize);
    fgets(bla,500,fin);  sscanf(bla,"%d",&ifo[i].ch1fileoffset);
    fgets(bla,500,fin);  sscanf(bla,"%d",&ifo[i].ch1doubleprecision);
    fgets(bla,500,fin);  sscanf(bla,"%d",&ifo[i].add2channels);
    
    fgets(bla,500,fin);  //Read the empty line
    
    fgets(bla,500,fin);  sscanf(bla,"%ld",&ifo[i].noiseGPSstart);
    fgets(bla,500,fin);  sscanf(bla,"%s",ifo[i].noisechannel);
    //fgets(bla,500,fin);  sscanf(bla,"%s",&ifo[i].noisefilepath);
    fgets(bla,500,fin);  sscanf(bla,"%s",subdir);
    sprintf(ifo[i].noisefilepath,"%s%s%s",datadir,"/",subdir);
    fgets(bla,500,fin);  sscanf(bla,"%s",ifo[i].noisefileprefix);
    fgets(bla,500,fin);  sscanf(bla,"%s",ifo[i].noisefilesuffix);
    fgets(bla,500,fin);  sscanf(bla,"%d",&ifo[i].noisefilesize);
    fgets(bla,500,fin);  sscanf(bla,"%d",&ifo[i].noisefileoffset);
    fgets(bla,500,fin);  sscanf(bla,"%d",&ifo[i].noisedoubleprecision);
    
  }
  fclose(fin);
  
}  //End of readDataInputfile
// ****************************************************************************************************************************************************  





// All parameters that are read in here should be(come) members of the runvar struct
// ****************************************************************************************************************************************************  
void readInjectionInputfile(struct runPar *run)
{
  int i;
  char bla[500];
  FILE *fin;
  
  if((fin = fopen(run->injectionFilename,"r")) == NULL) {
    printf("   Error reading injection input file: %s, aborting.\n\n\n",run->injectionFilename);
    exit(1);
  } else {
    printf("   Using injection input file: %s.\n",run->injectionFilename);
  }
  
  
  //Use and l for floats: %lf, %lg, etc, rather than %f, %g
  
  for(i=1;i<=2;i++) fgets(bla,500,fin);  //Read first 2 lines
  
  //General:
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->injectSignal);
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->injectionWaveform);
  fgets(bla,500,fin);  sscanf(bla,"%lf",&run->injectionSNR);
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->injRanSeed);
  
  //Get the number of injection parameters from the injectionWaveform
  if(run->injectSignal >= 1) {
    if(run->injectionWaveform==1) {
      printf("   Using Apostolatos, 1.5PN, 12-parameter waveform for the software injection.\n");
      run->nInjectPar=12;
    } else if(run->injectionWaveform==2) {
      printf("   Using LAL, 3.5PN, 12-parameter waveform for the software injection.\n");
      run->nInjectPar=12;
    } else if(run->injectionWaveform==3) {
      printf("   Using LAL, 3.5PN, 15-parameter waveform for the software injection.\n");
      run->nInjectPar=15;
    } else {
      printf("   Unknown waveform chosen as MCMC template: %d.   Available waveforms are:\n",run->injectionWaveform);
      printf("     1: Apostolatos, simple precession, 12 parameters\n");
      printf("     2: LAL, single spin, 12 parameters\n");
      printf("     3: LAL, double spin, 15 parameters\n");
      printf("   Please set injectionWaveform in %s to one of these values.\n\n",run->injectionFilename);
      exit(1);
    }
  }
  
  //Parameters:
  for(i=1;i<=5;i++) fgets(bla,500,fin);  //Read empty and comment lines
  
  for(i=0;i<run->nInjectPar;i++) {
    fscanf(fin,"%d %d %lf %d %lf %d %lf %lf",&run->injNumber[i],&run->injID[i],&run->injParValOrig[i],&run->injRanPar[i],&run->injSigma[i],&run->injBoundType[i],&run->injBoundLow[i],&run->injBoundUp[i]);
    fgets(bla,500,fin);  //Read rest of the line
    
    //printf("%d %d %lf %d %lf %d %lf %lf\n",run->injNumber[i],run->injID[i],run->injParValOrig[i],run->injRanPar[i],run->injSigma[i],run->injBoundType[i],run->injBoundLow[i],run->injBoundUp[i]);
    
    
    if(run->injNumber[i] != i+1) {
      printf("   Error reading injection input file %s:  parameter %d has number %d.\n   Aborting...\n\n",run->injectionFilename,i+1,run->injNumber[i]);
      exit(1);
    }
    
    if(run->parDef[run->injID[i]] != 1) {
      printf("\n\n   Error reading injection input file %s, parameter %d:\n     parameter ID %d is not defined.\n   Aborting...\n\n",
	     run->injectionFilename,run->injNumber[i],run->injID[i]);
      exit(1);
    }
    
    run->injRevID[run->injID[i]] = i;  //Reverse parameter ID
    
    
    // Get the desired injection boundaries:
    switch (run->injBoundType[i]) {
    case 1 :
      break;
    case 2 :
      if(run->injBoundLow[i] > 0.0 || run->injBoundUp[i] < 0.0) {
	printf("\n\n   Error reading injection input file %s, parameter %d (%s):\n     for injBoundType = 2, injBoundLow and injBoundUp must be <= 0 and >= 0 respectively.\n   Aborting...\n\n",
	       run->injectionFilename,run->injNumber[i],run->parAbrev[run->injID[i]]);
	exit(1);
      }
      run->injBoundLow[i] = run->injParValOrig[i] + run->injBoundLow[i];
      run->injBoundUp[i]  = run->injParValOrig[i] + run->injBoundUp[i];
      break;
    case 3 :
      if(run->injBoundLow[i] > 1.0 || run->injBoundUp[i] < 1.0) {
	printf("\n\n   Error reading injection input file %s, parameter %d (%s):\n     for injBoundType = 3, injBoundLow and injBoundUp must be <= 1 and >= 1 respectively.\n   Aborting...\n\n",
	       run->injectionFilename,run->injNumber[i],run->parAbrev[run->injID[i]]);
	exit(1);
      }
      run->injBoundLow[i] = run->injParValOrig[i] * run->injBoundLow[i];
      run->injBoundUp[i]  = run->injParValOrig[i] * run->injBoundUp[i];
      break;
    default :
      printf("\n\n   Error reading injection input file %s, parameter %d (%s):\n     %d is not a valid option for injBoundType.\n   Aborting...\n\n",
	     run->injectionFilename,run->injNumber[i],run->parAbrev[run->injID[i]],run->injBoundType[i]);
      exit(1);
    } //End switch
    
    
    
    // Check whether value for injRanPar is valid
    if(run->injRanPar[i] < 0 || run->injRanPar[i] > 2) {
      printf("\n\n   Error reading injection input file %s, parameter %d (%s):\n     %d is not a valid option for injRanPar.\n   Aborting...\n\n",
	     run->injectionFilename,run->injNumber[i],run->parAbrev[run->injID[i]],run->injRanPar[i]);
	exit(1);
    }      
    
    //Check whether the lower boundary < the upper
    if(run->injBoundLow[i] >= run->injBoundUp[i]) {
      printf("\n\n   Error reading injection input file %s, parameter %d (%s):\n     the lower boundary of the prior is larger than or equal to the upper boundary (%lf vs. %lf).\n   Aborting...\n\n",
	     run->injectionFilename,run->injNumber[i],run->parAbrev[run->injID[i]],run->injBoundLow[i],run->injBoundUp[i]);
      exit(1);
    }
    
    //Check whether  lower boundary <= injection value <= upper boundary
    if(run->injParValOrig[i] < run->injBoundLow[i] || run->injParValOrig[i] > run->injBoundUp[i]) {
      printf("\n\n   Error reading injection input file %s, parameter %d (%s):\n     the injection value (%lf) lies outside the prior range (%lf - %lf).\n   Aborting...\n\n",
	     run->injectionFilename,run->injNumber[i],run->parAbrev[run->injID[i]], run->injParValOrig[i], run->injBoundLow[i], run->injBoundUp[i]);
      exit(1);
    }
    
  } //End for
  
  
  
  setRandomInjectionParameters(run);    //Copy the injection parameters from injParValOrig to injParVal, and randomise where wanted
  prior_tc_mean = run->injParVal[2];    //CHECK prior_tc_mean is (still) used everywhere.  This value must be overwritten by the 'best' value in readParameterInputfile() which is called next, in the case of no SW injection
  
  
  fclose(fin);
  
  
  
  //Print injection parameters and prior ranges to screen:
  char StartStr[3][99];
  strcpy(StartStr[0],"Injection value");
  strcpy(StartStr[1],"Random value near injection value");
  strcpy(StartStr[2],"Random value from prior");
  
  printf("\n   Software-injection parameters:\n      Nr: Name:           Injection value:     Obtained:\n");
  for(i=0;i<run->nInjectPar;i++) {
    //printf("      %2d  %-11s     %15.4lf     %15.4lf %15.4lf     %-25s\n",run->injNumber[i],run->parAbrev[run->injID[i]],run->injParVal[i],
    //	   run->injBoundLow[i],run->injBoundUp[i],  StartStr[run->injRanPar[i]]);
    if(run->injRanPar[i]==0) {
      printf("      %2d  %-11s     %15.4lf      Taken from the value set in %s\n",run->injNumber[i],run->parAbrev[run->injID[i]],run->injParVal[i],
	     run->injectionFilename);
    } else if(run->injRanPar[i]==1) {
      printf("      %2d  %-11s     %15.4lf      Drawn randomly from a Gaussian distribution with centre  %lf  and width  %lf\n",run->injNumber[i],
      	     run->parAbrev[run->injID[i]],run->injParVal[i],run->injParValOrig[i],run->injSigma[i]);
    } else if(run->injRanPar[i]==2) {
      printf("      %2d  %-11s     %15.4lf      Drawn randomly from a uniform distribution  %14.4lf - %-14.4lf\n",run->injNumber[i],run->parAbrev[run->injID[i]]
	     ,run->injParVal[i],run->injBoundLow[i],run->injBoundUp[i]);
    }
  }
  printf("\n");
  
  
}  //End of readInjectionInputfile
// ****************************************************************************************************************************************************  













// All parameters that are read in here should be(come) members of the runvar struct
// ****************************************************************************************************************************************************  
void readParameterInputfile(struct runPar *run)
{
  int i,iInj;
  char bla[500];
  FILE *fin;
  
  if((fin = fopen(run->parameterFilename,"r")) == NULL) {
    printf("   Error reading parameter input file: %s, aborting.\n\n\n",run->parameterFilename);
    exit(1);
  } else {
    printf("   Using parameter input file: %s.\n",run->parameterFilename);
  }
  
  
  //Use and l for floats: %lf, %lg, etc, rather than %f, %g
  
  for(i=1;i<=2;i++) fgets(bla,500,fin);  //Read first 2 lines
  
  //Priors:
  fgets(bla,500,fin); fgets(bla,500,fin);  //Read the empty and comment line
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->priorSet);
  fgets(bla,500,fin);  sscanf(bla,"%d",&run->offsetMCMC);
  fgets(bla,500,fin);  sscanf(bla,"%lf",&run->offsetX);
  
  
  
  //Parameters:
  for(i=1;i<=5;i++) fgets(bla,500,fin);  //Read empty and comment lines
  
  for(i=0;i<run->nMCMCpar;i++) {
    fscanf(fin,"%d %d %lf %d %d %lf %d %lf %lf",&run->parNumber[i],&run->parID[i],&run->parBestVal[i],&run->parFix[i],&run->parStartMCMC[i],&run->parSigma[i],&run->priorType[i],&run->priorBoundLow[i],&run->priorBoundUp[i]);
    fgets(bla,500,fin);  //Read rest of the line
    
    //printf("%d:  %d %d %lf %d %lf %d %lf %lf\n",i,run->parNumber[i],run->parID[i],run->parBestVal[i],run->parStartMCMC[i],run->parSigma[i],
    //   run->priorType[i],run->priorBoundLow[i],run->priorBoundUp[i]);
    
    
    if(run->parNumber[i] != i+1) {
      printf("   Error reading parameter input file %s:  parameter %d has number %d.\n   Aborting...\n\n",run->parameterFilename,i+1,run->parNumber[i]);
      exit(1);
    }
    
    if(run->parDef[run->parID[i]] != 1) {
      printf("\n\n   Error reading parameter input file %s, parameter %d:\n     parameter ID %d is not defined.\n   Aborting...\n\n",
	     run->injectionFilename,run->parNumber[i],run->parID[i]);
      exit(1);
    }
    
    run->parRevID[run->parID[i]] = i;  //Reverse parameter ID
    
    
    
    // Get the desired boundary conditions:
    switch (run->priorType[i]) {
    case 11 :
      break;
    case 12 :
      if(run->priorBoundLow[i] > 0.0 || run->priorBoundUp[i] < 0.0) {
	printf("\n\n   Error reading parameter input file %s, parameter %d (%s):\n     for priorType = 12, priorBoundLow and priorBoundUp must be <= 0 and >= 0 respectively.\n   Aborting...\n\n",
	       run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]]);
	exit(1);
      }
      run->priorBoundLow[i] = run->parBestVal[i] + run->priorBoundLow[i];
      run->priorBoundUp[i]  = run->parBestVal[i] + run->priorBoundUp[i];
      break;
    case 13 :
      if(run->priorBoundLow[i] > 1.0 || run->priorBoundUp[i] < 1.0) {
	printf("\n\n   Error reading parameter input file %s, parameter %d (%s):\n     for priorType = 13, priorBoundLow and priorBoundUp must be <= 1 and >= 1 respectively.\n   Aborting...\n\n",
	       run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]]);
	exit(1);
      }
      run->priorBoundLow[i] = run->parBestVal[i] * run->priorBoundLow[i];
      run->priorBoundUp[i]  = run->parBestVal[i] * run->priorBoundUp[i];
      break;
    case 21 : 
      run->priorBoundLow[i] = 0.0;
      run->priorBoundUp[i]  = tpi;
      break;
    case 22 :
      run->priorBoundLow[i] = 0.0;
      run->priorBoundUp[i]  = pi;
      break;
    default :
      printf("\n\n   Error reading parameter input file %s, parameter %d (%s):\n     %d is not a valid option for priorType.\n   Aborting...\n\n",
	     run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]],run->priorType[i]);
      exit(1);
    } //End switch
    
    
    // Check whether value for fix is valid
    if(run->parFix[i] < 0 || run->parFix[i] > 2) {
      printf("\n\n   Error reading parameter input file %s, parameter %d (%s):\n     %d is not a valid option for parFix.\n   Aborting...\n\n",
	     run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]],run->parFix[i]);
	exit(1);
    }      
    
    // Check whether value for start is valid
    if(run->parStartMCMC[i] < 1 || run->parStartMCMC[i] > 5) {
      printf("\n\n   Error reading parameter input file %s, parameter %d (%s):\n     %d is not a valid option for parStartMCMC.\n   Aborting...\n\n",
	     run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]],run->parStartMCMC[i]);
	exit(1);
    }      
    
    //Check whether the lower prior boundary < the upper
    if(run->priorBoundLow[i] >= run->priorBoundUp[i]) {
      printf("\n\n   Error reading parameter input file %s, parameter %d (%s):\n     the lower boundary of the prior is larger than or equal to the upper boundary (%lf vs. %lf).\n   Aborting...\n\n",
	     run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]],run->priorBoundLow[i],run->priorBoundUp[i]);
      exit(1);
    }
    
    //Check whether  lower prior boundary <= best value <= upper boundary
    if(run->parBestVal[i] < run->priorBoundLow[i] || run->parBestVal[i] > run->priorBoundUp[i]) {
      printf("\n\n   Error reading parameter input file %s, parameter %d (%s):\n     the best value (%lf) lies outside the prior range (%lf - %lf).\n   Aborting...\n\n",
	     run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]], run->parBestVal[i], run->priorBoundLow[i], run->priorBoundUp[i]);
      exit(1);
    }
    
    //Check whether  lower prior boundary <= INJECTION value <= upper boundary
    iInj = run->injRevID[run->parID[i]];  //Get the index of this parameter in the injection set.  -1 if not available.
    if(iInj >= 0) {
      if(run->injParVal[iInj] < run->priorBoundLow[i] || run->injParVal[iInj] > run->priorBoundUp[i]) {
	printf("\n\n   Error reading parameter input file %s, parameter %d (%s):\n     the injection value (%lf) lies outside the prior range (%lf - %lf).\n   Aborting...\n\n",
	       run->parameterFilename, run->parNumber[i], run->parAbrev[run->parID[i]], run->injParVal[iInj], run->priorBoundLow[i], run->priorBoundUp[i]);
	exit(1);
      }
    } else {
      if(run->injectSignal != 0) printf("   Warning:  MCMC parameter %i (%s) does not occur in the injection template;  I cannot verify whether the injection value lies within the prior range.\n",
					run->parNumber[i],run->parAbrev[run->parID[i]]);
    }
    
  } //End for (i)
  
  
  if(run->injectSignal<=0) {
    prior_tc_mean = run->parBestVal[2];       //CHECK prior_tc_mean is (still) used everywhere.  This value overwrites the injection value from readInjectionInputfile() called earlier
    for(i=0;i<run->nMCMCpar;i++) run->injParVal[i] = run->parBestVal[i];   //CHECK Needed to avoid SegFault in the case of t_c
  }

  fclose(fin);
  
  
  
  //Print MCMC parameters and prior ranges to screen:
  char FixStr[3][99];
  strcpy(FixStr[0],"No, let it free");
  strcpy(FixStr[1],"Yes, to best value");
  strcpy(FixStr[2],"Yes, to injection");
  
  char StartStr[6][99];
  strcpy(StartStr[1],"From best value");
  strcpy(StartStr[2],"Randomly from Gaussian around best value");
  strcpy(StartStr[3],"From injection value");
  strcpy(StartStr[4],"Randomly from Gaussian around injection");
  strcpy(StartStr[5],"Randomly from prior");
  
  
  //Write parameter choice to screen:
  printf("\n   MCMC parameters:\n      Nr: Name:                Best value:     Prior:     min:            max:    Fix parameter?        Start chain:\n");
  for(i=0;i<run->nMCMCpar;i++) {
    printf("      %2d  %-11s     %15.4lf     %15.4lf %15.4lf     %-20s  %-45s\n",run->parNumber[i],run->parAbrev[run->parID[i]],run->parBestVal[i],
	   run->priorBoundLow[i],run->priorBoundUp[i],  FixStr[run->parFix[i]],StartStr[run->parStartMCMC[i]]);
  }
  printf("\n");
  
}  //End of readParameterInputfile
// ****************************************************************************************************************************************************  






//Read the input file for system (system-dependent) variables, e.g. mcmc.input.system
// ****************************************************************************************************************************************************  
void readSystemInputfile(struct runPar *run)
{
  int i;
  char bla[500];
  FILE *fin;
  
  if((fin = fopen(run->systemFilename,"r")) == NULL) {
    printf("   Error reading system file: %s, aborting.\n\n\n",run->systemFilename);
    exit(1);
  } else {
    printf("   Using system input file: %s.\n",run->parameterFilename);
  }
  
  //Use and l for floats: %lf, %lg, etc, rather than %f, %g
  for(i=1;i<=3;i++) { //Read first 3 lines
    fgets(bla,500,fin);
  }  

  //Data directory:
  fscanf(fin, "%s",datadir);
  
  fclose(fin);
}  //End of readSystemInputfile
// ****************************************************************************************************************************************************  





// ****************************************************************************************************************************************************  
void readInjectionXML(struct runPar *run)
{
  int NumInj=0,i=0,j=0;
  SimInspiralTable *injTable = NULL;
  double Sx=0.0,Sy=0.0,Sz=0.0,Sxy=0.0,S=0.0;
  
  printf("  Reading injection XML file %s, injection %d.\n",run->injXMLfilename,run->injXMLnr);
  NumInj = SimInspiralTableFromLIGOLw(&injTable,run->injXMLfilename,0,0);
  if(NumInj <= 0) {
    //fprintf(stderr,"\n  Error reading XML file %s.\n",run->injXMLfilename);
    fprintf(stderr,"   Aborting.\n\n");
    exit(1);
  }
  if(run->injXMLnr >= NumInj) {
    fprintf(stderr,"\n\n  ERROR: requested injection number %d larger than number of injections (%d) in XML file %s.\n  Aborting.\n\n",run->injXMLnr,NumInj-1,run->injXMLfilename);
    exit(1);
  }
  
  j=0;
  while(j<run->injXMLnr) {j++; injTable = injTable->next;}  // Select injection

  i=0;
  for(i=0;i<run->nInjectPar;i++) {
    
    //Time:
    if(run->injID[i] == 11) run->injParVal[i] = injTable->geocent_end_time.gpsSeconds+injTable->geocent_end_time.gpsNanoSeconds*1.0e-9;  // t_c in geocentre
    
    
    //Distance:
    if(run->injID[i] == 21) run->injParVal[i] = pow(injTable->distance,3.0);  // d_L^3
    if(run->injID[i] == 22) run->injParVal[i] = log(injTable->distance);      // log(d_L/Mpc)
    
    
    //Sky location:
    if(run->injID[i] == 31) run->injParVal[i] = injTable->longitude;          // RA
    if(run->injID[i] == 32) run->injParVal[i] = sin(injTable->latitude);      // sin Dec
    
    
    //Orbital/GW phase:
    if(run->injID[i] == 41) run->injParVal[i] = injTable->coa_phase;          // \phi_c
    
    
    //Binary orientation:
    if(run->injID[i] == 51) run->injParVal[i] = cos(injTable->inclination);   // cos(i)
    if(run->injID[i] == 52) run->injParVal[i] = injTable->polarization;       // \psi
    if(run->injID[i] == 53) {
      fprintf(stderr,"\n  *** You're reading injection data from an XML file, while using the Apostolatos J_0.  This may lead to unexpected and unwanted results ***\n");
      //exit(1);
      run->injParVal[i] = sin(injTable->inclination);   // sin(\theta_J0), should probably not be used
    }
    if(run->injID[i] == 54) {
      fprintf(stderr,"\n  *** You're reading injection data from an XML file, while using the Apostolatos J_0.  This may lead to unexpected and unwanted results ***\n");
      //exit(1);
      run->injParVal[i] = injTable->polarization;       // \phi_J0, should probably not be used
    }
    
    
    //Mass:
    if(run->injID[i] == 61) run->injParVal[i] = injTable->mchirp;             // Mc
    if(run->injID[i] == 62) run->injParVal[i] = injTable->eta;                // \eta
    if(run->injID[i] == 63) run->injParVal[i] = injTable->mass1;              // M1
    if(run->injID[i] == 64) run->injParVal[i] = injTable->mass2;              // M2
    
    
    //Spin 1:
    Sx = injTable->spin1x;
    Sy = injTable->spin1y;
    Sz = injTable->spin1z;
    Sxy = sqrt(Sx*Sx+Sy*Sy);
    S = sqrt(Sxy*Sxy+Sz*Sz);
    if(run->injID[i] == 71) run->injParVal[i] = S;                            // Spin1 magnitude (a_1)
    if(run->injID[i] == 72) {
      run->injParVal[i] = Sz/S;                                               // cos(\theta_1): angle between S1 and L (at t_c(?)))
      if(fabs(S)<1.0e-9) run->injParVal[i] = 1.0;                             // In case S1=0
    }
    if(run->injID[i] == 73) run->injParVal[i] = atan2(Sy,Sx);                 // \phi_1: phase of spin in orbital plane (at t_c(?))
    
    if(run->injID[i] == 75) run->injParVal[i] = Sx;                           // Spin1_x magnitude
    if(run->injID[i] == 76) run->injParVal[i] = Sy;                           // Spin1_y magnitude
    if(run->injID[i] == 77) run->injParVal[i] = Sz;                           // Spin1_z magnitude
    
    
    //Spin 2:
    Sx = injTable->spin2x;
    Sy = injTable->spin2y;
    Sz = injTable->spin2z;
    Sxy = sqrt(Sx*Sx+Sy*Sy);
    S = sqrt(Sxy*Sxy+Sz*Sz);
    if(run->injID[i] == 81) run->injParVal[i] = S;                            // Spin2 magnitude (a_2)
    if(run->injID[i] == 82) {
      run->injParVal[i] = Sz/S;                                               // cos(\theta_2): angle between S2 and L (at t_c(?)))
      if(fabs(S)<1.0e-9) run->injParVal[i] = 1.0;                             // In case S2=0
    }
    if(run->injID[i] == 83) run->injParVal[i] = atan2(Sy,Sx);                 // \phi_2: phase of spin in orbital plane (at t_c(?))
    
    if(run->injID[i] == 85) run->injParVal[i] = Sx;                           // Spin2_x magnitude
    if(run->injID[i] == 86) run->injParVal[i] = Sy;                           // Spin2_y magnitude
    if(run->injID[i] == 87) run->injParVal[i] = Sz;                           // Spin2_z magnitude
    
    
    //Merger, ringdown, ...:
    //if(run->injID[i] == 91) run->injParVal[i] = injTable->;                   // 
    
    printf("      %2d  %-11s     %15.4lf\n",run->injNumber[i],run->parAbrev[run->injID[i]],run->injParVal[i]);
    
  } // i (injectPar)
  
  
  run->lowFrequencyCutInj = injTable->f_lower;  // May be 0.0!

  printf("\n");
  
} // End void readInjectionXML()
// ****************************************************************************************************************************************************  











// ****************************************************************************************************************************************************  
void setRandomInjectionParameters(struct runPar *run)  
// Get random values for the injection parameters.
{
  int i=0;
  gsl_rng *ran;
  double ranGauss = 0.0, ranUnif=0.0, db=0.0;
  ran = gsl_rng_alloc(gsl_rng_mt19937);  // GSL random-number seed
  
  // Manually select a random seed, *** USE ONLY FOR TESTING ***
  if(1==2 && run->injRanSeed == 0) {
    printf("\n  *** SELECTING RANDOM SEED ***  This should only be done while testing!!! setRandomInjectionParameters() \n\n");
    run->injRanSeed = 0;
    setseed(&run->injRanSeed);
    printf("  Seed: %d\n", run->injRanSeed);
  }
  
  gsl_rng_set(ran, run->injRanSeed);     // Set seed
  
  for(i=0;i<run->nInjectPar;i++) {
    ranGauss = gsl_ran_gaussian(ran,run->injSigma[i]);                                    //Make sure you always draw the same number of random variables
    ranUnif = gsl_rng_uniform(ran);                                                       //Make sure you always draw the same number of random variables
    if(run->injRanPar[i]==0) {                  
      run->injParVal[i] = run->injParValOrig[i];                                          //Keep the suggested value
    } else if(run->injRanPar[i]==1) {                                                     
      run->injParVal[i] = run->injParValOrig[i] + ranGauss;                               //Draw random number from Gaussian with width ranGauss
      run->injParVal[i] = max(run->injParVal[i],run->injBoundLow[i]);                     //Stick to the boundary, rather than redrawing to keep number of random numbers constant
      run->injParVal[i] = min(run->injParVal[i],run->injBoundUp[i]);
    } else if(run->injRanPar[i]==2) {
      db = run->injBoundUp[i]-run->injBoundLow[i];                                        //Width of the range
      run->injParVal[i] = run->injBoundLow[i] + ranUnif*db;                               //Draw random number from uniform range with width db
    }
  }
  
  gsl_rng_free(ran);
}
// ****************************************************************************************************************************************************  





// ****************************************************************************************************************************************************  
void setParameterNames(struct runPar * run)
{
  //Set 01: time
  strcpy(run->parAbrev[11], "t_c");
  strcpy(run->parAbrv[11], "t_c");
  run->parDef[11] = 1;
  strcpy(run->parAbrev[12], "t_40");
  strcpy(run->parAbrv[12], "t_40");
  run->parDef[12] = 1;
  
  //Set 02: distance
  strcpy(run->parAbrev[21], "d^3");
  strcpy(run->parAbrv[21], "d^3");
  run->parDef[21] = 1;
  strcpy(run->parAbrev[22], "log(d)");
  strcpy(run->parAbrv[22], "logD");
  run->parDef[22] = 1;
  
  //Set 03: sky position
  strcpy(run->parAbrev[31], "R.A.");
  strcpy(run->parAbrv[31], "RA");
  run->parDef[31] = 1;
  strcpy(run->parAbrev[32], "sin(dec)");
  strcpy(run->parAbrv[32], "sdec");
  run->parDef[32] = 1;
  
  //Set 04: phase
  strcpy(run->parAbrev[41], "phi_orb");
  strcpy(run->parAbrv[41], "phio");
  run->parDef[41] = 1;
  
  //Set 05: orientation
  strcpy(run->parAbrev[51], "cos(i)");
  strcpy(run->parAbrv[51], "cosi");
  run->parDef[51] = 1;
  strcpy(run->parAbrev[52], "psi");
  strcpy(run->parAbrv[52], "psi");
  run->parDef[52] = 1;
  strcpy(run->parAbrev[53], "sin th_J0");
  strcpy(run->parAbrv[53], "thJ0");
  run->parDef[53] = 1;
  strcpy(run->parAbrev[54], "phi_J0");
  strcpy(run->parAbrv[54], "phJ0");
  run->parDef[54] = 1;
  
  //Set 06: mass
  strcpy(run->parAbrev[61], "Mc");
  strcpy(run->parAbrv[61], "Mc");
  run->parDef[61] = 1;
  strcpy(run->parAbrev[62], "eta");
  strcpy(run->parAbrv[62], "eta");
  run->parDef[62] = 1;
  strcpy(run->parAbrev[63], "M1");
  strcpy(run->parAbrv[63], "M1");
  run->parDef[63] = 1;
  strcpy(run->parAbrev[64], "M2");
  strcpy(run->parAbrv[64], "M2");
  run->parDef[64] = 1;
  
  //Set 07: spin1
  strcpy(run->parAbrev[71], "a_spin1");
  strcpy(run->parAbrv[71], "asp1");
  run->parDef[71] = 1;
  strcpy(run->parAbrev[72], "cs th_sp1");
  strcpy(run->parAbrv[72], "ths1");
  run->parDef[72] = 1;
  strcpy(run->parAbrev[73], "phi_spin1");
  strcpy(run->parAbrv[73], "phs1");
  run->parDef[73] = 1;
  
  strcpy(run->parAbrev[75], "S1_x");
  strcpy(run->parAbrv[75], "S1x");
  run->parDef[75] = 1;
  strcpy(run->parAbrev[76], "S1_y");
  strcpy(run->parAbrv[76], "s1y");
  run->parDef[76] = 1;
  strcpy(run->parAbrev[77], "S1_z");
  strcpy(run->parAbrv[77], "S1z");
  run->parDef[77] = 1;
  
  //Set 08: spin2
  strcpy(run->parAbrev[81], "a_spin2");
  strcpy(run->parAbrv[81], "asp2");
  run->parDef[81] = 1;
  strcpy(run->parAbrev[82], "cs th_sp2");
  strcpy(run->parAbrv[82], "ths2");
  run->parDef[82] = 1;
  strcpy(run->parAbrev[83], "phi_spin2");
  strcpy(run->parAbrv[83], "phs2");
  run->parDef[83] = 1;

  strcpy(run->parAbrev[85], "S2_x");
  strcpy(run->parAbrv[85], "S2x");
  run->parDef[85] = 1;
  strcpy(run->parAbrev[86], "S2_y");
  strcpy(run->parAbrv[86], "s2y");
  run->parDef[86] = 1;
  strcpy(run->parAbrev[87], "S2_z");
  strcpy(run->parAbrv[87], "S2z");
  run->parDef[87] = 1;
  
  
  //Set 09: merger, ringdown
  //strcpy(run->parAbrev[], "");
  //run->parDef[] = 1;
  
  //Set:
  //strcpy(run->parAbrev[], "");
  //run->parDef[] = 1;
  
}
// ****************************************************************************************************************************************************  






// Set the global variables.
// Many of these are now in the input file or unused.
// This routine should eventually contain mathematical and (astro)physical constants only
// ****************************************************************************************************************************************************  
void setconstants()
{
  tempi = 0; //A global variable that determines the current chain (temperature) in the temperature ladder
  
  // Mathematical constants:
  pi   = 3.141592653589793;   // pi
  tpi  = 6.283185307179586;   // 2 pi
  mtpi = 6.283185307179586e6; // Large multiple of 2 pi (2 megapi)
  
  // Define some physical constants:
  G    = 6.67259e-11;         // 6.674215e-11; */ /* gravity constant (SI)
  c    = 299792458.0;         // speed of light (m/s)
  
  Ms   = 1.9889194662e30;     // solar mass (kg)
  Mpc  = 3.08568025e22;       // metres in a Mpc  (LAL: 3.0856775807e22)
  Mpcs = 1.029272137e14;      // seconds in a Mpc  (Mpc/c)
}
// ****************************************************************************************************************************************************  








// ****************************************************************************************************************************************************  
void getInjectionParameters(struct parset *par, int nInjectionPar, double *injParVal)  //Set the parameters to the 'injection values'
{
  int i=0;
  for(i=0;i<nInjectionPar;i++) {
    par->par[i]      = injParVal[i];
  }
  
  //These should all disappear:
  par->mc       = injParVal[0];                    // Chirp mass
  par->eta      = injParVal[1];                    // mass ratio
  par->tc       = injParVal[2];                    // coalescence time
  par->longi    = injParVal[6];                    // RA
  par->sinlati  = injParVal[7];                    // sin latitude
  par->sinthJ0  = injParVal[9];                    // sin Theta_J0 ~ latitude, pi/2 = NP
  par->phiJ0    = injParVal[10];                   // Phi_J0 ~ azimuthal
  
  par->loctc    = NULL;
  par->localti  = NULL;
  par->locazi   = NULL;
  par->locpolar = NULL;
}
// ****************************************************************************************************************************************************  


// ****************************************************************************************************************************************************  
void getStartParameters(struct parset *par, struct runPar run)  //Set the parameters for the 12-parameter spinning template to the starting values for the MCMC chain
{
  
  int i=0;
  for(i=0;i<run.nMCMCpar;i++) {
    par->par[i]      = run.parBestVal[i];
  }
  
  par->loctc    = NULL;
  par->localti  = NULL;
  par->locazi   = NULL;
  par->locpolar = NULL;
}


//Allocate memory for the vectors in the struct parset
void allocparset(struct parset *par, int networksize)
{
  par->loctc    = NULL;
  par->localti  = NULL;
  par->locazi   = NULL;
  par->locpolar = NULL;
  
  par->loctc    = (double*)calloc(networksize,sizeof(double));
  par->localti  = (double*)calloc(networksize,sizeof(double));
  par->locazi   = (double*)calloc(networksize,sizeof(double));
  par->locpolar = (double*)calloc(networksize,sizeof(double));
}
// ****************************************************************************************************************************************************  


//Deallocate the vectors in the struct parset
// ****************************************************************************************************************************************************  
void freeparset(struct parset *par)
{
  free(par->loctc);         par->loctc        = NULL;
  free(par->localti);       par->localti      = NULL;
  free(par->locazi);        par->locazi       = NULL;
  free(par->locpolar);      par->locpolar     = NULL;
}
// ****************************************************************************************************************************************************  







