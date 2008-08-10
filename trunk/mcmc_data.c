#include <mcmc.h>



// *** Routines to handle IFOs and data I/O ***





// *** Routines that handle IFOs ***

void set_ifo_data(struct runpar run, struct interferometer ifo[])
//Set all the data for all IFOs that may be used
{
  int i=0,numberofdatasets = 5;
  //Description of the data sets below
  char datadescriptions[10][99];
  sprintf(datadescriptions[1],"Gaussian, stationary noise (GPS ~700006000)");
  sprintf(datadescriptions[2],"clean S5 data (GPS ~846226044)");
  sprintf(datadescriptions[3],"playground trigger data (GPS ~845348295)");
  sprintf(datadescriptions[4],"glitchy data (GPS ~846471090)");
  sprintf(datadescriptions[5],"NINJA");
  //sprintf(datadescriptions[],"");
  
  //run.selectdata = max(min(run.selectdata,numberofdatasets),1);
  if(run.selectdata < 1 || run.selectdata > numberofdatasets) {
    printf("\n\n   Unknown dataset %d selected.  Please set the parameter  selectdata  in  %s  to a value between 1 and %d: \n",run.selectdata, run.infilename,numberofdatasets);
    for(i=1;i<=numberofdatasets;i++) {
      printf("     %3d:  %s.\n",i,datadescriptions[i]);
    }
    printf("\n\n");
    exit(0);
  }
  
  // Print selected stretch of data:
  printf("   Data used: %s.\n",datadescriptions[run.selectdata]);
  
  
  
  // HANFORD, H1:
  sprintf(ifo[0].name, "Hanford");
  ifo[0].lati     = (  46.45/180.0)*pi;
  ifo[0].longi    = (-119.41/180.0)*pi;
  ifo[0].rightarm = (  36.80/180.0)*pi;
  ifo[0].leftarm  = ( 126.80/180.0)*pi;
  ifo[0].radius_eqt  = 6378137.0;       /* WGS 84 */
  ifo[0].radius_pole = 6356752.314;       
  ifo[0].lowCut  =   lowfrequencycut;  //The other two detectors may use this low,highCut
  ifo[0].highCut =   highfrequencycut;  //Define lower and upper limits of overlap integral
  ifo[0].before_tc = databeforetc;   //The other two detectors may use this before,after_tc
  ifo[0].after_tc = dataaftertc;   //Define data segment: [t_c-before_tc, t_c+after_tc]
  
  
  if(run.selectdata == 1) {
    // Gaussian, stationary noise
    sprintf(ifo[0].ch1name,       "H1:STRAIN"); 
    sprintf(ifo[0].ch1filepath,   datadir);
    sprintf(ifo[0].ch1fileprefix, "HL-SIM-");
    sprintf(ifo[0].ch1filesuffix, "-6000.gwf");
    ifo[0].ch1filesize   = 6000; 
    ifo[0].ch1fileoffset = 4000;  //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[0].ch1doubleprecision = 0;
    ifo[0].add2channels    = 0;  //0, unless you want to read a signal from file
    
    //This seems to be needed for synthetic data only...
    sprintf(ifo[0].ch2name,       "H1:STRAIN_INSP_INJ_ONLY"); 
    sprintf(ifo[0].ch2filepath,   datadir);
    sprintf(ifo[0].ch2fileprefix, "HL-SIM-");
    sprintf(ifo[0].ch2filesuffix, "-6000.gwf");
    ifo[0].ch2filesize   = 6000; 
    ifo[0].ch2fileoffset = 4000;    //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[0].ch2doubleprecision = 0;
    
    ifo[0].noiseGPSstart   = 700006000;
    sprintf(ifo[0].noisechannel,    "H1:STRAIN");
    sprintf(ifo[0].noisefilepath,   datadir);
    sprintf(ifo[0].noisefileprefix, "HL-SIM-");
    sprintf(ifo[0].noisefilesuffix, "-6000.gwf");
    ifo[0].noisefilesize   = 6000; 
    ifo[0].noisefileoffset = 4000;   //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[0].noisedoubleprecision = 0;
  }
  
  
  if(run.selectdata == 2) {
    // Clean S5 data 2
    sprintf(ifo[0].ch1name,       "H1:LSC-STRAIN");
    sprintf(ifo[0].ch1filepath,   datadir);
    sprintf(ifo[0].ch1fileprefix, "H-H1_RDS_C03_L2-");
    sprintf(ifo[0].ch1filesuffix, "-384.gwf");
    ifo[0].ch1filesize   = 384;
    ifo[0].ch1fileoffset = 242;  //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[0].ch1doubleprecision = 0;
    ifo[0].add2channels    = 0;  //0, unless you want to read a signal from file
    
    ifo[0].noiseGPSstart   = 846226044;
    sprintf(ifo[0].noisechannel,    "H1:LSC-STRAIN");
    sprintf(ifo[0].noisefilepath,   datadir);
    sprintf(ifo[0].noisefileprefix, "H-H1_RDS_C03_L2-");
    sprintf(ifo[0].noisefilesuffix, "-384.gwf");
    ifo[0].noisefilesize   = 384;
    ifo[0].noisefileoffset = 242;  //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[0].noisedoubleprecision = 0;
  }
  
  if(run.selectdata == 3) {
    // Playground trigger 845348295
    sprintf(ifo[0].ch1name,       "H1:LSC-STRAIN");
    sprintf(ifo[0].ch1filepath,   datadir);
    sprintf(ifo[0].ch1fileprefix, "H-H1_RDS_C03_L2-");
    sprintf(ifo[0].ch1filesuffix, "-128.gwf");
    ifo[0].ch1filesize   = 128;
    ifo[0].ch1fileoffset = 53;  //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[0].ch1doubleprecision = 0;
    ifo[0].add2channels    = 0;  //0, unless you want to read a signal from file
    
    ifo[0].noiseGPSstart   = 845348534;
    sprintf(ifo[0].noisechannel,    "H1:LSC-STRAIN");
    sprintf(ifo[0].noisefilepath,   datadir);
    sprintf(ifo[0].noisefileprefix, "H-H1_RDS_C03_L2-");
    sprintf(ifo[0].noisefilesuffix, "-128.gwf");
    ifo[0].noisefilesize   = 128;
    ifo[0].noisefileoffset = 53;  //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[0].noisedoubleprecision = 0;
  }
  
  if(run.selectdata == 4) {
    // glitchy data 
    sprintf(ifo[0].ch1name,       "H1:LSC-STRAIN");
    sprintf(ifo[0].ch1filepath,   datadir);
    sprintf(ifo[0].ch1fileprefix, "H-H1_RDS_C03_L2-");
    sprintf(ifo[0].ch1filesuffix, "-128.gwf");
    ifo[0].ch1filesize   = 128;
    ifo[0].ch1fileoffset = 49;
    ifo[0].ch1doubleprecision = 0;
    ifo[0].add2channels    = 0;  //0, unless you want to read a signal from file
    
    ifo[0].noiseGPSstart   = 846471090;
    sprintf(ifo[0].noisechannel,    "H1:LSC-STRAIN");
    sprintf(ifo[0].noisefilepath,   datadir);
    sprintf(ifo[0].noisefileprefix, "H-H1_RDS_C03_L2-");
    sprintf(ifo[0].noisefilesuffix, "-128.gwf");
    ifo[0].noisefilesize   = 128;
    ifo[0].noisefileoffset = 49;  //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[0].noisedoubleprecision = 0;
  }


  if(run.selectdata == 5) {
    // NINJA data set
    
    sprintf(ifo[0].ch1name,       "H1:STRAIN"); 
    sprintf(ifo[0].ch1filepath,   datadir);
    sprintf(ifo[0].ch1fileprefix, "H1-NINJADATA-");
    sprintf(ifo[0].ch1filesuffix, "-1024.gwf");
    ifo[0].ch1filesize   = 1024; 
    ifo[0].ch1fileoffset = 743;  //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[0].ch1doubleprecision = 0;
    ifo[0].add2channels    = 0;  //0, unless you want to read a signal from file
    
    ifo[0].noiseGPSstart   = 894457315; //894378400;  //894457315;
    sprintf(ifo[0].noisechannel,    "H1:STRAIN");
    sprintf(ifo[0].noisefilepath,   datadir);
    sprintf(ifo[0].noisefileprefix, "H1-NINJADATA-");
    sprintf(ifo[0].noisefilesuffix, "-1024.gwf");
    ifo[0].noisefilesize   = 1024; 
    ifo[0].noisefileoffset = 743;   //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[0].noisedoubleprecision = 0;
  }
  
  

  
  
  
  
  // LIVINGSTON, L1
  sprintf(ifo[1].name, "Livingston");
  ifo[1].lati     = (  30.56/180.0)*pi;
  ifo[1].longi    = ( -90.77/180.0)*pi;
  ifo[1].rightarm = ( 108.00/180.0)*pi;
  ifo[1].leftarm  = ( 198.00/180.0)*pi;
  ifo[1].radius_eqt  = 6378137.0;       /* WGS 84 */
  ifo[1].radius_pole = 6356752.314;       
  ifo[1].lowCut  =  ifo[0].lowCut;
  ifo[1].highCut = ifo[0].highCut;
  ifo[1].before_tc = ifo[0].before_tc;
  ifo[1].after_tc = ifo[0].after_tc;
  
  if(run.selectdata == 1) {
    // Gaussian, stationary noise
    sprintf(ifo[1].ch1name,       "L1:STRAIN");
    sprintf(ifo[1].ch1filepath,   datadir);
    sprintf(ifo[1].ch1fileprefix, "HL-SIM-");
    sprintf(ifo[1].ch1filesuffix, "-6000.gwf");
    ifo[1].ch1filesize   = 6000;
    ifo[1].ch1fileoffset = 4000;  //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[1].ch1doubleprecision = 0;
    ifo[1].add2channels    = 0;  //0, unless you want to read a signal from file
    
    //This seems to be needed for synthetic data only...
    sprintf(ifo[1].ch2name,       "L1:STRAIN_INSP_INJ_ONLY");
    sprintf(ifo[1].ch2filepath,   datadir);
    sprintf(ifo[1].ch2fileprefix, "HL-SIM-");
    sprintf(ifo[1].ch2filesuffix, "-6000.gwf");
    ifo[1].ch2filesize   = 6000;
    ifo[1].ch2fileoffset = 4000;  //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[1].ch2doubleprecision = 0;
    ifo[1].noiseGPSstart   = 700007000;
  
    sprintf(ifo[1].noisechannel,    "L1:STRAIN");
    sprintf(ifo[1].noisefilepath,   datadir);
    sprintf(ifo[1].noisefileprefix, "HL-SIM-");
    sprintf(ifo[1].noisefilesuffix, "-6000.gwf");
    ifo[1].noisefilesize   = 6000;
    ifo[1].noisefileoffset = 4000;  //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[1].noisedoubleprecision = 0;
  }
    
  if(run.selectdata == 2) {
    //Clean S5 data 2
    sprintf(ifo[1].ch1name,       "L1:LSC-STRAIN");
    sprintf(ifo[1].ch1filepath,   datadir);
    sprintf(ifo[1].ch1fileprefix, "L-L1_RDS_C03_L2-");
    sprintf(ifo[1].ch1filesuffix, "-384.gwf");
    ifo[1].ch1filesize   = 384;
    ifo[1].ch1fileoffset = 262;  //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[1].ch1doubleprecision = 0;
    ifo[1].add2channels    = 0;  //0, unless you want to read a signal from file
    
    ifo[1].noiseGPSstart   = 846226064;
    sprintf(ifo[1].noisechannel,    "L1:LSC-STRAIN");
    sprintf(ifo[1].noisefilepath,   datadir);
    sprintf(ifo[1].noisefileprefix, "L-L1_RDS_C03_L2-");
    sprintf(ifo[1].noisefilesuffix, "-384.gwf");
    ifo[1].noisefilesize   = 384;
    ifo[1].noisefileoffset = 262;  //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[1].noisedoubleprecision = 0;
  }
  
  if(run.selectdata == 3) {
    // Playground trigger 845348295
    sprintf(ifo[1].ch1name,       "L1:LSC-STRAIN");
    sprintf(ifo[1].ch1filepath,   datadir);
    sprintf(ifo[1].ch1fileprefix, "L-L1_RDS_C03_L2-");
    sprintf(ifo[1].ch1filesuffix, "-128.gwf");
    ifo[1].ch1filesize   = 128;
    ifo[1].ch1fileoffset = 92;  //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[1].ch1doubleprecision = 0;
    ifo[1].add2channels    = 0;  //0, unless you want to read a signal from file
    
    ifo[1].noiseGPSstart   = 845348573;
    sprintf(ifo[1].noisechannel,    "L1:LSC-STRAIN");
    sprintf(ifo[1].noisefilepath,   datadir);
    sprintf(ifo[1].noisefileprefix, "L-L1_RDS_C03_L2-");
    sprintf(ifo[1].noisefilesuffix, "-128.gwf");
    ifo[1].noisefilesize   = 128;
    ifo[1].noisefileoffset = 92;  //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[1].noisedoubleprecision = 0;
  }
  
  if(run.selectdata == 4) {
    // glitchy data
    sprintf(ifo[1].ch1name,       "L1:LSC-STRAIN");
    sprintf(ifo[1].ch1filepath,   datadir);
    sprintf(ifo[1].ch1fileprefix, "L-L1_RDS_C03_L2-");
    sprintf(ifo[1].ch1filesuffix, "-128.gwf");
    ifo[1].ch1filesize   = 128;
    ifo[1].ch1fileoffset = 125;
    ifo[1].ch1doubleprecision = 0;
    ifo[1].add2channels    = 0;  //0, unless you want to read a signal from file
    
    ifo[1].noiseGPSstart   = 846471038;
    sprintf(ifo[1].noisechannel,    "L1:LSC-STRAIN");
    sprintf(ifo[1].noisefilepath,   datadir);
    sprintf(ifo[1].noisefileprefix, "L-L1_RDS_C03_L2-");
    sprintf(ifo[1].noisefilesuffix, "-128.gwf");
    ifo[1].noisefilesize   = 128;
    ifo[1].noisefileoffset = 125;  //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[1].noisedoubleprecision = 0;
  }
  
  if(run.selectdata == 5) {
    // NINJA data set
    
    sprintf(ifo[1].ch1name,       "L1:STRAIN"); 
    sprintf(ifo[1].ch1filepath,   datadir);
    sprintf(ifo[1].ch1fileprefix, "L1-NINJADATA-");
    sprintf(ifo[1].ch1filesuffix, "-1024.gwf");
    ifo[1].ch1filesize   = 1024; 
    ifo[1].ch1fileoffset = 743;  //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[1].ch1doubleprecision = 0;
    ifo[1].add2channels    = 0;  //0, unless you want to read a signal from file
    
    ifo[1].noiseGPSstart   = 894457315; //894378400;  //894457315;
    sprintf(ifo[1].noisechannel,    "L1:STRAIN");
    sprintf(ifo[1].noisefilepath,   datadir);
    sprintf(ifo[1].noisefileprefix, "L1-NINJADATA-");
    sprintf(ifo[1].noisefilesuffix, "-1024.gwf");
    ifo[1].noisefilesize   = 1024; 
    ifo[1].noisefileoffset = 743;   //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[1].noisedoubleprecision = 0;
  }
  
  
  // PISA, Virgo.  You might want to set the parameter tukeywin from 0.05 to 0.15 when including Virgo data
  sprintf(ifo[2].name, "Pisa");
  ifo[2].lati     = (  43.63/180.0)*pi;
  ifo[2].longi    = (  10.50/180.0)*pi;
  ifo[2].rightarm = ( 341.50/180.0)*pi;
  ifo[2].leftarm  = (  71.50/180.0)*pi;
  ifo[2].radius_eqt  = 6378137.0;    /*6378000.0;*/     /* WGS 84 */
  ifo[2].radius_pole = 6356752.314;  /*6378000.0;*/
  ifo[2].lowCut  =  ifo[0].lowCut;
  ifo[2].highCut = ifo[0].highCut;
  ifo[2].before_tc = ifo[0].before_tc;
  ifo[2].after_tc = ifo[0].after_tc;
  
  if(run.selectdata == 1) {
    // Gaussian, stationary noise
    sprintf(ifo[2].ch1name,       "V1:noise");
    sprintf(ifo[2].ch1filepath,   datadir);
    sprintf(ifo[2].ch1fileprefix, "V-");
    sprintf(ifo[2].ch1filesuffix, "-6000.gwf");
    ifo[2].ch1filesize   = 6000; 
    ifo[2].ch1fileoffset = 4000;   //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[2].ch1doubleprecision = 0;
    ifo[2].add2channels    = 0;   //0, unless you want to read a signal from file
    
    sprintf(ifo[2].ch2name,       "V1:STRAIN_INSP_INJ_ONLY"); 
    sprintf(ifo[2].ch2filepath,   datadir);
    sprintf(ifo[2].ch2fileprefix, "HL-SIM-");
    sprintf(ifo[2].ch2filesuffix, "-6000.gwf");
    ifo[2].ch2filesize   = 6000; 
    ifo[2].ch2fileoffset = 4000; 
    ifo[2].ch2doubleprecision = 0;
    ifo[2].noiseGPSstart   = 700008000;
    
    sprintf(ifo[2].noisechannel,    "V1:noise");
    sprintf(ifo[2].noisefilepath,   datadir);
    sprintf(ifo[2].noisefileprefix, "V-");
    sprintf(ifo[2].noisefilesuffix, "-6000.gwf");
    ifo[2].noisefilesize   = 6000; 
    ifo[2].noisefileoffset = 4000;   //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[2].noisedoubleprecision = 0;
    
    
    //Use LIGO noise for Virgo (keep tukeywin = 0.05)
    if(1==1) {
      printf("   Using LIGO noise for Virgo\n");
      sprintf(ifo[2].ch1name,       "L1:STRAIN"); 
      sprintf(ifo[2].ch1filepath,   datadir);
      sprintf(ifo[2].ch1fileprefix, "HL-SIM-");
      sprintf(ifo[2].ch1filesuffix, "-6000.gwf");
      ifo[2].ch1filesize   = 6000; 
      ifo[2].ch1fileoffset = 4000; 
      ifo[2].ch1doubleprecision = 0;
      ifo[2].add2channels    = 1; 
      
      sprintf(ifo[2].ch2name,       "L1:STRAIN_INSP_INJ_ONLY"); 
      sprintf(ifo[2].ch2filepath,   datadir);
      sprintf(ifo[2].ch2fileprefix, "HL-SIM-");
      sprintf(ifo[2].ch2filesuffix, "-6000.gwf");
      ifo[2].ch2filesize   = 6000; 
      ifo[2].ch2fileoffset = 4000; 
      ifo[2].ch2doubleprecision = 0;
      ifo[2].noiseGPSstart   = 700008000;
      
      sprintf(ifo[2].noisechannel,    "L1:STRAIN");
      sprintf(ifo[2].noisefilepath,   datadir);
      sprintf(ifo[2].noisefileprefix, "HL-SIM-");
      sprintf(ifo[2].noisefilesuffix, "-6000.gwf");
      ifo[2].noisefilesize   = 6000; 
      ifo[2].noisefileoffset = 4000; 
      ifo[2].noisedoubleprecision = 0;
    }
  }

  if(run.selectdata == 5) {
    // NINJA data set
    
    sprintf(ifo[2].ch1name,       "V1:STRAIN");
    sprintf(ifo[2].ch1filepath,   datadir);
    sprintf(ifo[2].ch1fileprefix, "V1-NINJADATA-");
    sprintf(ifo[2].ch1filesuffix, "-1024.gwf");
    ifo[2].ch1filesize   = 1024;
    ifo[2].ch1fileoffset = 743;  //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[2].ch1doubleprecision = 0;
    ifo[2].add2channels    = 0;  //0, unless you want to read a signal from file
    
    ifo[2].noiseGPSstart   = 894457315; //894378400;  //894457315;
    sprintf(ifo[2].noisechannel,    "V1:STRAIN");
    sprintf(ifo[2].noisefilepath,   datadir);
    sprintf(ifo[2].noisefileprefix, "V1-NINJADATA-");
    sprintf(ifo[2].noisefilesuffix, "-1024.gwf");
    ifo[2].noisefilesize   = 1024;
    ifo[2].noisefileoffset = 743;   //If the Frame filename ends in: -839366009-128.gwf, fileoffset = mod(839366009,128)
    ifo[2].noisedoubleprecision = 0;
  }
}



void ifoinit(struct interferometer **ifo, int networksize)
/* Determines interferometer arm (unit-) vectors (and others)                   */
/* given position (lat/long) and arm angles.                                    */
/* Vectors refer to the (right-handed) earth                                    */ 
/* coordinate system spanned by the three vectors:                              */
/*   x) from geocenter to intersection of greenwich meridian with equator plane */
/*   y) from geocenter to intersection of 90E meridian with equator plane       */
/*   z) from geocenter to north pole                                            */
{
  double merinormal[3];  /* normal vector of meridian plane         */
  char latchar[2];
  char longchar[2];
  int i,j;
  double f,logf;
  double flattening, eccentricitySQ, curvatureradius;
  for(i=0; i<networksize; ++i){
    ifo[i]->index = i;
    /* some text output... */
    if(intscrout==1) printf(" | Interferometer %d: `%s'", ifo[i]->index+1, ifo[i]->name);
    if((ifo[i]->lati)  < 0.0) sprintf(latchar, "S");
    else sprintf(latchar, "N");
    if((ifo[i]->longi) < 0.0) sprintf(longchar, "W");
    else sprintf(longchar, "E");
    if(intscrout==1) printf(" at  %1.0f*%2.1f'%s  %1.0f*%2.1f'%s  (%3.0f/%3.0f)\n",
           floor(fabs(ifo[i]->lati*r2d)),  (fabs(ifo[i]->lati*r2d)-floor(fabs(ifo[i]->lati*r2d)))*60.0, latchar,
           floor(fabs(ifo[i]->longi*r2d)), (fabs(ifo[i]->longi*r2d)-floor(fabs(ifo[i]->longi*r2d)))*60.0, longchar,
           360.0 - ifo[i]->rightarm*r2d, 360.0 - ifo[i]->leftarm*r2d);
    if(ifo[i]->ch1doubleprecision) if(intscrout==1) printf(" | frame file precision: double (64 bit)\n"); 
    else if(intscrout==1) printf(" | frame file precision: float (32 bit)\n"); 
    //if(intscrout==1) printf(" | frequency range: %.0f to %.0f Hz.\n", ifo[i]->lowCut, ifo[i]->highCut);
    if(intscrout==1) printf(" | initialising vectors etc...");
    
    
    //Longitude: East is positive
    
    // Place arms on equator plane, so that its designated N-S-direction is aligned with its meridian plane:
    ifo[i]->rightvec[0]  = -cos(ifo[i]->longi + ifo[i]->rightarm);
    ifo[i]->rightvec[1]  = -sin(ifo[i]->longi + ifo[i]->rightarm);
    ifo[i]->rightvec[2]  = 0.0;
    
    ifo[i]->leftvec[0]   = -cos(ifo[i]->longi + ifo[i]->leftarm);
    ifo[i]->leftvec[1]   = -sin(ifo[i]->longi + ifo[i]->leftarm);
    ifo[i]->leftvec[2]   = 0.0;
    
    ifo[i]->normalvec[0] = 0.0;
    ifo[i]->normalvec[1] = 0.0;
    ifo[i]->normalvec[2] = 1.0;
    
    // The following vector is rightarm + 90deg and usually (but not necessarily) identical to the left arm (leftvec):
    ifo[i]->orthoarm[0]  = -cos(ifo[i]->longi + ifo[i]->rightarm + 0.5*pi);
    ifo[i]->orthoarm[1]  = -sin(ifo[i]->longi + ifo[i]->rightarm + 0.5*pi);
    ifo[i]->orthoarm[2]  = 0.0;
    
    // Determine normal vector of meridian plane (i.e., the vector that points E(?) when standing at the equator at longi):
    merinormal[0] = cos(ifo[i]->longi - 0.5*pi);
    merinormal[1] = sin(ifo[i]->longi - 0.5*pi);
    merinormal[2] = 0.0;
    
    // The three vectors:                                                     
    //   x) from geocenter to intersection of ifo meridian with equator plane
    //   y) from geocenter to north pole
    //   z) the above normal vector (merinormal(?))
    // again form another (orthonormal) right-handed system.
    // Now turn all arms clockwise around the normal vector of meridian plane, to account for the latitude of the detectors:
    rotate(ifo[i]->rightvec,  pi/2.0 - ifo[i]->lati, merinormal);
    rotate(ifo[i]->leftvec,   pi/2.0 - ifo[i]->lati, merinormal);
    rotate(ifo[i]->orthoarm,  pi/2.0 - ifo[i]->lati, merinormal);
    rotate(ifo[i]->normalvec, pi/2.0 - ifo[i]->lati, merinormal);
    
    // Initialise the ifo position (!NOT! unit-) vector:
    coord2vec(sin(ifo[i]->lati), ifo[i]->longi, ifo[i]->positionvec);
    if(ifo[i]->radius_eqt < ifo[i]->radius_pole) printf("  CHECK EARTH MODEL RADII !!  ");
    flattening      = (ifo[i]->radius_eqt - ifo[i]->radius_pole) / ifo[i]->radius_eqt;
    eccentricitySQ  = flattening*(2.0-flattening);  /* (squared eccentricity) */
    curvatureradius = ifo[i]->radius_eqt / sqrt(1.0-eccentricitySQ*pow(sin(ifo[i]->lati),2.0));
    ifo[i]->positionvec[0] *= curvatureradius;
    ifo[i]->positionvec[1] *= curvatureradius;
    ifo[i]->positionvec[2] *= curvatureradius*(1.0-eccentricitySQ);
    if(intscrout==1) printf(" ok.\n");
    // printf("== normalvec (%1.2f, %1.2f, %1.2f)\n", ifo[i]->normalvec[0],ifo[i]->normalvec[1],ifo[i]->normalvec[2]);
    // if(intscrout==1) printf(" : f=%f  e^2=%f  v=%f \n", flattening, eccentricitySQ, curvatureradius);
    
    
    printf("   Reading noise...\n");
    noisePSDestimate(ifo[i]);
    
    printf("   Reading data...\n");
    dataFT(ifo,i,networksize);
    
    /* initialise array of different powers of Fourier frequencies */
    /* corresponding to the elements of `ifo[i]->dataTrafo':       */
    /* First loop to determine index bounds & range:               */
    ifo[i]->lowIndex = 0; ifo[i]->highIndex=0;
    for(j=1; j<ifo[i]->FTsize; ++j){
      f = (((double)j)/((double)ifo[i]->FTsize*2.0)) * ((double) ifo[i]->samplerate);
      if((ifo[i]->lowIndex==0)  && (f>=ifo[i]->lowCut)) ifo[i]->lowIndex = j;
      if((ifo[i]->highIndex==0) && (f>ifo[i]->highCut)) ifo[i]->highIndex = j-1;
      /* ...so `lowIndex' and `highIndex' are the extreme indexes WITHIN frequency band */
    }
    ifo[i]->indexRange = ifo[i]->highIndex - (ifo[i]->lowIndex - 1);
    /* `lowIndex' and `highIndex' are the indices analogous to `lowCut' and `highCut' */
    /* but in order to access the respective elements of `raw_dataTrafo'.             */

    ifo[i]->freqpowers = ((double**) malloc(sizeof(double*) * ifo[i]->indexRange));
    ifo[i]->noisePSD   = ((double*) malloc(sizeof(double) * ifo[i]->indexRange));
    ifo[i]->dataTrafo  = ((fftw_complex*) malloc(sizeof(fftw_complex) * ifo[i]->indexRange));
    for(j=0; j<ifo[i]->indexRange; ++j){
      ifo[i]->freqpowers[j] = (double*) malloc(8*sizeof(double));
      f = (((double)(j+ifo[i]->lowIndex))/((double)ifo[i]->FTsize*2.0)) * ((double) ifo[i]->samplerate);
      logf = log(f);
      ifo[i]->freqpowers[j][0] = f;                     /* the `plain' frequency */
      ifo[i]->freqpowers[j][1] = exp((-7.0/6.0)*logf);  /* used in cosine chirp  */
      ifo[i]->freqpowers[j][2] = exp((-5.0/3.0)*logf);  /* used in `psi'         */
      ifo[i]->freqpowers[j][3] = exp(  (-1.0)  *logf);  /* used in `psi'         */
      ifo[i]->freqpowers[j][4] = exp((-2.0/3.0)*logf);  /* used in `psi'         */
      ifo[i]->freqpowers[j][5] = exp((-1.0/3.0)*logf);  /* used in `psi'         */
      ifo[i]->freqpowers[j][6] = logf;                  /* used in `psi'         */
      ifo[i]->freqpowers[j][7] = 2.0*pi*f;              
      ifo[i]->noisePSD[j]      = interpol_log_noisePSD(f,ifo[i]);
      //although smoothing was done for log noise, we store real noise on output
      ifo[i]->noisePSD[j] = exp(ifo[i]->noisePSD[j]);

      ifo[i]->dataTrafo[j]     = ifo[i]->raw_dataTrafo[j+ifo[i]->lowIndex];
    }
    if(intscrout==1) printf(" | %d Fourier frequencies within operational range %.0f--%.0f Hz.\n", 
           ifo[i]->indexRange, ifo[i]->lowCut, ifo[i]->highCut);
    if(i<networksize-1)
      if(intscrout==1) printf(" | --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --\n");
  }
}


void ifodispose(struct interferometer *ifo)
{
  int i;
  if(intscrout==1) printf(" | Interferometer %d `%s' is taken offline.\n", ifo->index, ifo->name);
  free(ifo->raw_noisePSD);       ifo->raw_noisePSD = NULL;
  fftw_free(ifo->raw_dataTrafo); ifo->raw_dataTrafo = NULL;
  free(ifo->noisePSD);           ifo->noisePSD = NULL;
  free(ifo->dataTrafo);          ifo->dataTrafo = NULL;
  for(i=0; i<ifo->indexRange; ++i) free(ifo->freqpowers[i]);
  free(ifo->freqpowers);         ifo->freqpowers = NULL;  
  fftw_destroy_plan(ifo->FTplan);
  fftw_free(ifo->FTin);          ifo->FTin = NULL;
  fftw_free(ifo->FTout);         ifo->FTout = NULL;
  free(ifo->FTwindow);           ifo->FTwindow = NULL;
}






// *** Routines that do data I/O and data handling ***

double *filter(int *order, int samplerate, double upperlimit)
/* Computes FIR filter coefficients.                                */
/* The filter has a total of N=(order+order-1) coefficients         */
/* that are symmetric, i.e. coef[k] = coef[N-k].                    */
/* For details & filter properties see downsampling function below. */
{
  /*-- specify filter characteristics: --*/
  int ncoef      = 129;  /* number of filter coefficients... 129 should be sufficient   */
  int totalcoef  = ncoef+ncoef-1;   /* total number of coefficients                    */
  double desired[2] = {1.0, 0.0};      /* desired gain                                    */
  double weights[2] = {1.0, 1.0};      /* weight for `loss' in pass- & stopband           */
  double transitionbandwidth=0.0125;  //suggested by Christian Roever via 07/30/08 e-mail
  //place transition bandwidth half-way between upper edge of pass band, which is
  //(upperlimit/samplerate) in relative units, and new Nyquist frequency, which is
  //0.5/downsamplefactor in relative units.
  //band vector contains 0, end of pass band, end of transition band, 0.5 (old Nyquist)
  //if there's not enough room for a transition band, we have a problem!
  if(0.5/downsamplefactor - upperlimit/((double)samplerate) < transitionbandwidth){
	printf(" Problem in filter() function while downsampling!\n");
	printf(" New Nyquist frequency after downsampling by %g is %g\n",
			downsamplefactor, ((double)samplerate)/downsamplefactor/2.0);
	printf(" Desired upper limit is %g\n", upperlimit);
	printf(" This doesn't leave enough room for a transition band of relative width %g\n",
		transitionbandwidth);
	printf(" Aborting!\n");
	exit(1);
  }
  double endpassband= upperlimit/((double)samplerate) +
	(0.5/downsamplefactor - upperlimit/((double)samplerate) - transitionbandwidth)/2.0;
  double endtransitionband=0.5/downsamplefactor - 
	(0.5/downsamplefactor - upperlimit/((double)samplerate) - transitionbandwidth)/2.0;
  double bands[4]   = {0.0, endpassband, endtransitionband, 0.5};
  double *coef;                        /* vector of coefficients (symmetric)               */
  coef = (double*) malloc(sizeof(double)*totalcoef);
  /*-- determine filter coefficients: --*/
  remez(coef, totalcoef, 2, bands, desired, weights, BANDPASS);
  *order = ncoef;
  return coef;
}


double *downsample(double data[], int *datalength, double filtercoef[], int ncoef)
/* Downsamples a time series by factor downsamplefactor by first low-pass   */
/* filtering it using a finite-impulse-response (FIR) filter */
/* and then thinning the data.                               */
/* Filter coefficients are determined using the              */
/* `Parks-McClellan' or `Remez exchange' algorithm.          */
/* Resulting data vector is shorter than original            */
/* Returned vector is allocated using `fftw_malloc()' and    */
/* thus must be freed again using `fftw_free()'.             */
{
     int flength = *datalength-2*(ncoef-1);
     int tlength = (int)ceil(((double)flength)/downsamplefactor);
  double *thinned;      /* vector of filtered & thinned data */
     int i,j,k;

  /*-- filter & thin data: --*/
  thinned = (double*) fftw_malloc(sizeof(double) * tlength);
  k = 0;
  for(i=ncoef-1; i<*datalength-ncoef; i+=(int)downsamplefactor) {
    thinned[k] = filtercoef[ncoef-1]*data[i];
    for(j=1; j<ncoef; ++j)
      thinned[k] += filtercoef[(ncoef-1)+j]*(data[i-j]+data[i+j]);
    ++k;
  }

  /*-- return results --*/
  *datalength = tlength;
  return thinned;
}



double hann(int j, int N)
/* `hann window' for windowing data  */
/* j = 0, ..., N-1                   */
{
  return 0.5*(1.0-cos(((double)j/(double)N)*2.0*pi));
}

double tukey(int j, int N, double r)
/* Tukey window... for r=0 equal to rectangular window, for r=1 equal to Hann window. */
/* ( 0 < r < 1 denotes the fraction of the window in which it behaves sinusoidal)     */
/* j = 0, ..., N-1                                                                    */
{
  double win = 1.0;
  if(((double)j) > (((double)N)/2.0)) j = N-j;
  if(((double)j) < (r*(((double)N)/2.0)))
    win = 0.5*(1.0-cos(((2.0*pi)/r)*(((double)j)/((double)N))));
  return win;
}




void dataFT(struct interferometer *ifo[], int i, int networksize)
// Computes the Fourier Transform for the specified range of the specified Frame (".gwf") file,
// after adding up the two (signal & noise) channels, or injecting a waveform template into the noise.
// Also takes care of preparing FT stuff  (ifo[i]->FTplan, ->FTin, ->FTout, ...).
{
  if(MvdSdebug) printf("  DataFT\n");
  struct FrFile *iFile=NULL;                // Frame file(s)
  struct FrVect *svect=NULL, *nvect=NULL;   // data vectors (signal & noise)
  int           N;                          // size of input
  double        *raw;                       // downsampling input
  int           j, ncoef;
  double        *filtercoef;
  long          filestart;
  char          filenames[1000]="";
  int           filecount = 0;
  double        from, to, delta;
  double        *injection;
  
  // `from' and `to' are determined so that the range specified by `before_tc' and `after_tc'
  // falls into the flat part of the (Tukey-) window:                                        
  from  = floor(prior_tc_mean - ifo[i]->before_tc - (ifo[i]->before_tc+ifo[i]->after_tc) * 0.5 * (tukeywin/(1.0-tukeywin)));
  to    =  ceil(prior_tc_mean + ifo[i]->after_tc  + (ifo[i]->before_tc+ifo[i]->after_tc) * 0.5 * (tukeywin/(1.0-tukeywin)));
  delta = (to) - (from);
  if(intscrout==1) printf(" | investigated time range : from %.1f to %.1f (%.1f seconds)\n", from, to, delta);
  
  // Starting time of first(!) Frame file to be read:
  filestart = (((((long)(from))-ifo[i]->ch1fileoffset) / ifo[i]->ch1filesize) * ifo[i]->ch1filesize) + ifo[i]->ch1fileoffset;
  //if(intscrout==1) printf(" | chirp file(s) to be read:\n");
  
  // Assemble the filename character string:
  while (((double)filestart) < to){
    if(filecount == 0)
      sprintf(filenames, "%s/%s%ld%s", ifo[i]->ch1filepath, ifo[i]->ch1fileprefix, (long)filestart, ifo[i]->ch1filesuffix);  // Fill in filename etc. for first file
    else
      sprintf(filenames, "%s %s/%s%ld%s", filenames, ifo[i]->ch1filepath, ifo[i]->ch1fileprefix, (long)filestart, ifo[i]->ch1filesuffix);  // Append filename etc. for following files
    //if(intscrout==1) printf(" |   %s%ld%s\n", ifo[i]->framefileprefix, (long)filestart, ifo[i]->framefilesuffix);
    filestart += ifo[i]->ch1filesize;
    filecount += 1;
  }
  
  // Open frame file(s):
  if(intscrout==1) printf(" | opening %d chirp data file(s)... \n", filecount);
  iFile = FrFileINew(filenames);
  if(iFile == NULL) {
    printf("\n\n   ERROR opening data file: %s (channel 1), aborting.\n\n\n",filenames);
    exit(1);
  }
  if(intscrout==1) printf(" | %s\n",filenames);
  // Read 1st channel (noise or noise+signal):
  if(ifo[i]->ch1doubleprecision) //Seems precision is single
    nvect = FrFileIGetVectD(iFile, ifo[i]->ch1name, from, delta);
  else
    nvect = FrFileIGetVectF(iFile, ifo[i]->ch1name, from, delta);
  if(nvect == NULL) {
    printf("\n\n   ERROR reading data file: %s (channel 1), aborting.\n\n\n",filenames);
    exit(1);
  }

  FrFileIEnd(iFile);
  
  N = nvect->nData;
  ifo[i]->samplerate = (int)(1.0 / (nvect->dx[0]) + 0.5);  // Add 0.5 for correct truncation/rounding
  if(intscrout==1) printf(" | original sampling rate: %d Hz\n", ifo[i]->samplerate);
  
  // Inject the signal into the noise
  if(inject) {
    if(intscrout==1) printf(" :  injecting signal:\n");
    
    // Define injection parameters:
    struct parset injectpar;
    gettrueparameters(&injectpar);
    double root = sqrt(0.25-injectpar.eta);
    double fraction = (0.5-root) / (0.5+root);
    double inversefraction = 1.0/fraction;
    double m1 = injectpar.mc * (pow(1+fraction,0.2) / pow(fraction,0.6));
    double m2 = injectpar.mc * (pow(1+inversefraction,0.2) / pow(inversefraction,0.6));
    injectpar.loctc           = NULL;
    injectpar.localti         = NULL;
    injectpar.locazi          = NULL;
    injectpar.locpolar        = NULL;
    injectpar.loctc    = (double*)calloc(networksize,sizeof(double));
    injectpar.localti  = (double*)calloc(networksize,sizeof(double));
    injectpar.locazi   = (double*)calloc(networksize,sizeof(double));
    injectpar.locpolar = (double*)calloc(networksize,sizeof(double));
    
    if(intscrout==1) {
      printf(" :   m1 = %.1f Mo,  m2 = %.1f Mo  (Mc = %.3f Mo,  eta = %.4f)\n", m1, m2, injectpar.mc, injectpar.eta);
      printf(" :   tc = %.4f s,  dist = %.1f Mpc\n", injectpar.tc, exp(injectpar.logdl));
      printf(" :   ra = %.2f h,  dec = %.2f deg  (GMST = %.2f h)\n",(rightAscension(injectpar.longi,GMST(injectpar.tc))/pi)*12.0, (asin(injectpar.sinlati)/pi)*180.0, (GMST(injectpar.tc)/pi)*12.0);
      printf(" :   phase = %.2f rad\n", injectpar.phase);
    }
    ifo[i]->FTstart = from; // Temporary setting so `parupdate()' works properly
    
    localpar(&injectpar, ifo, networksize);
    
    if(intscrout==1) {
      printf(" :   local parameters:\n");
      printf(" :   tc           = %.5f s\n",injectpar.loctc[i]+from);
      printf(" :   altitude     = %.2f rad\n",injectpar.localti[i]);
      printf(" :   azimuth      = %.2f rad\n",injectpar.locazi[i]);
      printf(" :   polarisation = %.2f rad\n",injectpar.locpolar[i]);
    }
    
    
    // Initialise further components:
    
    // Generate template:
    injection = malloc(sizeof(double) * N);
    double *tempinj = ifo[i]->FTin;
    double tempfrom = ifo[i]->FTstart;
    int tempN = ifo[i]->samplesize;
    ifo[i]->FTin = injection;
    ifo[i]->FTstart = from;
    ifo[i]->samplesize = N;
    template(&injectpar,ifo,i);
    ifo[i]->FTin = tempinj;
    ifo[i]->FTstart = tempfrom;
    ifo[i]->samplesize = tempN;
    
    // Write high-sampled injection signal to disc
    if(writesignal && 1==2){
      char filename[100]="";
      sprintf(filename,"%s-injection.dat",ifo[i]->name);
      FILE *dump = fopen(filename,"w");
      fprintf(dump,"%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n","m1","m2","mc","eta","tc","dl","lat","lon","phase","spin","kappa","thJ0","phJ0","alpha");
      fprintf(dump,"%12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	      injectpar.m1,injectpar.m2,injectpar.mc,injectpar.eta,injectpar.tc,exp(injectpar.logdl),asin(injectpar.sinlati)*r2d,injectpar.longi*r2d,injectpar.phase,
              injectpar.spin,injectpar.kappa,injectpar.sinthJ0,injectpar.phiJ0,injectpar.alpha);
      fprintf(dump,"t Ht\n");
      for(j=0; j<N; ++j)
	fprintf(dump, "%9.9f %.6e\n", from+(((double)j)/((double) (ifo[i]->samplerate))), injection[j]);  //*******************************************************************************************
      fclose(dump);
      if(intscrout==1) printf(" : (signal written to file)\n");
    }
    
    freeparset(&injectpar);
  } // if(inject)
  else if(ifo[i]->add2channels) { // Read 2nd channel (signal only)  (if not doing a software injection)
    filestart = (((((long)(from))-ifo[i]->ch2fileoffset) / ifo[i]->ch2filesize) * ifo[i]->ch2filesize) + ifo[i]->ch2fileoffset;
    
    // Assemble the filename character string:
    sprintf(filenames, " "); filecount = 0;
    while (((double)filestart) < to){
      if(filecount == 0)
        sprintf(filenames, "%s/%s%ld%s", ifo[i]->ch2filepath, ifo[i]->ch2fileprefix, (long)filestart, ifo[i]->ch2filesuffix);  // Fill in filename etc. for first file
      else
        sprintf(filenames, "%s %s/%s%ld%s", filenames, ifo[i]->ch2filepath, ifo[i]->ch2fileprefix, (long)filestart, ifo[i]->ch2filesuffix);  // Append filename etc. for following files:
      filestart += ifo[i]->ch2filesize;
      filecount += 1;
    }
    
    // Open file:
    iFile = FrFileINew(filenames);
    if(iFile == NULL) {
      printf("\n\n   ERROR opening data file: %s (channel 2), aborting.\n\n\n",filenames);
      exit(1);
    }
    
    if(ifo[i]->ch2doubleprecision) {
      svect = FrFileIGetVectD(iFile, ifo[i]->ch2name, from, delta);
    } else {
      svect = FrFileIGetVectF(iFile, ifo[i]->ch2name, from, delta);
    }
    if(svect == NULL) {
      printf("\n\n   ERROR reading data file: %s (channel 2), aborting.\n\n\n",filenames);
      exit(1);
    }
    FrFileIEnd(iFile);
  } //End if not doing a software injection
  
  
  // Allocate memory for transform input:
  raw  = malloc(sizeof(double) * N);
  // Fill in values:
  for(j=0; j<N; ++j) raw[j] = nvect->dataF[j];
  
  // Add channels (noise plus signal):
  if(inject){
    for(j=0; j<N; ++j) raw[j] += injection[j];
  } else if(ifo[i]->add2channels) {
    for(j=0; j<N; ++j) raw[j] += svect->dataF[j];
  }
  
  
  // Release the FrVect objects:
  FrVectFree(svect);
  FrVectFree(nvect);
  if(inject) free(injection);
  
  
  int screwcount = 0;
  for(j=0; j<N; ++j)
    if(!(raw[j]<HUGE_VAL)) ++screwcount;
  if(screwcount>0){
    printf(" : %d missing data points in DATA file(s) !!\n",screwcount);
    printf(" : (maybe the precision is incorrect)\n");
  }
 
  // Downsample (by factor downsamplefactor):    *** changes value of N ***
  if(intscrout==1) printf(" | downsampling... \n");

/*NINJA*/
  filtercoef = filter(&ncoef, ifo[i]->samplerate, ifo[i]->highCut);
  ifo[i]->FTin = downsample(raw, &N, filtercoef, ncoef);

  ifo[i]->FTstart = from + ((double)(ncoef-1))/((double)(ifo[i]->samplerate));
  ifo[i]->deltaFT = delta - ((double)((ncoef-1)*2))/((double)(ifo[i]->samplerate));

  ifo[i]->samplesize = N;
  ifo[i]->FTsize = (N/2)+1;
  ifo[i]->samplerate = (int)((double)ifo[i]->samplerate/downsamplefactor);
  free(raw);
  free(filtercoef);

/*
  ifo[i]->FTin = (double*) fftw_malloc(sizeof(double)*N);
  for (j=0; j<N; ++j)
    ifo[i]->FTin[j] = raw[j];
  ifo[i]->FTstart = from;
  ifo[i]->deltaFT = delta;
  ifo[i]->samplesize = N;
  ifo[i]->FTsize = (N/2)+1;  
  free(raw);
*/ 

  // Window input data with a Tukey window:
  ifo[i]->FTwindow = malloc(sizeof(double) * N);
  for(j=0; j<N; ++j){
    ifo[i]->FTwindow[j] =  tukey(j, N, tukeywin);
    ifo[i]->FTin[j] *= ifo[i]->FTwindow[j];
  }
  
  
  // Write windowed, time-domain data (signal + noise) to disc
  if(writesignal && 1==1){
    char filename[1000]="";
    sprintf(filename, "%s-data.dat", ifo[i]->name);  //Write in current dir
    FILE *dump = fopen(filename,"w");
    
    // Get true signal parameters and write them to the header
    struct parset par;
    gettrueparameters(&par);
    fprintf(dump,"%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n","m1","m2","mc","eta","tc","dl","lat","lon","phase","spin","kappa","thJ0","phJ0","alpha");
    fprintf(dump,"%12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	    par.m1,par.m2,par.mc,par.eta,par.tc,exp(par.logdl),asin(par.sinlati)*r2d,par.longi*r2d,par.phase,par.spin,par.kappa,par.sinthJ0,par.phiJ0,par.alpha);
    fprintf(dump,"       GPS time (s)         H(t)\n");
    for(j=0; j<N; ++j)
      fprintf(dump, "%9.9f %.6e\n", ifo[i]->FTstart+(((double)j)/((double) (ifo[i]->samplerate))), ifo[i]->FTin[j]);
    fclose(dump);
    if(intscrout) printf(" : (signal written to file)\n");
    freeparset(&par);
  }
  
  
  // Allocate memory for Fourier-transform output:
  ifo[i]->FTout = fftw_malloc(sizeof(fftw_complex) * (ifo[i]->FTsize));  
  
  // Contruct a FFTW plan:
  ifo[i]->FTplan = fftw_plan_dft_r2c_1d(N, ifo[i]->FTin, ifo[i]->FTout, FFTW_ESTIMATE);
  //ifo[i]->FTplan = fftw_plan_dft_r2c_1d(N, ifo[i]->FTin, ifo[i]->FTout, FFTW_MEASURE);  //This must be done before initialisation of FTin and could optimise the FFT
  
  // Compute the FFT:
  if(intscrout==1) printf(" | performing data Fourier transform (%.1f s at %d Hz)... ",delta, ifo[i]->samplerate);
  fftw_execute(ifo[i]->FTplan);
  if(ifo[i]->FTout == NULL){
    printf("\n\n   ERROR performing Fourier transform: %s, aborting.\n\n\n",filenames);
    exit(1);
  }
  else if(intscrout==1) printf("ok.\n");
  
  
  // Normalise transform (divide by sampling rate, see Mark's code, line 184):
  for(j=0; j<ifo[i]->FTsize; ++j) ifo[i]->FTout[j] /= (double)ifo[i]->samplerate;
  
  // Copy to 'raw_dataTrafo':
  ifo[i]->raw_dataTrafo = fftw_malloc(sizeof(fftw_complex) * (ifo[i]->FTsize));  
  for(j=0; j<ifo[i]->FTsize; ++j) ifo[i]->raw_dataTrafo[j] = ifo[i]->FTout[j];
  
  
  
  
  // Write data FFT to disc (i.e., FFT or amplitude spectrum of signal+noise)
  if(writesignal && 1==1){
    double f=0.0;
    double complex tempvar=0.0;
    char filename[1000]="";
    sprintf(filename, "%s-dataFFT.dat", ifo[i]->name);  //Write in current dir
    FILE *dump1 = fopen(filename,"w");
    
    // Get true signal parameters and write them to the header
    struct parset par;
    gettrueparameters(&par);
    fprintf(dump1,"%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n","m1","m2","mc","eta","tc","dl","lat","lon","phase","spin","kappa","thJ0","phJ0","alpha");
    fprintf(dump1,"%12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	    par.m1,par.m2,par.mc,par.eta,par.tc,exp(par.logdl),asin(par.sinlati)*r2d,par.longi*r2d,par.phase,par.spin,par.kappa,par.sinthJ0,par.phiJ0,par.alpha);
    fprintf(dump1,"       f (Hz)    real(H(f))    imag(H(f))\n");
    
    double fact1a = ((double)ifo[i]->samplerate) / (2.0*(double)ifo[i]->FTsize);
    double fact1b = sqrt(2.0)*2.0/ifo[i]->deltaFT;  //Extra factor of sqrt(2) to get the numbers right with the outside world
    // Loop over the Fourier frequencies 
    for(j=1; j<ifo[i]->FTsize; ++j){
      f = fact1a * ((double)(j+ifo[i]->lowIndex));
      //if(f>0.9*ifo[i]->lowCut) fprintf(dump1, "%9.9f %.6e\n", log10(f), log10(2.0*cabs( ifo[i]->raw_dataTrafo[j] )/ifo[i]->deltaFT )  );
      tempvar = fact1b * ifo[i]->raw_dataTrafo[j];
      if(f>0.9*ifo[i]->lowCut) fprintf(dump1, "%13.6e %13.6e %13.6e\n", f, creal(tempvar), cimag(tempvar) );  //Save the real and imaginary parts of the data FFT
    }
    fclose(dump1);
    if(intscrout) printf(" : (data FFT written to file)\n");
  }
  
  
}
// End dataFT()



void noisePSDestimate(struct interferometer *ifo)
/* returns a (smoothed) estimate of the log- Power Spectral Density.  */
/* data is split into K segments of M seconds,                        */
/* and K-1 overlapping segments of length 2M are eventually           */
/* windowed and transformed                                           */
  
{
  struct FrFile  *iFile=NULL;  /* Frame File                           */
  struct FrVect   *vect=NULL; 
  double           *raw=NULL;  /* FFTW vectors etc.                    */
  double            *in=NULL;
  fftw_complex     *out=NULL;
  fftw_plan           FTplan;
  double           *PSD=NULL;  /* vector containing PSD                */
  double          *sPSD=NULL;  /* vector containing smoothed PSD       */
  double           *win=NULL;  /* window                               */
  double      Mseconds=  8.0;  /* M expressed in seconds               */
  double      Nseconds=256.0;  /* total seconds                        */
  int  K=(int)(Nseconds/Mseconds);  /* number of 8-second-segments          */
  double wss=0.0,log2=log(2.0);  /* squared & summed window coefficients  etc.*/
  int     i, j, M, N, dummyN;
  int             samplerate;
  int           lower, upper;  /* indices of lower & upper frequency bounds in FT vector */
  double             nyquist;  /* the critical nyquist frequency       */
  int          smoothrange=0;  /* number of samples to left & right that are averaged (=0!!)*/
  int               PSDrange;
  double                 sum;
  int                 FTsize;
  double         *filtercoef;
  int                  ncoef; 
  char    filenames[2000];
  long             filestart;
  int            filecount=0;
  
  
  /* starting time of first(!) frame file to be read: */
  filestart = (((ifo->noiseGPSstart-ifo->noisefileoffset) / ifo->noisefilesize) * ifo->noisefilesize) + ifo->noisefileoffset;
  /* Assemble the filename character string: */
  while (((double)filestart) < (((double)ifo->noiseGPSstart)+Nseconds)){
    if(filecount == 0) /* fill in filename for first file: */
      sprintf(filenames,"%s/%s%ld%s",ifo->noisefilepath,ifo->noisefileprefix,(long)filestart,ifo->noisefilesuffix);
    else /* append filename for following files: */
      sprintf(filenames,"%s %s/%s%ld%s",filenames,ifo->noisefilepath,ifo->noisefileprefix,(long)filestart,ifo->noisefilesuffix);
    filestart += ifo->noisefilesize;
    filecount += 1;
  }
  
  /*-- open frame file --*/
  if(intscrout==1) printf(" | opening %d noise data file(s)... \n",filecount);
  if(intscrout==1) printf(" | %s\n",filenames);
  iFile = FrFileINew(filenames);
  if(iFile == NULL) {
    printf("\n\n   ERROR opening noise data file: %s, aborting.\n\n\n",filenames);
    exit(1);
  }
  /*else if(intscrout==1) printf("ok.\n");*/
  
  if(intscrout==1) printf(" | estimating noise PSD... ");
  /*-- read first two bits (2M seconds) --*/
  /*-- access (noise) channel           --*/
  if(ifo->noisedoubleprecision)
    vect = FrFileIGetVectD(iFile, ifo->noisechannel, ((double)ifo->noiseGPSstart), Mseconds*2);
  else
    vect = FrFileIGetVectF(iFile, ifo->noisechannel, ((double)ifo->noiseGPSstart), Mseconds*2);
  if(vect == NULL) {
    printf("\n\n   ERROR reading noise data file: %s, aborting.\n\n\n",filenames);
    exit(1);
  }
  N = vect->nData; /* length of filtered & downsampled data (not yet!) */
  M = (int)(N/2.0);
  samplerate = (int)(1.0 / (vect->dx[0]) +0.5);
  /* add 0.5 for correct truncation/rounding */
  
  /*-- copy data to vector `raw' --*/
  raw = (double*)malloc(sizeof(double)*N);
  for(i=0; i<N; ++i)
    raw[i] = vect->dataF[i];
  
  int screwcount = 0;
  for(i=0; i<N; ++i)
    if(!(raw[i]<HUGE_VAL))
      ++screwcount;
  
  /*-- DOWNSAMPLE (by factor downsamplefactor)           --*/
  /*-- !! changes value of `N' as well !! --*/
  /*NINJA*/ 
  filtercoef = filter(&ncoef, samplerate, ifo->highCut);
  in = downsample(raw, &N, filtercoef, ncoef);
  FTsize = (N/2)+1;
  samplerate = (int)((double)samplerate/downsamplefactor);
/*
  in = (double*) fftw_malloc(sizeof(double)*N);
  for (i=0; i<N; ++i)
    in[i] = raw[i];
  FTsize = (N/2)+1;    
*/

  nyquist      = ((double)samplerate)/2.0;
  lower        = (int)(floor((ifo->lowCut/nyquist)*(FTsize-1)));
  upper        = (int)(ceil((ifo->highCut/nyquist)*(FTsize-1)));
  ifo->PSDsize = upper-lower;
  PSDrange     = ifo->PSDsize+2*smoothrange;

  /*-- allocate memory for window: --*/
  win = (double *) malloc(sizeof(double) * N);
  /*-- compute windowing coefficients: --*/
  for(i=0; i<N; ++i) {
    win[i] = hann(i, N);
    wss += win[i]*win[i];
  }
  /*-- normalise window coefs as in Mark's/Nelson's code (see line 695): --*/
  for(i=0; i<N; ++i) {
    win[i] /= sqrt(wss * ((double)(K-1)) * ((double)samplerate));
  }
  /*-- apply to data: --*/
  for(i=0; i<N; ++i) {
    in[i] *= win[i];
  }
  /*-- put last M values in first M places for following iterations: --*/
  /*-- (`vect' vectors in later iterations are only half as long)    --*/
  for(i=0; i<M; ++i)
    vect->dataF[i] = vect->dataF[i+M];
  /*-- allocate memory for transform output --*/
  out = fftw_malloc(sizeof(fftw_complex) * FTsize);  
  /*-- contruct a transform plan --*/ 
  FTplan = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
  /*-- (`FFTW_MEASURE' option not appropriate here.) --*/
  /*-- transform --*/
  fftw_execute(FTplan);
  if(out == NULL){
    printf("\n\n   ERROR performing noise Fourier transform: %s, aborting.\n\n\n",filenames);
    exit(1);
  }
  fftw_free(in); 
  /*-- allocate & initialize PSD vector --*/
  PSD = (double *) malloc(PSDrange*sizeof(double));
  for(i=0; i<PSDrange; ++i)
    PSD[i] = pow(cabs(out[(lower-smoothrange)+i]), 2.0);

  // Read segments 3 to K:
  for(j=3; j<=K; ++j){ 
    /*-- copy first half of data from previous iteration (M seconds) --*/
    for(i=0; i<M; ++i) 
      raw[i] = vect->dataF[i];
    /*-- read 2nd half of data (again, M seconds) --*/
    FrVectFree(vect);
    if(ifo->noisedoubleprecision)
      vect = FrFileIGetVectD(iFile, ifo->noisechannel, ((double)ifo->noiseGPSstart)+((double)(j-1))*Mseconds, Mseconds);
    else
      vect = FrFileIGetVectF(iFile, ifo->noisechannel, ((double)ifo->noiseGPSstart)+((double)(j-1))*Mseconds, Mseconds);
    if(vect == NULL) {
      printf("\n : error accessing noise channel!\n");
    }
    
    // Copy 2nd half of data:
    for(i=0; i<M; ++i) 
      raw[i+M] = vect->dataF[i];
    
    for(i=0; i<(2*M); ++i)
      if(!(raw[i]<HUGE_VAL))
        ++screwcount;
    
    // Downsample:
    /*NINJA*/ 
    dummyN = 2*M;
    in = downsample(raw, &dummyN, filtercoef, ncoef);
/*
    in = (double*) fftw_malloc(sizeof(double)*N);
    for (i=0; i<N; ++i)
      in[i] = raw[i];
*/

    // Window data:
    for(i=0; i<N; ++i)
      in[i] *= win[i];
    
    // Execute FT:
    fftw_destroy_plan(FTplan); // Previous `in'-vector was freed in the meantime
    FTplan = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
    fftw_execute(FTplan);
    if(out == NULL){
      printf("\n : error performing (noise) Fourier transform!\n");
    }
    fftw_free(in); 
    
    // Add to PSD vector
    for(i=0; i<PSDrange; ++i)
      PSD[i] += pow(cabs(out[(lower-smoothrange)+i]),2.0);
    ///*xyz*/
    //for(i=0; i<FTsize; ++i)
    //  outPSD[i] += pow(cabs(out[i]), 2.0);
    ///*xyz*/
  }
  FrVectFree(vect);
  FrFileIEnd(iFile);
  fftw_destroy_plan(FTplan);
  fftw_free(out);
  free(win);
  free(filtercoef);
  
  
  // `PSD' now contains the squared FT magnitudes, summed over (K-1) segments
  
  // Normalise & log the summed PSD's:
  for(i=0; i<PSDrange; ++i)
    PSD[i] = log(PSD[i]) + log2;
  
  if(intscrout==1) printf("ok.\n");
  /*if(intscrout==1) printf(" | averaged over %d overlapping segments of %1.0f s each (%.0f s total).\n", 
    K-1, Mseconds*2, Nseconds);*/
  if(screwcount>0){
    printf(" : %d missing data points in NOISE file(s) !!\n",screwcount);
    printf(" : (maybe the precision is incorrect)\n");
  }
  
  // Smooth PSD:
  sPSD = (double *) malloc((ifo->PSDsize)*sizeof(double));
  for(i=smoothrange; i<(PSDrange-smoothrange); ++i) {
    sum = 0.0;
    for(j=-smoothrange; j<=smoothrange; ++j)
      sum += PSD[i+j];
    sPSD[i-smoothrange] = sum / (2.0*smoothrange+1.0);
  }
  if(smoothrange>0) {
    if(intscrout==1) printf(" | and a range of +/- %0.2f Hz.\n", (((double)smoothrange)/((double)FTsize))*nyquist);
  }
  
  // PSD estimation finished
  free(PSD);
  ifo->raw_noisePSD = sPSD;
  free(raw);
}


double interpol_log_noisePSD(double f, struct interferometer *ifo)
// Returns linearly interpolated (log-) noise PSD.
{
  double dblindex = (((f-ifo->lowCut)/(ifo->highCut-ifo->lowCut)) * (double)ifo->PSDsize);
  int lowindex    = (int)dblindex;   // (truncated!)
  int highindex   = lowindex + 1;
  double weight1  = ((double)highindex) - dblindex;
  double weight2  = dblindex - ((double)lowindex);
  double result   = weight1*ifo->raw_noisePSD[lowindex] + weight2*ifo->raw_noisePSD[highindex];
  return result;
}
