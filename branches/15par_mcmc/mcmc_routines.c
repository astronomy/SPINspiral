#include <mcmc.h>


//***  Different routines for the spinning MCMC code  ***



//*** MASSES ***//

double massratio(double m1, double m2)
/* mass ratio (eta) for given masses */
{
  return (m1*m2)/pow(m1+m2,2.0);
}

void mc2masses(double mc, double eta, double *m1, double *m2)
/* total mass (mt) for given mass ratio & chirp mass */
{
  double root = sqrt(0.25-eta);
  double fraction = (0.5+root) / (0.5-root);
  double inversefraction = (0.5-root) / (0.5+root);
  *m1 = mc * (pow(1+fraction,0.2) / pow(fraction,0.6));
  *m2 = mc * (pow(1+inversefraction,0.2) / pow(inversefraction,0.6));
}

double chirpmass(double m1, double m2)
/* chirp mass (mc) for given masses */
{
  return pow(m1*m2,0.6) / pow(m1+m2,0.2);
}





//*** GMT and RIGHT ASCENSION ***

double GMST(double GPSsec)
// Derives the `Greenwich Mean Sidereal Time' (in radians!) from GPS time (in seconds).
// (see K.R.Lang(1999), p.80sqq.)
{
  double GPS_Jan1st2000midnight    = 630720013.0;
  double leapseconds = 32.0;  /* at Jan 1st 2000 */
  double seconds, days, centuries, secCurrentDay, result;
  if (GPSsec > (GPS_Jan1st2000midnight + 189388800.0)) leapseconds += 1.0; // One more leapsecond after 2005/'06
  if (GPSsec < 630720013.0) printf(" : WARNING: GMST's before 1.1.2000 are inaccurate! \n");
  
  // Time since Jan 1st 2000 (0:00h):
  seconds       = (GPSsec - GPS_Jan1st2000midnight) + (leapseconds - 32.0);
  days          = floor(seconds/86400.0) - 0.5;
  secCurrentDay = fmod(seconds, 86400.0);
  centuries     = days / 36525.0;
  result  = 24110.54841+(centuries*(8640184.812866+centuries*(0.093104+centuries*6.2e-6)));
  result += secCurrentDay * 1.002737909350795; // UTC day is 1.002 * MST day
  result  = fmod(result/86400.0,1.0);
  result *= 2.0*pi;
  return result;
}

double rightAscension(double longi, double GMST)
// Derives right ascension (in radians!) from longitude given GMST (radians).
{
  double result = fmod(longi + GMST + mtpi,tpi);  //Bring it between 0 and 2pi
  return result;
}

double longitude(double rightAscension, double GMST)
// Derives longitude from right ascension (radians), given GMST (radians).  Actually this is something like the Greenwich hour angle of the corresponding RA
{
  double result = fmod(rightAscension - GMST + mtpi,tpi);
  return result;
}






//*** VECTORS ***

double dotproduct(double vec1[3], double vec2[3])
// Dot product of two vectors
{
  return vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
}

void facvec(double vec1[3], double fac, double vec2[3])
// Multiply a vector by a factor (scalar): vec2 = fac*vec1
{
  vec2[0] = fac*vec1[0];
  vec2[1] = fac*vec1[1];
  vec2[2] = fac*vec1[2];
}

void addvec(double vec1[3], double vec2[3], double result[3])
// Add two vectors result = vec1 + vec2
{
  result[0] = vec1[0] + vec2[0];
  result[1] = vec1[1] + vec2[1];
  result[2] = vec1[2] + vec2[2];
}


void normalise(double vec[3])
/*   vec  :=  vec / |vec|   */
{
  double length = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
  vec[0] /= length;
  vec[1] /= length;
  vec[2] /= length;
}

void crossproduct(double vec1[3], double vec2[3], double result[3])
/* cross product (cartesian p., outer p.) of two vectors: */
/*    -  x*y is orthogonal to x and y.                    */
/*    -  |x*y| = |x| * |y| * sin(angle(x,y)).             */
/*    -  x, y, and x*y form a right-handed system.        */
/* note:  x*y  !=  y*x  (not commutative)                 */
{
  result[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
  result[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
  result[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
}

void rotate(double x[3], double angle, double axis[3])
/* rotates vector x clockwise around `axis'               */
/* (looking along axis while it is pointing towards you). */
/*   !!  `axis' must be a UNIT VECTOR  !!                 */
{
  int i, j;
  double cosa = cos(-angle);
  double sina = sin(-angle);
  double omcosa = 1.0 - cosa;
  double R[3][3] = {{cosa+axis[0]*axis[0]*omcosa, 
                     axis[0]*axis[1]*omcosa-axis[2]*sina,
                     axis[0]*axis[2]*omcosa+axis[1]*sina},
                    {axis[1]*axis[0]*omcosa+axis[2]*sina,
                     cosa+axis[1]*axis[1]*omcosa,
                     axis[1]*axis[2]*omcosa-axis[0]*sina},
                    {axis[2]*axis[0]*omcosa-axis[1]*sina,
                     axis[2]*axis[1]*omcosa+axis[0]*sina,
                     cosa+axis[2]*axis[2]*omcosa}};
  double result[3] = {0.0, 0.0, 0.0};
  for (i=0; i<3; ++i)
    for (j=0; j<3; ++j)
      result[i] += R[i][j]*x[j];
  for (i=0; i<3; ++i) x[i] = result[i];
}

int righthanded(double x[3], double y[3], double z[3])
/* Determines whether vectors x,y & z constitute a right-handed system */
/* by checking the sign of the triple product or det(x,y,z).           */
{
  double spatprodukt =   x[0]*y[1]*z[2] + y[0]*z[1]*x[2] + z[0]*x[1]*y[2]
                       - z[0]*y[1]*x[2] - x[0]*z[1]*y[2] - y[0]*x[1]*z[2];
  int result = (spatprodukt >= 0.0) ? 1 : 0;
  return result;
}

void orthoproject(double x[3], double vec1[3], double vec2[3])
/* Determines the orthogonal projection of vector x onto the span of */
/* the two ORTHONORMAL (!) vectors vec1 and vec2.                    */
{
  double sprod1 = x[0]*vec1[0] + x[1]*vec1[1] + x[2]*vec1[2];
  double sprod2 = x[0]*vec2[0] + x[1]*vec2[1] + x[2]*vec2[2];
  x[0] = sprod1*vec1[0] + sprod2*vec2[0];
  x[1] = sprod1*vec1[1] + sprod2*vec2[1];
  x[2] = sprod1*vec1[2] + sprod2*vec2[2];
}

double angle(double x[3], double y[3])
/* Determines the angle between vectors x & y. */
{
  double sprod = x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
  double absx  = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  double absy  = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
  return acos(sprod/(absx*absy));
}

void coord2vec(double sinlati, double longi, double x[3])
/* Turns geographical coordinates (latitude & longitude) into a vector -     */
/* Result is a unit (!) vector referring to the (right-handed) coordinate    */
/* system spanned by the three vectors pointing from geocenter to:           */
/*   x) intersection of greenwich meridian and equator                       */
/*   y) intersection of 90 deg. East meridian and equator                    */
/*   z) north pole                                                           */
{
  double coslati = sqrt(1.0-sinlati*sinlati);
  x[0] = cos(longi) * coslati;  /* `Greenwich'  */
  x[1] = sin(longi) * coslati;  /* `Ganges'     */
  x[2] = sinlati;               /* `North Pole' */
}

void vec2coord(double x[3], double *sinlati, double *longi)
/* Derive geographical coordinates from a vector (see also previous function). */
{
  double greenwich[3] = {1.0, 0.0, 0.0};
  double ganges[3]    = {0.0, 1.0, 0.0};
  double northpole[3] = {0.0, 0.0, 1.0};
  double dummy[3]     = {x[0], x[1], x[2]};
  *sinlati = sin(0.5*pi - angle(northpole, x));
  orthoproject(dummy, greenwich, ganges);
  *longi = angle(greenwich, dummy);
  if (righthanded(greenwich,northpole,dummy))
    *longi *= -1.0;
}







void setmcmcseed(struct runpar *run)
//If run->mcmcseed==0, set it using the system clock
{
  struct timeval time;
  struct timezone tz;
  gettimeofday(&time, &tz);
  if(run->mcmcseed==0) run->mcmcseed = time.tv_usec;
}


void setseed(int *seed)
//If seed==0, set (randomise) it using the system clock
{
  struct timeval time;
  struct timezone tz;
  gettimeofday(&time, &tz);
  if(*seed==0) *seed = time.tv_usec;
}


