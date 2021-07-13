#include <math.h>

void tptPtpre ( double sinphi, double cosphi, int mount,
                double ae2nm[3][3], int jbp, double vw[3],
                double *xa, double *ya, double *za,
                double *az, double *el, double *zd,
                double *xe, double *ye, double *ze,
                double *ha, double *dec )
/*
**  - - - - - - - - -
**   t p t P t p r e
**  - - - - - - - - -
**
**  As a preliminary to evaluating a pointing term, transform the
**  working vector into various useful forms.
**
**  Given (arguments):
**     sinphi    double       sine of telescope latitude (Note 1)
**     cosphi    double       cosine of telescope latitude (Note 1)
**     mount     int          mount type (Note 2)
**     ae2nm     double[3][3] [Az,El] to nominal mount matrix (Note 2)
**     jbp       int          "below pole" flag (Note 3)
**     vw        double[3]    [-h,d] vector for calculations (Note 4)
**
**  Returned (arguments):
**     xa,ya,za  double      Az/El direction cosines (Note 2)
**     az        double      mount azimuth (radians, S=0,E=90; Note 2)
**     el        double      mount elevation (radians; Note 2)
**     zd        double      zenith distance (radians)
**     xe,ye,ze  double      -HA/Dec direction cosines
**     ha        double      hour angle (radians)
**     dec       double      declination (radians)
**
**  Notes:
**
**  1) The functions of latitude, sinphi and cosphi, are used as supplied,
**     without checks for range and consistency.
**
**  2) The "mount" argument specifies the type of mount:
**
**     * If mount=1, the mount is an equatorial and the basic model is
**       -IH, ID, FLOP, -CH, -NP, MA, ME.
**
**     * If mount=2, the mount is an altazimuth and the basic model is
**       IA, IE, FLOP, CA, NPAE, AW, AN.
**
**     * If mount=3, the mount is a generalized gimbal and the basic
**       model is IA, IE, FLOP, CA, NPAE, AW, AN.  In this case, all
**       the terms except FLOP are in the mount's nominal frame, as given
**       by the matrix ae2nm.
**
**     For a generalized gimbal mount, the arguments az and el receive
**     the mount roll and pitch rather than the true azimuth and
**     elevation.  The zd argument is always the true ZD.
**
**  3) The "jbp" flag specifies what state the mount is in:
**
**     * If jbp=0, the mount is in the normal state where the mechanical
**       declination or elevation is in the range +/- 90 deg.
**
**     * If jbp=1, the mount is in the other state, where the mechanical
**       declination or elevation is outside the range +/- 90 deg and the
**       mechanical hour angle or azimuth is 180 deg away from the
**       corresponding celestial value.  This state is variously called
**       "below the pole", "other side of the pier", "beyond the zenith"
**       etc.
**
**  4) The "vw" [-h,d] vector is the direction used to calculate the
**     correction.
**
**  Last revision:   8 November 1999
**
**  Copyright 1999 P.T.Wallace.  All rights reserved.
*/

#define DPI 3.1415926535897932384626433832795028841971693993751
#define DPI2 1.5707963267948966192313216916397514420985846996876
#define EQUAT 1
#define GIMBAL 3

{
   double xe_, ye_, ze_, r, ha_, dec_, zd_, xa_, ya_, za_, az_, el_,
          x, y, z;



/* Working [-HA,Dec] direction cosines. */
   xe_ = vw[0];
   ye_ = vw[1];
   ze_ = vw[2];

/* To Spherical. */
   r = sqrt ( xe_ * xe_ + ye_ * ye_ );
   ha_ = ( r != 0.0 ) ? atan2 ( -ye_, xe_ ) : 0.0;
   dec_ = ( ze_ != 0.0 ) ? atan2 ( ze_, r ) : 0.0;

/* Rotate into [Az,El]. */
   xa_ = sinphi * xe_ - cosphi * ze_;
   ya_ = ye_;
   za_ =  cosphi * xe_ + sinphi * ze_;

/* To spherical. */
   r = sqrt ( xa_ * xa_ + ya_ * ya_ );
   az_ = ( r != 0.0 ) ? atan2 ( ya_, xa_ ) : 0.0;
   el_ = ( za_ != 0.0 ) ? atan2 ( za_, r ) : 0.0;

/* Zenith distance. */
   zd_ = DPI2 - el_;

/* For a generalized gimbal, substitute nominal mount [roll,pitch]. */
   if ( mount == GIMBAL ) {
      x = xa_;
      y = ya_;
      z = za_;
      xa_ = ae2nm[0][0] * x + ae2nm[0][1] * y + ae2nm[0][2] * z;
      ya_ = ae2nm[1][0] * x + ae2nm[1][1] * y + ae2nm[1][2] * z;
      za_ = ae2nm[2][0] * x + ae2nm[2][1] * y + ae2nm[2][2] * z;
      r = sqrt ( xa_ * xa_ + ya_ * ya_ );
      az_ = ( r != 0.0 ) ? atan2 ( ya_, xa_ ) : 0.0;
      el_ = ( za_ != 0.0 ) ? atan2 ( za_, r ) : 0.0;
   }

/* If "below the pole", use the mechanical HA/Dec or Az/El. */
   if ( jbp ) {
      if ( mount == EQUAT ) {
         ha_ += DPI;
         dec_ = ( ( sinphi < 0.0 ) ? -DPI : DPI ) - dec_;
      } else {
         az_ += DPI;
         el_ = DPI - el_;
      }
   }

/* Return the results. */
   *xa = xa_;
   *ya = ya_;
   *za = za_;
   *az = az_;
   *el = el_;
   *zd = zd_;
   *xe = xe_;
   *ye = ye_;
   *ze = ze_;
   *ha = ha_;
   *dec = dec_;

   return;
}
