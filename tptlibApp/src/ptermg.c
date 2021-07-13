#include <math.h>
#include <string.h>
#include <ctype.h>

void tptPtpre ( double, double, int, double[3][3], int, double[3],
                double*, double*, double*, double*, double*, double*,
                double*, double*, double*, double*, double* );
int tptPtappl ( char, double, double, double, int, double[3][3],
                int, int, double[3], double[8] );
int tptScomp ( char*, char* );

int tptPtermg ( char *cname, char *coform, double cval, double auxil[3],
                double sinphi, double cosphi, int mount,
                double ae2nm[3][3], int jbp, int method,
                double vw[3], double vc[3], double pmodel[7] )
/*
**  - - - - - - - - - -
**   t p t P t e r m g
**  - - - - - - - - - -
**
**  Compute the effect of one GENERIC term in the pointing model.
**
**  The effects of STANDARD and LOCALLY IMPLEMENTED terms are computed
**  in the routines tptPterms and tptPterml respectively.
**
**  Given (arguments):
**     cname    char[][9]    name of generic term (Note 1)
**     coform   char[][9]    format of generic term (supplied by tptTerm)
**     cval     double       coefficient value (radians)
**     auxil    double[3]    auxiliary readings if any
**     sinphi   double       sine of telescope latitude (Note 2)
**     cosphi   double       cosine of telescope latitude (Note 2)
**     mount    int          mount type (Note 3)
**     ae2nm    double[3][3] [Az,El] to nominal mount matrix (Note 3)
**     jbp      int          "below pole" flag (Note 4)
**     method   int          1 = add correction, else = subtract correction
**     vw       double[3]    working [-h,d] vector for calculations (Note 5)
**     vc       double[3]    current [-h,d] vector (Note 5)
**     pmodel   double[7]    basic pointing model (Note 6)
**
**  Returned (arguments):
**     vc       double[3]  updated
**     pmodel   double[7]  updated
**
**  Returned (function value):
**              int        status: -ve = error
**
**  Called:  tptPtpre, tptScomp, tptPtappl
**
**  Notes:
**
**  1) It is up to the caller to make sure that the term name "cname" is
**     valid.  If it is not, no adjustment is made and an error status
**     is returned.
**
**  2) The functions of latitude, sinphi and cosphi, are used as supplied,
**     without checks for range and consistency.
**
**  3) The "mount" argument specifies the type of mount:
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
**  4) The "jbp" flag specifies what state the mount is in:
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
**  5) The "vw" [-h,d] vector is the direction used to calculate the
**     correction.  The "vc" [-h,d] vector is the direction to which
**     the correction is being applied.
**
**  6) The array "pmodel" is a list of seven coefficients which
**     comprise a basic pointing model:
**
**                                               term names
**        element    meaning               equatorial     other
**
**         [0,1]   index errors              -IH,ID       IA,IE
**          [2]    vertical error             FLOP        FLOP
**          [3]    collimation error          -CH          CA
**          [4]    nonperpendicularity        -NP         NPAE
**         [5,6]   misalignment              -MA,ME       AW,AN
**
**     This 7-term basic model combines the contributions from all the
**     terms in the actual, operational, model.  The recipient, for
**     example a telescope control system, need implement only this
**     7-term model (once and for all, efficiently and with proper
**     handling of awkward cases) in order to correct for all the terms
**     in whatever operational model is in use at the time.  Note also
**     that the present routine need be invoked only occasionally, to
**     to keep track of the effects of the changing telescope position.
**     Between updates, the 7-term model will be an adequate summary of
**     the operational model for the current region of sky.
**
**  7) The purpose of the present routine is to compute the effect of
**     one generic pointing term at the current position.
**
**     To simplify the coding of the different terms, the routine makes a
**     preliminary call to the tptPtpre function in order to make available
**     the telescope position in four different systems:
**
**        spherical HA/Dec
**        Cartesian -HA/Dec
**        spherical Az/El
**        Cartesian Az/El
**
**     The correction for each term is implemented by making one (or at
**     most two) calls to the tptPtappl function.
**
**  8) The following sign conventions are used:
**
**        Azimuth - zero for due south, and 90 degrees for due east.
**        (Beware!  Not the surveyor's convention, the one used by
**        TPOINT for its external interfaces!)
**
**        Elevation - the horizon is elevation zero, and the zenith
**        is elevation 90 degrees.
**
**        HA and EW corrections - positive when the corrected
**        telescope position is west of the uncorrected position.
**
**        Dec corrections - positive when the corrected telescope
**        position is north of the uncorrected position.
**
**        Azimuth and LR corrections - positive when the
**        corrected telescope position is to the left of the
**        uncorrected position as seen by someone standing at
**        the telescope looking at the sky.
**
**        Zenith distance corrections - positive when the corrected
**        telescope position is below the uncorrected position.
**
**        Elevation corrections - positive when the corrected
**        telescope position is above the uncorrected position.
**
**     The terms are reversed in sign when adjusting star positions
**     to fit telescope positions, so that the coefficients stay
**     unchanged in sign whichever of the two options (star adjusted
**     to fit telescope, or telescope adjusted to fit star) is
**     being used.
**
**  9) This routine is designed to be usable in other applications,
**     with or without local enhancements.  It is self-contained, with
**     all information exchanged through the argument list;  no
**     non-standard header files are included.
**
**  Last revision:   15 September 2000
**
**  Copyright 2000 P.T.Wallace.  All rights reserved.
*/

{
   double xa, ya, za, az, el, zd, xe, ye, ze, ha, dec, w,
          s, c, q, w1, w2, ef;
   int i1, i2, i, phora, c1, c2, f1, f2;
   char type;


/*
** ---------------------------------------------
** Transform working Az/El vector to other forms
** ---------------------------------------------
*/

/* Spherical and Cartesian forms of Az/El and HA/Dec. */
   tptPtpre ( sinphi, cosphi, mount, ae2nm, jbp, vw,
              &xa, &ya, &za, &az, &el, &zd, &xe, &ye, &ze, &ha, &dec );

/* Parallactic angle. */
   s = sin ( ha ) * cosphi;
   c = cos ( dec ) * sinphi - cos ( ha ) * ze * cosphi;
   q = ( s != 0.0 || c != 0.0 ) ? atan2 ( s, c ) : 0.0;

/*
** ------------------------------------------
** Evaluate adjustment for this term in model
** ------------------------------------------
*/

/* Polynomial, harmonic or auxiliary? */
   phora = (int) coform [ 0 ];
   if ( phora == 'P' ) {


   /* ---------- */
   /* Polynomial */
   /* ---------- */

   /* Extract adjustment type. */
      type = cname [ 1 ];

   /* Extract elements for each format separately. */
      if ( ! tptScomp ( coform, "PRCICI" ) ) {
         c1 = (int) cname [ 2 ];
         i1 = ( (int) cname [ 3 ] ) - '0';
         c2 = (int) cname [ 4 ];
         i2 = ( (int) cname [ 5 ] ) - '0';
      } else if ( ! tptScomp ( coform, "PRCIC" ) ) {
         c1 = (int) cname [ 2 ];
         i1 = ( (int) cname [ 3 ] ) - '0';
         c2 = cname [ 4 ];
         i2 = 1;
      } else if ( ! tptScomp ( coform, "PRCI" ) ) {
         c1 = (int) cname [ 2 ];
         i1 = ( (int) cname [ 3 ] ) - '0';
         c2 = 'A';
         i2 = 0;
      } else if ( ! tptScomp ( coform, "PRC" ) ) {
         c1 = (int) cname [ 2 ];
         i1 = 1;
         c2 = 'A';
         i2 = 0;
      } else {

      /* Illegal format (impossible). */
         c1 = '?';
         i1 = 1;
         c2 = '?';
         i2 = 1;
      }

   /* First coordinate. */
      if ( c1 == 'H' ) {
         w = ha;
      } else if ( c1 == 'D' ) {
         w = dec;
      } else if ( c1 == 'A' ) {
         w = az;
      } else if ( c1 == 'Z' ) {
         w = zd;
      } else if ( c1 == 'E' ) {
         w = el;
      } else if ( c1 == 'Q' ) {
         w = q;
      } else {
         w = 0.0;
      }

   /* Raise to nominated power. */
      for ( w1 = 1.0, i = 0; i < i1; w1 *= w, i++ );

   /* Second coordinate. */
      if ( c2 == 'H' ) {
         w = ha;
      } else if ( c2 == 'D' ) {
         w = dec;
      } else if ( c2 == 'A' ) {
         w = az;
      } else if ( c2 == 'Z' ) {
         w = zd;
      } else if ( c2 == 'E' ) {
         w = el;
      } else if ( c2 == 'Q' ) {
         w = q;
      } else {
         w = 0.0;
      }

   /* Raise to nominated power. */
      for ( w2 = 1.0, i = 0; i < i2; w2 *= w, i++ );

   } else if ( phora == 'H' ) {

   /* -------- */
   /* Harmonic */
   /* -------- */

   /* Extract adjustment type. */
      type = cname [ 1 ];

   /* Extract elements for each format separately. */
      if ( ! tptScomp ( coform, "HRFCIFCI" ) ) {
         f1 = (int) cname [ 2 ];
         c1 = (int) cname [ 3 ];
         i1 = ( (int) cname [ 4 ] ) - '0';
         f2 = (int) cname [ 5 ];
         c2 = (int) cname [ 6 ];
         i2 = ( (int) cname [ 7 ] ) - '0';
      } else if ( ! tptScomp ( coform, "HRFCIFC" ) ) {
         f1 = (int) cname [ 2 ];
         c1 = (int) cname [ 3 ];
         i1 = ( (int) cname [ 4 ] ) - '0';
         f2 = (int) cname [ 5 ];
         c2 = (int) cname [ 6 ];
         i2 = 1;
      } else if ( ! tptScomp ( coform, "HRFCFCI" ) ) {
         f1 = (int) cname [ 2 ];
         c1 = (int) cname [ 3 ];
         i1 = 1;
         f2 = (int) cname [ 4 ];
         c2 = (int) cname [ 5 ];
         i2 = ( (int) cname [ 6 ] ) - '0';
      } else if ( ! tptScomp ( coform, "HRFCFC" ) ) {
         f1 = (int) cname [ 2 ];
         c1 = (int) cname [ 3 ];
         i1 = 1;
         f2 = (int) cname [ 4 ];
         c2 = (int) cname [ 5 ];
         i2 = 1;
      } else if ( ! tptScomp ( coform, "HRFCIII" ) ) {
         f1 = (int) cname [ 2 ];
         c1 = (int) cname [ 3 ];
         i1 = 100 * ( (int) cname [ 4 ] - '0' ) +
               10 * ( (int) cname [ 5 ] - '0' ) +
                    ( (int) cname [ 6 ] - '0' );
         f2 = 'C';
         c2 = 'H';
         i2 = 0;
      } else if ( ! tptScomp ( coform, "HRFCII" ) ) {
         f1 = cname [ 2 ];
         c1 = cname [ 3 ];
         i1 = 10 * ( (int) cname [ 4 ] - '0' ) +
                   ( (int) cname [ 5 ] - '0' );
         f2 = 'C';
         c2 = 'H';
         i2 = 0;
      } else if ( ! tptScomp ( coform, "HRFCI" ) ) {
         f1 = cname [ 2 ];
         c1 = cname [ 3 ];
         i1 = ( (int) cname [ 4 ] ) - '0';
         f2 = 'C';
         c2 = 'H';
         i2 = 0;
      } else if ( ! tptScomp ( coform, "HRFC" ) ) {
         f1 = cname [ 2 ];
         c1 = cname [ 3 ];
         i1 = 1;
         f2 = 'C';
         c2 = 'H';
         i2 = 0;
      } else if ( ! tptScomp ( coform, "HRFCIIFC" ) ) {
         f1 = cname [ 2 ];
         c1 = cname [ 3 ];
         i1 = 10 * ( (int) cname [ 4 ] - '0' ) +
                   ( (int) cname [ 5 ] - '0' );
         f2 = cname [ 6 ];
         c2 = cname [ 7 ];
         i2 = 1;
      } else if ( ! tptScomp ( coform, "HRFCFCII" ) ) {
         f1 = cname [ 2 ];
         c1 = cname [ 3 ];
         i1 = 1;
         f2 = cname [ 4 ];
         c2 = cname [ 5 ];
         i2 = 10 * ( (int) cname [ 6 ] - '0' ) +
                   ( (int) cname [ 7 ] - '0' );
      } else {

      /* Illegal format (impossible). */
         f1 = 'S';
         c1 = '?';
         i1 = 0;
         f2 = 'S';
         c2 = '?';
         i2 = 0;
      }

   /* First coordinate. */
      if ( c1 == 'H' ) {
         w = ha;
      } else if ( c1 == 'D' ) {
         w = dec;
      } else if ( c1 == 'A' ) {
         w = az;
      } else if ( c1 == 'Z' ) {
         w = zd;
      } else if ( c1 == 'E' ) {
         w = el;
      } else if ( c1 == 'Q' ) {
         w = q;
      } else {
         w = 0.0;
      }

   /* Multiple. */
      w *= (double) i1;

   /* Sin or cos. */
      if ( f1 == 'S' ) {
         w1 = sin ( w );
      } else {
         w1 = cos ( w );
      }

   /* Second coordinate. */
      if ( c2 == 'H' ) {
         w = ha;
      } else if ( c2 == 'D' ) {
         w = dec;
      } else if ( c2 == 'A' ) {
         w = az;
      } else if ( c2 == 'Z' ) {
         w = zd;
      } else if ( c2 == 'E' ) {
         w = el;
      } else if ( c2 == 'Q' ) {
         w = q;
      } else {
         w = 0.0;
      }

   /* Multiple. */
      w *= (double) i2;

   /* Sin or cos. */
      if ( f2 == 'S' ) {
         w2 = sin ( w );
      } else {
         w2 = cos ( w );
      }

   } else if ( phora == 'A' ) {

   /* --------- */
   /* Auxiliary */
   /* --------- */

   /* Extract adjustment type. */
      type = cname [ 2 ];

   /* Pick the specified auxiliary reading. */
      w1 = auxil [ (int) cname [ 1 ] - '1' ];
      w2 = 1.0;

   } else {

   /* Illegal format (impossible). */
      return 1;
   }

/* Adjustment function value. */
   ef = w1 * w2;

/*
** --------------------
** Apply the adjustment
** --------------------
*/

   return tptPtappl ( type, ef * cval, sinphi, cosphi, mount, ae2nm,
                      jbp, method, vc, pmodel );
}
