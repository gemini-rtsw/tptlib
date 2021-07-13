#include <math.h>

void tptPtpre ( double, double, int, double[3][3], int, double[3],
                double*, double*, double*, double*, double*, double*,
                double*, double*, double*, double*, double* );
int tptPtappl ( char, double, double, double, int, int, int,
                double[3], double[8] );

int tptPterml ( int nterm, double cval, double auxil[3],
                double sinphi, double cosphi, int mount,
                double ae2nm[3][3], int jbp, int method,
                double vw[3], double vc[3], double pmodel[7] )
/*
**  - - - - - - - - - -
**   t p t P t e r m l
**  - - - - - - - - - -
**
**  Compute the effect of one LOCALLY IMPLEMENTED term in the pointing
**  model.
**
**  The effects of STANDARD and GENERIC terms are computed in the
**  routines tptPterms and tptpPtermg respectively.
**
**  n.b.  Changes in this routine must be coordinated with changes
**        in the routine tptTrmloc.
**
**  *** THIS IS AN EXAMPLE PROGRAM, NEVER CALLED ***
**  (because the supplied tptTrmloc doesn't define any term names)
**
**  Given (arguments):
**     nterm    int          coefficient number within local set (Note 1)
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
**  Called:  tptPtpre, tptPtappl
**
**  Notes:
**
**  1) It is up to the caller to make sure that the argument "nterm" is
**     within the allowed range, i.e. 1 to the number of local terms
**     implemented.  If it is not, no adjustment is made and an error
**     status is returned.
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
**     one local pointing term at the current position.
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
**  9) This routine is analogous to the tptPterms routine, but is
**     intended for the implementation of local terms.  It has an
**     identical argument list, through which everything passes, and
**     includes no non-standard header files.
**
**  Last revision:   15 September 2000
**
**  Copyright 2000 P.T.Wallace.  All rights reserved.
*/

#define TINY 1e-20

{
#if 0
   double xa, ya, za, az, el, zd, xe, ye, ze, ha, dec, ef1, ef2, w,
          h2, h3, h5, d2, sh, s2h, s3h, s4h, h12, d3, d6, ch;
   char type1, type2;
   int j;



/*
** ---------------------------------------------
** Transform working Az/El vector to other forms
** ---------------------------------------------
*/


   tptPtpre ( sinphi, cosphi, mount, jbp, vw,
              &xa, &ya, &za, &az, &el, &zd, &xe, &ye, &ze, &ha, &dec );

/*
** ------------------------------------------
** Evaluate adjustment for this term in model
** ------------------------------------------
*/

/* Disable second component by default. */
   type2 = (char) '\0';

/* Deal with the the specified term. */
   switch ( nterm ) {

/*--------------------------------------------------------------------
**
** BEND: example affecting just HA
**
** Change in HA equal to coefficient BEND.
**
** In equatorial spherical coordinates:
**    dHA  = +BEND
**    dDec = 0
*/
   case 1:
      ef1 = 1.0;
      type1 = 'H';
      break;

/*--------------------------------------------------------------------
**
** SLIP: example affecting just Dec
**
** Change in Dec equal to coefficient SLIP.
**
** In equatorial spherical coordinates:
**    dHA  = 0
**    dDec = +SLIP
*/
   case 2:
      ef1 = 1.0;
      type1 = 'D';
      break;

/*--------------------------------------------------------------------
**
** JUDDER:  example using auxiliary reading
**
** Change in HA equal to coefficient JUDDER multiplied by first
** auxiliary reading.
**
** In equatorial spherical coordinates:
**    dHA  = +JUDDER*AUX1
**    dDec = 0
*/
   case 3:
      ef1 = auxil[0];
      type1 = 'H';
      break;

/*--------------------------------------------------------------------
**
** SHAKE: example using rotation of the coordinate system
**
** Rotation about EW axis equal to coefficient SHAKE.
**
** In equatorial spherical coordinates (approximately):
**    dHA  = +SHAKE*sin(HA)*tan(DEC)
**    dDec = +SHAKE*cos(HA)
*/
   case 4:
      ef1 = 1.0;
      type1 = 'N';
      break;

/*--------------------------------------------------------------------
**
** RATTLE: example affecting just elevation
**
** Change in ZD equal to coefficient RATTLE multiplied by sin(ZD).
**
** In altazimuth spherical coordinates:
**    dAz = 0
**    dEl = +RATTLE*cos(EL)
**
** In equatorial spherical coordinates (approximately):
**    dHA  = -RATTLE*cos(PHI)*sin(HA)/cos(DEC)
**    dDec = -RATTLE*(cos(PHI)*cos(HA)*sin(DEC)-sin(PHI)*cos(DEC))
*/
   case 5:
      ef1 = sin ( zd );
      type1 = 'Z';
      break;

/*--------------------------------------------------------------------
**
** ROLL: example using EW and NS corrections
**
** Change EW on the sky equal to coefficient ROLL multiplied by
** -sin(Dec), and change in Dec equal to the same coefficient
** multiplied by cos(Dec).
**
** In equatorial spherical coordinates (approximately):
**    dHA  = -ROLL*tan(DEC)
**    dDec = +ROLL*cos(DEC)
*/
   case 6:
      ef1 = - ze;
      type1 = 'X';
      ef2 = sin ( dec );
      type2 = 'D';
      break;

/*--------------------------------------------------------------------
**
** CRUNCH: example using sideways and up/down corrections
**
** Change LR on the sky equal to coefficient GRAUNCH multiplied
** by sin(El), and change in El equal to the same coefficient
** multiplied by cos(El).
**
** In altazimuth spherical coordinates (approximately):
**    dAz = +CRUNCH*tan(EL)
**    dEl = +CRUNCH*cos(EL)
*/
   case 7:
      ef1 = ze;
      type1 = 'S';
      ef2 = cos ( el );
      type2 = 'E';
      break;

/*--------------------------------------------------------------------
**
** SLOUCH: just another example
**
** Change in Az equal to coefficient SLOUCH.
**
** In Altazimuth spherical coordinates:
**    dAz = +SLOUCH
**    dEl = 0
*/
   case 8:
      ef1 = 1.0;
      type1 = 'A';
      break;

/*--------------------------------------------------------------------
**
** WANDER: just another example
**
** Change in HA equal to coefficient WANDER multiplied by sin(HA).
**
** In equatorial spherical coordinates:
**    dHA  = +WANDER*sin(HA)
**    dDec = 0
*/
   case 9:
      ef1 = sin ( ha );
      type1 = 'H';
      break;

/*--------------------------------------------------------------------
**
** WOBBLE: just another example
**
** Change in HA/Dec nonperpendicularity equal to coefficient WOBBLE
** multiplied by cos(HA).
**
** In equatorial spherical coordinates (approximately):
**    dHA  = +WOBBLE*cos(HA)*tan(DEC)
**    dDec = 0
*/
   case 10:
      ef1 = cos ( ha );
      type1 = 'P';
      break;

/*--------------------------------------------------------------------
**
** DRIFT: just another example
**
** Change in left-right polar axis misalignment equal to coefficient
** DRIFT multiplied by cos(Dec).
**
** In equatorial spherical coordinates (approximately):
**    dHA  = -DRIFT*cos(HA)*sin(DEC)
**    dDec = +DRIFT*sin(HA)*cos(DEC)
*/
   case 12:
      ef1 = cos ( dec );
      type1 = 'L';
      break;

/*--------------------------------------------------------------------*/

/* Unrecognized term number. */
   default:
      return -1;
   }

/*
** --------------------
** Apply the adjustment
** --------------------
*/


/* The first component. */
   if ( ( j = tptPtappl ( type1, ef1 * cval,
                          sinphi, cosphi, mount, ae2nm, jbp, method,
                          vc, pmodel ) ) ) return j;

/* The second component, if any. */
   if ( (int) type2 ) {
      if ( ( j = tptPtappl ( type2, ef2 * cval,
                             sinphi, cosphi, mount, ae2nm, jbp, method,
                             vc, pmodel ) ) ) return j;
   }
#endif

/* Success. */
   return 0;
}
