#include <math.h>

static void aiaie ( double, double, int, double[3][3],
                    double[3], double[3], double*, double* );
static void aihid ( double[3], double[3], double*, double* );
static void xyz2ab ( double, double, double, double*, double* );
static void ab2xyz ( double, double, double*, double*, double* );
static void ve2a ( double, double, int, double[3][3],
                   double, double, double, double*, double*, double* );
static void va2e ( double, double, int, double[3][3],
                   double, double, double, double[3] );

int tptPtappl ( char atype, double asize,
                double sinphi, double cosphi, int mount, double ae2nm[3][3],
                int jbp, int method, double vc[3], double pmodel[7] )
/*
**  - - - - - - - - - -
**   t p t P t a p p l
**  - - - - - - - - - -
**
**  Apply the effect of one term in the pointing model.
**
**  Given:
**     atype    int          the form of the adjustment (Note 1)
**     asize    double       the size of the adjustment (radians)
**     sinphi   double       sine of telescope latitude (Note 2)
**     cosphi   double       cosine of telescope latitude (Note 2)
**     mount    int          mount type (Note 3)
**     ae2nm    double[3][3] [Az,El] to nominal mount matrix (Note 3)
**     jbp      int          "below pole" flag (Note 4)
**     method   int          1 = add correction, else = subtract correction
**     vc       double[3]    current [-h,d] vector (Note 5)
**     pmodel   double[7]    basic pointing model (Note 6)
**
**  Returned (arguments):
**     vc       double[3]    updated
**     pmodel   double[7]    updated
**
**  Returned (function value):
**              int          status: -ve = error
**
**  Notes:
**
**  1) The argument "atype" is a character that selects in what form the
**     adjustment "asize" is to be applied:
**
**          atype    what is adjusted
**
**           'H'     hour angle
**           'X'     east-west collimation
**           'D'     declination
**           'U'     polar axis up/down
**           'L'     polar axis left/right
**           'P'     HA/Dec nonperpendicularity
**           'A'     mount azimuth
**           'S'     left-right collimation
**           'Z'     zenith distance (sky)
**           'E'     mount elevation
**           'N'     azimuth axis north/south
**           'W'     azimuth axis east/west
**           'V'     Az/El nonperpendicularity
**
**     For a generalized gimbal mount, "azimuth" and "elevation" mean
**     mount roll and pitch respectively.  Zenith distance always means
**     with respect to the zenith, not with respect to the mount pole.
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
**  5) The "vc" vector is the [-h,d] direction to which the correction is
**     being applied.
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
**  7) This routine, which employs no non-standard #includes, is suitable
**     for use in applications other than the TPOINT pointing-analysis
**     application itself, for example a telescope control system.
**
**  Called:  aiaie, xyz2ab, ab2xy
**
**  Last revision:   15 September 2000
**
**  Copyright 2000 P.T.Wallace.  All rights reserved.
*/

#define DPI    3.1415926535897932384626433832795028841971693993751
#define D2PI   6.2831853071795864769252867665590057683943387987502
#define EQUAT  1
#define GIMBAL 3
#define TINY   1e-20

/* Basic model, equatorial */
#define MIH  pmodel[0]
#define ID   pmodel[1]
#define FLOP pmodel[2]
#define MCH  pmodel[3]
#define MNP  pmodel[4]
#define MMA  pmodel[5]
#define ME   pmodel[6]

/* Basic model, altazimuth */
#define IA   pmodel[0]
#define IE   pmodel[1]
/*      FLOP pmodel[2] */
#define CA   pmodel[3]
#define NPAE pmodel[4]
#define AW   pmodel[5]
#define AN   pmodel[6]

{
   int isgim, i;
   double xold, yold, zold, adj, s, c, vnew[3], aold, bold, anew, bnew,
          r, xnew, ynew, znew, x, y, z, xa, ya, za;



/* Set the "generalized gimbal" flag appropriately. */
   isgim = ( mount == GIMBAL );

/* Copy the unadjusted [-HA,Dec] vector. */
   xold = vc[0];
   yold = vc[1];
   zold = vc[2];

/* Select the specified action. */
   switch ( (int) atype ) {

/*
** ------------------------
** Adjustment to hour angle
** ------------------------
*/

   case 'H':
   case 'h':

   /* Apply the adjustment to the given vector */
   /* ---------------------------------------- */

   /* The adjustment to -IH. */
      adj = method == 1 ? -asize : asize;

   /* Apply the adjustment. */
      s = sin ( adj );
      c = cos ( adj );
      vnew[0] = xold * c - yold * s;
      vnew[1] = xold * s + yold * c;
      vnew[2] = zold;

   /* Also apply the adjustment to the model */
   /* -------------------------------------- */

      if ( mount == EQUAT) {

      /* Equatorial case: apply the adjustment directly to -IH. */
         MIH -= asize;

      } else {

      /* Altaz/gimbal case: apply corresponding adjustments to IA and IE. */
         aiaie ( sinphi, cosphi, isgim, ae2nm, vc, vnew, &IA, &IE );
      }

      break;

/*
** -------------------------
** Adjustment to declination
** -------------------------
*/

   case 'D':
   case 'd':

   /* Apply the adjustment to the given vector */
   /* ---------------------------------------- */

   /* To spherical. */
      xyz2ab ( xold, yold, zold, &aold, &bold );

   /* The adjustment to ID. */
      adj = method == 1 ? asize : -asize;
      if ( jbp ) adj = -adj;

   /* Apply the adjustment. */
      anew = aold;
      bnew = bold + adj;

   /* The adjusted [-HA,Dec] vector. */
      ab2xyz ( anew, bnew, vnew, vnew+1, vnew+2 );

   /* Also apply the adjustment to the model */
   /* -------------------------------------- */

      if ( mount == EQUAT) {

      /* Equatorial case: apply the adjustment directly to ID. */
         ID += asize;

      } else {

      /* Altaz/gimbal case: apply corresponding adjustments to IA and IE. */
         aiaie ( sinphi, cosphi, isgim, ae2nm, vc, vnew, &IA, &IE );
      }

      break;

/*
** -----------------------------------------------------------------
** Adjustment to east-west collimation or HA/Dec nonperpendicularity
** -----------------------------------------------------------------
*/

   case 'X':
   case 'x':
   case 'P':
   case 'p':

   /* Apply the adjustment to the given vector */
   /* ---------------------------------------- */

   /* Functions of -HA. */
      r = sqrt ( xold * xold + yold * yold );
      if ( r > TINY ) {
         s = yold / r;
         c = xold / r;
      } else {
         s = 0.0;
         c = 1.0;
      }

   /* The adjustment to -CH or -NP. */
      adj = method == 1 ? -asize : asize;
      if ( jbp ) adj = -adj;

   /* Convert a -NP adjustment to -CH. */
      if ( (int) atype == 'P' || (int) atype == 'p' ) adj *= zold;

   /* Apply the adjustment. */
      vnew[0] = xnew = xold - adj * s;
      vnew[1] = ynew = yold + adj * c;
      vnew[2] = znew = zold;

   /* The adjusted and normalized [-HA,Dec] vector. */
      r = sqrt ( xnew * xnew + ynew * ynew + znew * znew );
      if ( r > TINY ) {
         vnew[0] /= r;
         vnew[1] /= r;
         vnew[2] /= r;
      }

   /* Also apply the adjustment to the model */
   /* -------------------------------------- */

      if ( mount == EQUAT) {

      /* Equatorial case: apply the adjustment directly
                                             to minus CH or minus NP. */
         if ( (int) atype == 'X' || (int) atype == 'x' ) {
            MCH -= asize;
         } else {
            MNP -= asize;
         }

      } else {

      /* Altaz/gimbal case: apply corresponding adjustments to IA and IE. */
         aiaie ( sinphi, cosphi, isgim, ae2nm, vc, vnew, &IA, &IE );
      }

      break;

/*
** --------------------------------
** Adjustment up-down to polar axis
** --------------------------------
*/

   case 'U':
   case 'u':

   /* Apply the adjustment to the given vector */
   /* ---------------------------------------- */

   /* The adjustment to ME. */
      adj = method ==1 ? asize : -asize;

   /* Apply the adjustment. */
      s = sin ( adj );
      c = cos ( adj );
      vnew[0] = xold * c - zold * s;
      vnew[1] = yold;
      vnew[2] = xold * s + zold * c;

   /* Also apply the adjustment to the model */
   /* -------------------------------------- */

      if ( mount == EQUAT) {

      /* Equatorial case: apply the adjustment directly to ME. */
         ME += asize;

      } else {

      /* Altaz/gimbal case: apply corresponding adjustments to IA and IE. */
         aiaie ( sinphi, cosphi, isgim, ae2nm, vc, vnew, &IA, &IE );
      }

      break;

/*
** -----------------------------------
** Left-right adjustment to polar axis
** -----------------------------------
*/

   case 'L':
   case 'l':

   /* Apply the adjustment to the given vector */
   /* ---------------------------------------- */

   /* The adjustment to MMA. */
      adj = method == 1 ? -asize : asize;

   /* Apply the adjustment. */
      s = sin ( adj );
      c = cos ( adj );
      vnew[0] = xold;
      vnew[1] = yold * c + zold * s;
      vnew[2] = - yold * s + zold * c;

   /* Also apply the adjustment to the model */
   /* -------------------------------------- */

      if ( mount == EQUAT) {

      /* Equatorial case: apply the adjustment directly to minus MA. */
         MMA -= asize;

      } else {

      /* Altaz/gimbal case: apply corresponding adjustments to IA and IE. */
         aiaie ( sinphi, cosphi, isgim, ae2nm, vc, vnew, &IA, &IE );
      }

      break;

/*
** -----------------------------
** Adjustment to zenith distance
** -----------------------------
**
** This means sky ZD even in the GIMBAL case.
*/

   case 'Z':
   case 'z':

   /* Apply the adjustment to the given vector */
   /* ---------------------------------------- */

   /* Rotate the given [-HA,Dec] vector into [Az,El] (sky, not mount). */
      ve2a ( sinphi, cosphi, 0, ae2nm, xold, yold, zold, &xa, &ya, &za );

   /* To spherical. */
      xyz2ab ( xa, ya, za, &aold, &bold );

   /* The adjustment. */
      adj = method == 1 ? -asize : asize;

   /* Apply the adjustment. */
      anew = aold;
      bnew = bold + adj;

   /* Back to Cartesian. */
      ab2xyz ( anew, bnew, &xa, &ya, &za );

   /* Return the adjusted [-HA,Dec] vector. */
      va2e ( sinphi, cosphi, 0, ae2nm, xa, ya, za, vnew );

   /* Also apply the adjustment to the model */
   /* -------------------------------------- */

   /* All cases: apply the adjustment directly to FLOP. */
      FLOP += asize;

      break;

/*
** ---------------------
** Adjustment to azimuth
** ---------------------
**
** In the GIMBAL case, "azimuth" means the mount's longitude coordinate.
*/

   case 'A':
   case 'a':

   /* Apply the adjustment to the given vector */
   /* ---------------------------------------- */

   /* Rotate the given [-HA,Dec] vector into [Az,El]. */
      ve2a ( sinphi, cosphi, isgim, ae2nm, xold, yold, zold, &xa, &ya, &za );

   /* The adjustment to IA. */
      adj = method == 1 ? asize : -asize;

   /* Apply the adjustment. */
      s = sin ( adj );
      c = cos ( adj );
      x = xa;
      y = ya;
      xa = x * c - y * s;
      ya = x * s + y * c;

   /* Return the adjusted [-HA,Dec] vector. */
      va2e ( sinphi, cosphi, isgim, ae2nm, xa, ya, za, vnew );

   /* Also apply the adjustment to the model */
   /* -------------------------------------- */

      if ( mount == EQUAT) {

      /* Equatorial case: apply corresponding adjustments to -IH and ID. */
         aihid ( vc, vnew, &MIH, &ID );

      } else {

      /* Altaz/gimbal case: apply the adjustment directly to IA. */
         IA += asize;
      }

      break;

/*
** -----------------------
** Adjustment to elevation
** -----------------------
**
** In the GIMBAL case, "elevation" means the mount's latitude coordinate.
*/

   case 'E':
   case 'e':

   /* Apply the adjustment to the given vector */
   /* ---------------------------------------- */

   /* Rotate the given [-HA,Dec] vector into [Az,El]. */
      ve2a ( sinphi, cosphi, isgim, ae2nm, xold, yold, zold, &xa, &ya, &za );

   /* To spherical. */
      xyz2ab ( xa, ya, za, &aold, &bold );

   /* The adjustment. */
      adj = method == 1 ? asize : -asize;
      if ( jbp ) adj = -adj;

   /* Apply the adjustment. */
      anew = aold;
      bnew = bold + adj;

   /* Back to Cartesian. */
      ab2xyz ( anew, bnew, &xa, &ya, &za );

   /* Return the adjusted [-HA,Dec] vector. */
      va2e ( sinphi, cosphi, isgim, ae2nm, xa, ya, za, vnew );

   /* Also apply the adjustment to the model */
   /* -------------------------------------- */

      if ( mount == EQUAT) {

      /* Equatorial case: apply corresponding adjustments to -IH and ID. */
         aihid ( vc, vnew, &MIH, &ID );

      } else {

      /* Altaz/gimbal case: apply the adjustment directly to IE. */
         IE += asize;
      }

      break;

/*
** -----------------------------------------------------------------
** Adjustment to left-right collimation or Az/El nonperpendicularity
** -----------------------------------------------------------------
**
** In the GIMBAL case, the Az/El axes are the mount axes themselves.
*/

   case 'S':
   case 's':
   case 'V':
   case 'v':

   /* Apply the adjustment to the given vector */
   /* ---------------------------------------- */

   /* Rotate the given [-HA,Dec] vector into [Az,El]. */
      ve2a ( sinphi, cosphi, isgim, ae2nm, xold, yold, zold, &xa, &ya, &za );

   /* Functions of azimuth (S=0,E=90). */
      r = sqrt ( xa * xa + ya * ya );
      if ( r > TINY ) {
         s = ya / r;
         c = xa / r;
      } else {
         s = 0.0;
         c = 1.0;
      }

   /* The adjustment to CA or NPAE. */
      adj = method == 1 ? asize : -asize;
      if ( jbp ) adj = -adj;

   /* Convert an NPAE adjustment to CA. */
      if ( (int) atype == 'V' || (int) atype == 'v' ) adj *= za;

   /* Apply the adjustment. */
      xa -= adj * s;
      ya += adj * c;

   /* Re-normalize. */
      r = sqrt ( xa * xa + ya * ya + za * za );
      if ( r > TINY ) {
         xa /= r;
         ya /= r;
         za /= r;
      }

   /* Return the adjusted [-HA,Dec] vector. */
      va2e ( sinphi, cosphi, isgim, ae2nm, xa, ya, za, vnew );

   /* Also apply the adjustment to the model */
   /* -------------------------------------- */

      if ( mount == EQUAT) {

      /* Equatorial case: apply corresponding adjustments to -IH and ID. */
         aihid ( vc, vnew, &MIH, &ID );

      } else {

      /* Altaz/gimbal case: apply the adjustment directly to CA or NPAE. */
         if ( (int) atype == 'S' || (int) atype == 's' ) {
            CA += asize;
         } else {
            NPAE += asize;
         }
      }

      break;

/*
** --------------------------------------
** Adjustment north-south to azimuth axis
** --------------------------------------
**
** In the GIMBAL case, "azimuth axis" means the mount's longitude axis,
** and "north south" means along the mount's prime meridian.
*/

   case 'N':
   case 'n':

   /* Apply the adjustment to the given vector */
   /* ---------------------------------------- */

   /* Rotate the given [-HA,Dec] vector into [Az,El]. */
      ve2a ( sinphi, cosphi, isgim, ae2nm, xold, yold, zold, &xa, &ya, &za );

   /* The adjustment to AN. */
      adj = method ==1 ? asize : -asize;

   /* Apply the adjustment. */
      s = sin ( adj );
      c = cos ( adj );
      x = xa;
      z = za;
      xa = x * c - z * s;
      za = x * s + z * c;

   /* Return the adjusted [-HA,Dec] vector. */
      va2e ( sinphi, cosphi, isgim, ae2nm, xa, ya, za, vnew );

   /* Also apply the adjustment to the model */
   /* -------------------------------------- */

      if ( mount == EQUAT) {

      /* Equatorial case: apply corresponding adjustments to -IH and ID. */
         aihid ( vc, vnew, &MIH, &ID );

      } else {

      /* Altaz/gimbal case: apply the same adjustment directly to AN. */
         AN += asize;
      }

      break;

/*
** ------------------------------------
** East-west adjustment to azimuth axis
** ------------------------------------
**
** In the GIMBAL case, "azimuth axis" means the mount's longitude axis,
** and "east-west" means at right-angles to the mount's prime meridian.
*/

   case 'W':
   case 'w':

   /* Apply the adjustment to the given vector */
   /* ---------------------------------------- */

   /* Rotate the given [-HA,Dec] vector into [Az,El]. */
      ve2a ( sinphi, cosphi, isgim, ae2nm, xold, yold, zold, &xa, &ya, &za );

   /* The adjustment to AW. */
      adj = method == 1 ? asize : -asize;

   /* Apply the adjustment. */
      s = sin ( adj );
      c = cos ( adj );
      y = ya;
      z = za;
      ya = y * c - z * s;
      za = y * s + z * c;

   /* Return the adjusted [-HA,Dec] vector. */
      va2e ( sinphi, cosphi, isgim, ae2nm, xa, ya, za, vnew );

   /* Also apply the adjustment to the model */
   /* -------------------------------------- */

      if ( mount == EQUAT) {

      /* Equatorial case: apply corresponding adjustments to -IH and ID. */
         aihid ( vc, vnew, &MIH, &ID );

      } else {

      /* Altaz/gimbal case: apply the adjustment directly to AW. */
         AW += asize;
      }

      break;

/* Unrecognized adjustment type. */
   default:
      return -1;
   }

/* Return the adjusted [-HA,Dec] vector. */
   for ( i = 0; i < 3; i++ ) {
      vc[i] = vnew[i];
   }

/* Finished. */
   return 0;
}

static void aiaie ( double sinphi, double cosphi,
                    int isgim, double ae2nm[3][3],
                    double vold[3], double vnew[3],
                    double* pia, double* pie )
/*
**  - - - - - -
**   a i a i e
**  - - - - - -
**
**  Take unadjusted and adjusted [-HA,Dec] vectors and apply
**  corresponding adjustments to mount [Az,El].
**
**  Given:
**     sinphi   double       sine of telescope latitude (Note 1)
**     cosphi   double       cosine of telescope latitude (Note 1)
**     isgim    int          true if generalized gimbal case (Note 2)
**     ae2nm    double[3][3] [Az,El] to nominal mount matrix (Note 2)
**     vold     double[3]    unadjusted [-HA,Dec] vector
**     vnew     double[3]    adjusted [-HA,Dec] vector
**     pia      double*      unadjusted mount azimuth index error
**     pie      double*      unadjusted mount elevation index error
**
**  Returned:
**     pia      double*      adjusted mount azimuth index error
**     pie      double*      adjusted mount elevation index error
**
**  Notes:
**
**  1) The functions of latitude, sinphi and cosphi, are used as supplied,
**     without checks for range and consistency.
**
**  2) The isgim argument signifies the mount is a "generalized
**     gimbal", rotated with respect to [Az,El] in accordance with
**     the matrix ae2nm.
**
**  Called:   ve2a, xyz2ab
**
**  Last revision:   13 November 1999
**
**  Copyright 1999 P.T.Wallace.  All rights reserved.
*/

{
   double xa, ya, za, aold, bold, anew, bnew;



/* Old [-HA,Dec] vector to mount Az and El. */
   ve2a ( sinphi, cosphi, isgim, ae2nm, vold[0], vold[1], vold[2],
          &xa, &ya, &za );
   xyz2ab ( xa, ya, za, &aold, &bold );

/* New [-HA,Dec] to spherical mount [Az,El]. */
   ve2a ( sinphi, cosphi, isgim, ae2nm, vnew[0], vnew[1], vnew[2],
          &xa, &ya, &za );
   xyz2ab ( xa, ya, za, &anew, &bnew );

/* Adjust IA and IE according to how much Az and El have changed. */
   *pia += anew - aold;
   *pie += bnew - bold;
}

static void aihid ( double vold[3], double vnew[3],
                    double* pmih, double* pid )
/*
**  - - - - - -
**   a i h i d
**  - - - - - -
**
**  Take unadjusted and adjusted [-HA,Dec] vectors and apply
**  corresponding adjustments to [-HA,Dec].
**
**  Given:
**     vold     double[3]    unadjusted [-HA,Dec] vector
**     vnew     double[3]    adjusted [-HA,Dec] vector
**     pmih     double*      unadjusted HA index error
**     pid      double*      unadjusted Dec index error
**
**  Returned:
**     pmih     double*      adjusted minus HA index error
**     pid      double*      adjusted Dec index error
**
**  Called:   xyz2ab
**
**  Last revision:   26 November 1999
**
**  Copyright 1999 P.T.Wallace.  All rights reserved.
*/

{
   double aold, bold, anew, bnew;


/* Unadjusted and adjusted [-HA,Dec] to spherical. */
   xyz2ab ( vold[0], vold[1], vold[2], &aold, &bold );
   xyz2ab ( vnew[0], vnew[1], vnew[2], &anew, &bnew );

/* Adjust -IH and ID according to how much -HA and Dec have changed. */
   *pmih += anew - aold;
   *pid += bnew - bold;
}

static void ve2a ( double sinphi, double cosphi,
                   int isgim, double ae2nm[3][3],
                   double xe, double ye, double ze,
                   double* xa, double* ya, double* za )
/*
**  - - - - -
**   v e 2 a
**  - - - - -
**
**  [-HA,Dec] vector to [Az,El] or mount [roll,pitch].
**
**  Given:
**     sinphi   double       sine of telescope latitude (Note 1)
**     cosphi   double       cosine of telescope latitude (Note 1)
**     isgim    int          true if generalized gimbal case (Note 2)
**     ae2nm    double[3][3] [Az,El] to nominal mount matrix (Note 2)
**     xe       double       [-HA,Dec] x
**     ye       double       [-HA,Dec] y
**     ze       double       [-HA,Dec] z
**
**  Returned:
**     xa       double*      [Az,El] or [roll,pitch] x
**     ya       double*      [Az,El] or [roll,pitch] y
**     za       double*      [Az,El] or [roll,pitch] z
**
**  Notes:
**
**  1) The functions of latitude, sinphi and cosphi, are used as supplied,
**     without checks for range and consistency.
**
**  2) The isgim argument signifies the mount is a "generalized
**     gimbal", rotated with respect to [Az,El] in accordance with
**     the matrix ae2nm.
**
**  Last revision:   13 November 1999
**
**  Copyright 1999 P.T.Wallace.  All rights reserved.
*/

{
   double x, y, z, xm, ym, zm;



/* [-HA,Dec] to [Az,El]. */
   x = sinphi * xe - cosphi * ze;
   y = ye;
   z = cosphi * xe + sinphi * ze;

/* For gimbal, [Az,El] to [roll,pitch]. */
   if ( isgim ) {
      xm = ae2nm[0][0] * x + ae2nm[0][1] * y + ae2nm[0][2] * z;
      ym = ae2nm[1][0] * x + ae2nm[1][1] * y + ae2nm[1][2] * z;
      zm = ae2nm[2][0] * x + ae2nm[2][1] * y + ae2nm[2][2] * z;
      x = xm;
      y = ym;
      z = zm;
   }

/* Return the result. */
   *xa = x;
   *ya = y;
   *za = z;
}

static void va2e ( double sinphi, double cosphi,
                   int isgim, double ae2nm[3][3],
                   double xa, double ya, double za,
                   double ve[3] )
/*
**  - - - - -
**   v a 2 e
**  - - - - -
**
**  [Az,El] or mount [roll,pitch] vector to [-HA,Dec].
**
**  Given:
**     sinphi   double       sine of telescope latitude (Note 1)
**     cosphi   double       cosine of telescope latitude (Note 1)
**     isgim    int          true if generalized gimbal case (Note 2)
**     ae2nm    double[3][3] [Az,El] to nominal mount matrix (Note 2)
**     xa       double       [Az,El] or [roll,pitch] x
**     ya       double       [Az,El] or [roll,pitch] y
**     za       double       [Az,El] or [roll,pitch] z
**
**  Returned:
**     ve       double[3]    [-HA,Dec] vector
**
**  Notes:
**
**  1) The functions of latitude, sinphi and cosphi, are used as supplied,
**     without checks for range and consistency.
**
**  2) The isgim argument signifies the mount is a "generalized
**     gimbal", rotated with respect to [Az,El] in accordance with
**     the matrix ae2nm.
**
**  Last revision:   13 November 1999
**
**  Copyright 1999 P.T.Wallace.  All rights reserved.
*/

{
   double x, y, z;



/* If gimbal, [roll,pitch] to [Az,El]. */
   if ( isgim ) {
      x = ae2nm[0][0] * xa + ae2nm[1][0] * ya + ae2nm[2][0] * za;
      y = ae2nm[0][1] * xa + ae2nm[1][1] * ya + ae2nm[2][1] * za;
      z = ae2nm[0][2] * xa + ae2nm[1][2] * ya + ae2nm[2][2] * za;
   } else {
      x = xa;
      y = ya;
      z = za;
   }

/* [Az,El] to [-HA,Dec]. */
   ve[0] = sinphi * x + cosphi * z;
   ve[1] = y;
   ve[2] = - cosphi * x + sinphi * z;
}

static void xyz2ab ( double x, double y, double z, double* a, double* b )
/*
**  - - - - - - -
**   x y z 2 a b
**  - - - - - - -
**
**  Cartesian to spherical.
**
**  Given:
**     x,y,z    double       righthanded Cartesian coordinates
**
**  Returned:
**     a,b      double       righthanded longitude and latitude
**
**  Last revision:   13 November 1999
**
**  Copyright 1999 P.T.Wallace.  All rights reserved.
*/

{
   double r;


   r = sqrt ( x * x + y * y );
   *a = ( r != 0.0 ) ? atan2 ( y, x ) : 0.0;
   *b = ( z != 0.0 ) ? atan2 ( z, r ) : 0.0;
}

static void ab2xyz ( double a, double b, double* x, double* y, double* z )
/*
**  - - - - - - -
**   a b 2 x y z
**  - - - - - - -
**
**  Spherical to Cartesian.
**
**  Given:
**     a,b      double     righthanded longitude and latitude
**
**  Returned:
**     x,y,z    double*    righthanded Cartesian coordinates
**
**  Last revision:   13 November 1999
**
**  Copyright 1999 P.T.Wallace.  All rights reserved.
*/

{
   double cb;


   cb = cos ( b );
   *x = cos ( a ) * cb;
   *y = sin ( a ) * cb;
   *z = sin ( b );
}
