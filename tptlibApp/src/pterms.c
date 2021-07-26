#include <math.h>

void tptPtpre ( double, double, int, double[3][3], int, double[3],
                double*, double*, double*, double*, double*, double*,
                double*, double*, double*, double*, double* );
int tptPtappl ( char, double, double, double, int, double[3][3], int, int,
                double[3], double[8] );

int tptPterms ( int nterm, double cval, double auxil[3],
                double sinphi, double cosphi, int mount,
                double ae2nm[3][3], int jbp, int method,
                double vw[3], double vc[3], double pmodel[7] )
/*
**  - - - - - - - - - -
**   t p t P t e r m s
**  - - - - - - - - - -
**
**  Compute the effect of one STANDARD term in the pointing model.
**
**  The effects of LOCALLY IMPLEMENTED and GENERIC terms are computed
**  in the routines tptPterml and tptPtermg.
**
**  n.b.  Changes in this routine must be coordinated with changes
**        in the routine tptTrmstd.
**
**  Given (arguments):
**     nterm    int          coefficient number within standard set (Note 1)
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
**     within the allowed range, i.e. 1 to the number of standard terms
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
**     one standard pointing term at the current position.
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
**     The routine also serves as a template for the tptPterml routine,
**     used for implementing locally-defined pointing terms.
**
**  Last revision:   15 September 2000
**
**  Copyright 2000 P.T.Wallace.  All rights reserved.
*/

#define TINY 1e-20

{
   double xa, ya, za, az, el, zd, xe, ye, ze, ha, dec, ef1, ef2 = 0.0, w,
          h2, h3, h5, d2, sh, s2h, s3h, h12, d3, d6, ch;
   char type1, type2;
   int j;



/*
** ---------------------------------------------
** Transform working Az/El vector to other forms
** ---------------------------------------------
*/

   tptPtpre ( sinphi, cosphi, mount, ae2nm, jbp, vw,
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
** IH: index error in HA
**
** Change in HA equal to coefficient IH.
**
** In equatorial spherical coordinates:
**    dHA  = +IH
**    dDec = 0
*/
   case 1:
      ef1 = 1.0;
      type1 = 'H';
      break;

/*--------------------------------------------------------------------
**
** ID: index error in Dec
**
** Change in Dec equal to coefficient ID.
**
** In equatorial spherical coordinates:
**    dHA  = 0
**    dDec = +ID
*/
   case 2:
      ef1 = 1.0;
      type1 = 'D';
      break;

/*--------------------------------------------------------------------
**
** ZH: AAT HA Z-gear effect in HA
**
** Change in HA equal to coefficient ZH multiplied by first
** auxiliary reading.
**
** In equatorial spherical coordinates:
**    dHA  = +ZH*AUX1
**    dDec = 0
*/
   case 3:
      ef1 = auxil[0];
      type1 = 'H';
      break;

/*--------------------------------------------------------------------
**
** NP: HA/Dec nonperpendicularity
**
** Change EW on the sky equal to coefficient NP multiplied by sin(Dec).
**
** For north declinations, +ve NP displaces the pointing axis to
** the west of where it would otherwise be.
**
** In equatorial spherical coordinates (approximately):
**    dHA  = +NP*tan(DEC)
**    dDec = 0
*/
   case 4:
      ef1 = 1.0;
      type1 = 'P';
      break;

/*--------------------------------------------------------------------
**
** CH: Dec/optical axis nonperpendicularity
**
** Change EW on the sky equal to coefficient CH.
**
** +ve CH displaces the pointing axis to the west of where it would
** otherwise be.
**
** In equatorial spherical coordinates (approximately):
**    dHA  = +CH/cos(DEC)
**    dDec = 0
*/
   case 5:
      ef1 = 1.0;
      type1 = 'X';
      break;

/*--------------------------------------------------------------------
**
** HF: AAT horseshoe flexure
**
** Change EW on the sky equal to coefficient HF multiplied by sin(HA).
**
** In equatorial spherical coordinates (approximately):
**    dHA  = +HF*sin(HA)/cos(DEC)
**    dDec = 0
*/
   case 6:
      ef1 = sin ( ha );
      type1 = 'X';
      break;

/*--------------------------------------------------------------------
**
** ZE: AAT HA Z-gear effect in polar axis elevation
**
** Rotation about EW axis equal to coefficient ZE multiplied by
** first auxiliary reading.
**
** In equatorial spherical coordinates (approximately):
**    dHA  = +ZE*AUX1*sin(HA)*tan(DEC)
**    dDec = +ZE*AUX1*cos(HA)
*/
   case 7:
      ef1 = auxil[0];
      type1 = 'U';
      break;

/*--------------------------------------------------------------------
**
** ME: polar axis elevation error
**
** Rotation about EW axis equal to coefficient ME.
**
** In the northern hemisphere, +ve ME means that the pole of the
** mounting is below the true (unrefracted) pole.  A mounting aligned
** the refracted pole will have -ve ME.
**
** In the southern hemisphere, +ve ME means that the pole of the
** mounting is above the true (unrefracted) pole.  A mounting aligned
** the refracted pole will have +ve ME.
**
** In equatorial spherical coordinates (approximately):
**    dHA  = +ME*sin(HA)*tan(DEC)
**    dDec = +ME*cos(HA)
*/
   case 8:
      ef1 = 1.0;
      type1 = 'U';
      break;

/*--------------------------------------------------------------------
**
** MA: polar axis error LR
**
** Rotation about axis through (HA=0,Dec=0) equal to coefficient MA.
**
** In the northern hemisphere, +ve MA means that the pole of the
** mounting is to the right of due north.
**
** In the southern hemisphere, +ve MA means that the pole of the
** mounting is to the right of due south.
**
** In equatorial spherical coordinates (approximately):
**    dHA  = -MA*cos(HA)*tan(DEC)
**    dDec = +MA*sin(HA)
*/
   case 9:
      ef1 = 1.0;
      type1 = 'L';
      break;

/*--------------------------------------------------------------------
**
** TF: tube flexure - sin Z law
**
** Change in ZD equal to coefficient TF multiplied by sin(ZD).
**
** In altazimuth spherical coordinates:
**    dAz = 0
**    dEl = -TF*cos(EL)
**
** In equatorial spherical coordinates (approximately):
**    dHA  = +TF*cos(PHI)*sin(HA)/cos(DEC)
**    dDec = +TF*(cos(PHI)*cos(HA)*sin(DEC)-sin(PHI)*cos(DEC))
*/
   case 10:
      ef1 = sin ( zd );
      type1 = 'Z';
      break;

/*--------------------------------------------------------------------
**
** TX: tube flexure - tan Z law
**
** Change in ZD equal to coefficient TX multiplied by tan(ZD).
**
** In altazimuth spherical coordinates:
**
**    dAz = 0
**    dEl = -TX*cot(EL)
**
** In equatorial spherical coordinates (approximately):
**
**    dHA  = +TX*cos(PHI)*sin(HA)/
**           ((sin(DEC)*sin(PHI)+cos(DEC)*cos(HA)*cos(PHI))*cos(DEC))
**
**    dDec = +TX*(cos(PHI)*cos(HA)*sin(DEC)-sin(PHI)*cos(DEC))/
**           (sin(DEC)*sin(PHI)+cos(DEC)*cos(HA)*cos(PHI))
*/
   case 11:
      w = cos ( zd );
      ef1 = sin ( zd ) / ( w != 0.0 ? w : TINY );
      type1 = 'Z';
      break;

/*--------------------------------------------------------------------
**
** FO: fork flexure
**
** Change in Dec equal to coefficient FO multiplied by cos(HA).
**
** In equatorial spherical coordinates:
**    dHA  = 0
**    dDec = +FO*cos(HA)
*/
   case 12:
      ef1 = cos ( ha );
      type1 = 'D';
      break;

/*--------------------------------------------------------------------
**
** IE: index error in elevation
**
** Change in El equal to coefficient IE.
**
** In altazimuth spherical coordinates:
**    dAz = 0
**    dEl = +IE
*/
   case 13:
      ef1 = 1.0;
      type1 = 'E';
      break;

/*--------------------------------------------------------------------
**
** IA: index error in azimuth
**
** Change in Az equal to coefficient IA.
**
** In Altazimuth spherical coordinates (n.b. S=0, E=90 !):
**    dAz = +IA
**    dEl = 0
*/
   case 14:
      ef1 = 1.0;
      type1 = 'A';
      break;

/*--------------------------------------------------------------------
**
** HCES: HA centring error (sin component)
**
** Change in HA equal to coefficient HCES multiplied by sin(HA).
**
** In equatorial spherical coordinates:
**    dHA  = +HCES*sin(HA)
**    dDec = 0
*/
   case 15:
      ef1 = sin ( ha );
      type1 = 'H';
      break;

/*--------------------------------------------------------------------
**
** HCEC: HA centring error (cos component)
**
** Change in HA equal to coefficient HCEC multiplied by cos(HA).
**
** In equatorial spherical coordinates:
**    dHA  = +HCEC*cos(HA)
**    dDec = 0
*/
   case 16:
      ef1 = cos ( ha );
      type1 = 'H';
      break;

/*--------------------------------------------------------------------
**
** DCES: Dec centring error (sin component)
**
** Change in Dec equal to coefficient DCES multiplied by sin(Dec).
**
** In equatorial spherical coordinates:
**    dHA  = 0
**    dDec = +DCES*sin(DEC)
*/
   case 17:
      ef1 = ze;
      type1 = 'D';
      break;

/*--------------------------------------------------------------------
**
** DCEC: Dec centring error (cos component)
**
** Change in Dec equal to coefficient DCEC multiplied by cos(Dec).
**
** In equatorial spherical coordinates:
**    dHA  = 0
**    dDec = +DCEC*cos(DEC)
*/
   case 18:
      ef1 = cos ( dec );
      type1 = 'D';
      break;

/*--------------------------------------------------------------------
**
** DNP: Dynamic nonperpendicularity
**
** Change EW on the sky equal to coefficient DNP multiplied
** by sin(HA)*sin(Dec).
**
** In equatorial spherical coordinates (approximately):
**    dHA  = +DNP*sin(HA)*tan(DEC)
**    dDec = 0
*/
   case 19:
      ef1 = sin ( ha );
      type1 = 'P';
      break;

/*--------------------------------------------------------------------
**
** X2HC: cos (2H) term EW
**
** Change EW on the sky equal to coefficient X2HC multiplied by cos(2*HA).
**
** In equatorial spherical coordinates (approximately):
**    dHA  = +X2HC*cos(2*HA)/cos(DEC)
**    dDec = 0
*/
   case 20:
      ef1 = cos ( ha + ha );
      type1 = 'X';
      break;

/*--------------------------------------------------------------------
**
** HGES: 36-min gear error in HA - sin
**
** Change in HA equal to coefficient HGES multiplied by sin(40*HA).
**
** In equatorial spherical coordinates:
**    dHA  = +HGES*sin(40*HA)
**    dDec = 0
*/
   case 21:
      ef1 = sin ( 40.0 * ha );
      type1 = 'H';
      break;

/*--------------------------------------------------------------------
**
** HGEC: 36-min gear error in HA - cos
**
** Change in HA equal to coefficient HGEC multiplied by cos(40*HA).
**
** In equatorial spherical coordinates:
**    dHA  = +HGEC*cos(40*HA)
**    dDec = 0
*/
   case 22:
      ef1 = cos ( 40.0 * ha );
      type1 = 'H';
      break;

/*--------------------------------------------------------------------
**
** DGES: 9-deg gear error in Dec - sin
**
** Change in Dec equal to coefficient DGES multiplied by sin(40*Dec).
**
** In equatorial spherical coordinates:
**    dHA  = 0
**    dDec = +DGES*sin(40*Dec)
*/
   case 23:
      ef1 = sin ( 40.0 * dec );
      type1 = 'D';
      break;

/*--------------------------------------------------------------------
**
** DGEC: 9-deg gear error in Dec - cos
**
** Change in Dec equal to coefficient DGEC multiplied by cos(40*Dec).
**
** In equatorial spherical coordinates:
**    dHA  = 0
**    dDec = +DGEC*cos (40*Dec)
*/
   case 24:
      ef1 = cos ( 40.0 * dec );
      type1 = 'D';
      break;

/*--------------------------------------------------------------------
**
** TFP: AAT tube flexure - non-Hooke's-law term
**
** Change in ZD equal to coefficient TFP multiplied by a
** canonical polynomial in ZD.  The polynomial currently
** favoured is simply a Z**4 term.  The sign is chosen to
** match past practice;  +ve TFP values mean the tube
** droops.
*/
   case 25:
      w = zd * zd;
      ef1 = w * w;
      type1 = 'Z';
      break;

/*--------------------------------------------------------------------
**
** CA: collimation error in azimuth
**
** Change LR on the sky equal to coefficient CA.
**
** In altazimuth spherical coordinates (n.b. S=0, E=90 !):
**    dAz = +CA/cos(El)
**    dEl = 0
** (for small CA and not close to the zenith).
*/
   case 26:
      ef1 = 1.0;
      type1 = 'S';
      break;

/*--------------------------------------------------------------------
**
** CD4A: AAT coude 4 collimation error A component
**
** Change in HA equal to coefficient CD4A, and change in Dec
** equal to the same coefficient multiplied by sin(Dec).
**
** In equatorial spherical coordinates:
**    dHA  = +CD4A
**    dDec = +CD4A*sin(DEC)
*/
   case 27:
      ef1 = 1.0;
      type1 = 'H';
      ef2 = ze;
      type2 = 'D';
      break;

/*--------------------------------------------------------------------
**
** CD4B: AAT coude 4 collimation error B component
**
** Change EW on the sky equal to coefficient CD4B multiplied by
** -sin(Dec), and change in Dec equal to the same coefficient
** multiplied by cos(Dec).
**
** In equatorial spherical coordinates (approximately):
**    dHA  = -CD4B*tan(DEC)
**    dDec = +CD4B*cos(DEC)
*/
   case 28:
      ef1 = - ze;
      type1 = 'X';
      ef2 = cos ( dec );
      type2 = 'D';
      break;

/*--------------------------------------------------------------------
**
** CD5A: AAT coude 5 collimation error A component
**
** Change EW on the sky equal to coefficient CD5A multiplied by
** -sin(HA+Dec), and change in Dec equal to the same coefficient
** multiplied by cos(HA+Dec).
**
** In equatorial spherical coordinates (approximately):
**    dHA  = -CD5A*sin(HA+Dec)/cos(Dec)
**    dDec = +CD5A*cos(HA+Dec)
*/
   case 29:
      w = ha + dec;
      ef1 = - sin (w);
      type1 = 'X';
      ef2 = cos ( w );
      type2 = 'D';
      break;

/*--------------------------------------------------------------------
**
** CD5B: AAT coude 5 collimation error B component
**
** Change EW on the sky equal to coefficient CD5B multiplied by
** cos(HA+Dec), and change in Dec equal to the same coefficient
** multiplied by sin(HA+Dec).
**
** In equatorial spherical coordinates (approximately):
**    dHA  = -CD5B*cos(HA+Dec)/cos(Dec)
**    dDec = +CD5B*sin(HA+Dec)
*/
   case 30:
      w = ha + dec;
      ef1 = cos (w);
      type1 = 'X';
      ef2 = sin ( w );
      type2 = 'D';
      break;

/*--------------------------------------------------------------------
**
** AN: NS misalignment of azimuth axis
**
** Rotation about a horizontal EW axis equal to coefficient AN: the
** pole of the mounting is north of the zenith.
**
** In altazimuth spherical coordinates (n.b. S=0, E=90 !):
**    dAz = +AN*sin(Az)*tan(El)
**    dEl = +AN*cos(Az)
*/
   case 31:
      ef1 = 1.0;
      type1 = 'N';
      break;

/*--------------------------------------------------------------------
**
** AW: EW misalignment of azimuth axis
**
** Rotation about a horizontal NS axis equal to coefficient AW: the
** pole of the mounting is west of the zenith.
**
** In altazimuth spherical coordinates (n.b. S=0, E=90 !):
**    dAz = -AW*cos(Az)*tan(El)
**    dEl = +AW*sin(Az)
*/
   case 32:
      ef1 = 1.0;
      type1 = 'W';
      break;

/*--------------------------------------------------------------------
**
** NPAE: nonperpendicularity of azimuth and elevation axes
**
** Change LR on the sky equal to coefficient NPAE
** multiplied by sin(El).
**
** In altazimuth spherical coordinates (n.b. S=0, E=90 !):
**    dAz = +NPAE*tan(El)
**    dEl = 0
*/
   case 33:
      ef1 = 1.0;
      type1 = 'V';
      break;

/*--------------------------------------------------------------------
**
** ACES: azimuth centring error (sin component)
**
** Change in Az equal to coefficient ACES multiplied by sin(Az).
**
** In altazimuth spherical coordinates (n.b. S=0, E=90 !):
**    dAz = +ACES*sin(AZ)
**    dEl = 0
*/
   case 34:
      ef1 = sin ( az );
      type1 = 'A';
      break;

/*--------------------------------------------------------------------
**
** ACEC: azimuth centring error (cos component)
**
** Change in Az equal to coefficient ACEC multiplied by cos(Az).
**
** In altazimuth spherical coordinates (n.b. S=0, E=90 !):
**    dAz = +ACEC*cos(Az)
**    dEl = 0
*/
   case 35:
      ef1 = cos ( az );
      type1 = 'A';
      break;

/*--------------------------------------------------------------------
**
** ECES: elevation centring error (sin component)
**
** Change in El equal to coefficient ECES multiplied by sin(El).
**
** In altazimuth spherical coordinates:
**    dAz = 0
**    dEl = +ECES*sin(EL)
*/
   case 36:
      ef1 = za;
      type1 = 'E';
      break;

/*--------------------------------------------------------------------
**
** ECEC: elevation centring error (cos component)
**
** Change in El equal to coefficient ECEC multiplied by cos (El).
**
** In altazimuth spherical coordinates:
**    dAz = 0
**    dEl = +ECEC*cos(El)
*/
   case 37:
      ef1 = cos ( el );
      type1 = 'E';
      break;

/*--------------------------------------------------------------------
**
** NRX: Nasmyth rotator displacement (vertical)
**
** Change LR on the sky equal to coefficient NRX multiplied
** by cos(El), and change in El equal to the same coefficient
** multiplied by -sin(El).
**
** In altazimuth spherical coordinates (n.b. S=0, E=90 !):
**    dAz = +NRX
**    dEl = -NRX*sin(El)
*/
   case 38:
      ef1 = cos ( el );
      type1 = 'S';
      ef2 = - za;
      type2 = 'E';
      break;

/*--------------------------------------------------------------------
**
** NRY: Nasmyth rotator displacement (horizontal)
**
** Change LR on the sky equal to coefficient NRY multiplied
** by sin(El), and change in El equal to the same coefficient
** multiplied by cos(El).
**
** In altazimuth spherical coordinates (n.b. S=0, E=90 !):
**    dAz = +NRY*tan(El)
**    dEl = +NRY*cos(El)
*/
   case 39:
      ef1 = za;
      type1 = 'S';
      ef2 = cos ( el );
      type2 = 'E';
      break;

/*--------------------------------------------------------------------
**
** AUX1H: HA change given by auxiliary reading
**
** Change in HA equal to coefficient AUX1H multiplied by first
** auxiliary reading.
**
** In equatorial spherical coordinates:
**    dHA  = +AUX1H*AUX1
**    dDec = 0
*/
   case 40:
      ef1 = auxil[0];
      type1 = 'H';
      break;

/*--------------------------------------------------------------------
**
** AUX1X: EW change given by auxiliary reading
**
** Change EW on the sky equal to coefficient AUX1X multiplied by
** first auxiliary reading.
**
** In equatorial spherical coordinates (approximately):
**    dHA  = +AUX1X*AUX1/cos(Dec)
**    dDec = 0
*/
   case 41:
      ef1 = auxil[0];
      type1 = 'X';
      break;

/*--------------------------------------------------------------------
**
** AUX1D: Dec change given by auxiliary reading
**
** Change in Dec equal to coefficient AUX1D multiplied by first
** auxiliary reading.
**
** In equatorial spherical coordinates:
**    dHA  = 0
**    dDec = AUX1D*AUX1
*/
   case 42:
      ef1 = auxil[0];
      type1 = 'D';
      break;

/*--------------------------------------------------------------------
**
** AUX1A: Az change given by auxiliary reading
**
** Change in Az equal to coefficient AUX1A multiplied by first
** auxiliary reading.
**
** In altazimuth spherical coordinates (n.b. S=0, E=90 !):
**    dAz = +AUX1A*AUX1
**    dEl = 0
*/
   case 43:
      ef1 = auxil[0];
      type1 = 'A';
      break;

/*--------------------------------------------------------------------
**
** AUX1S: LR change given by auxiliary reading
**
** Change LR on the sky equal to coefficient AUX1S multiplied by
** first auxiliary reading.
**
** In altazimuth spherical coordinates (n.b. S=0, E=90 !):
**    dAz = +AUX1S*AUX1/cos(EL)
**    dEl = 0
*/
   case 44:
      ef1 = auxil[0];
      type1 = 'S';
      break;

/*--------------------------------------------------------------------
**
** AUX1E: El change given by auxiliary reading
**
** Change in El equal to coefficient AUX1E multiplied by first
** auxiliary reading.
**
** In altazimuth spherical coordinates:
**    dAz = 0
**    dEl = +AUX1E*AUX1
*/
   case 45:
      ef1 = auxil[0];
      type1 = 'E';
      break;

/*--------------------------------------------------------------------
**
** AUX2H: HA change given by auxiliary reading
**
** Change in HA equal to coefficient AUX1H multiplied by second
** auxiliary reading.
**
** In equatorial spherical coordinates:
**    dHA  = +AUX2H*AUX2
**    dDec = 0
*/
   case 46:
      ef1 = auxil[1];
      type1 = 'H';
      break;

/*--------------------------------------------------------------------
**
** AUX2X: EW change given by auxiliary reading
**
** Change EW on the sky equal to coefficient AUX2X multiplied by
** second auxiliary reading.
**
** In equatorial spherical coordinates (approximately):
**    dHA  = +AUX2X*AUX2/cos(Dec)
**    dDec = 0
*/
   case 47:
      ef1 = auxil[1];
      type1 = 'X';
      break;

/*--------------------------------------------------------------------
**
** AUX2D: Dec change given by auxiliary reading
**
** Change in Dec equal to coefficient AUX2D multiplied by second
** auxiliary reading.
**
** In equatorial spherical coordinates:
**    dHA  = 0
**    dDec = AUX2D*AUX2
*/
   case 48:
      ef1 = auxil[1];
      type1 = 'D';
      break;

/*--------------------------------------------------------------------
**
** AUX2A: Az change given by auxiliary reading
**
** Change in Az equal to coefficient AUX2A multiplied by second
** auxiliary reading.
**
** In altazimuth spherical coordinates (n.b. S=0, E=90 !):
**    dAz = +AUX2A*AUX2
**    dEl = 0
*/
   case 49:
      ef1 = auxil[1];
      type1 = 'A';
      break;

/*--------------------------------------------------------------------
**
** AUX2S: LR change given by auxiliary reading
**
** Change LR on the sky equal to coefficient AUX2S multiplied by
** second auxiliary reading.
**
** In altazimuth spherical coordinates (n.b. S=0, E=90 !):
**    dAz = +AUX2S*AUX2/cos(EL)
**    dEl = 0
*/
   case 50:
      ef1 = auxil[1];
      type1 = 'S';
      break;

/*--------------------------------------------------------------------
**
** AUX2E: El change given by auxiliary reading
**
** Change in El equal to coefficient AUX2E multiplied by second
** auxiliary reading.
**
** In altazimuth spherical coordinates:
**    dAz = 0
**    dEl = +AUX2E*AUX2
*/
   case 51:
      ef1 = auxil[1];
      type1 = 'E';
      break;

/*--------------------------------------------------------------------
**
** DAF: Declination axis flexure
**
** Downward sag of cantilevered declination axis (e.g. English
** cross-axis and German equatorial mountings) proportional to
** the sine of the zenith distance of the declination axis, the
** latter emerging from the polar axis towards the west when the
** telescope points at the meridian).
**
** In equatorial spherical coordinates (approximately):
**    dHA  = -DAF*(cos(PHI)*cos(HA)+sin(PHI)*tan(DEC))
**    dDec = 0
*/
   case 52:
      ef1 = - ( cosphi * cos ( ha ) * cos ( dec ) + sinphi * ze );
      type1 = 'X';
      break;

/*--------------------------------------------------------------------
**
** HFX: AAT horseshoe flexure EW
**
** Change EW equal to coefficient HFX multiplied by a canonical
** expression in HA and Dec.
*/
   case 53:
      h2 = ha + ha;
      h3 = h2 + ha;
      h5 = h3 + h2;
      d2 = dec + dec;
      sh = sin ( ha );
      s2h = sin ( h2 );
      s3h = sin ( h3 );
      ef1 =   0.219 * sin ( h5 )
           + 0.368 * cos ( h5 + ha )
           + 0.195 * sin ( h5 + h5 + h5 )
       + (   1.279 * s2h
           - 1.718 * sin ( h2 + h2 )
           + 1.202 * cos ( h5 ) ) * ze
           - 1.928 * s3h * cos ( dec )
       + (   3.662 * sh
           + 0.381 * cos ( h3 ) ) * sin ( d2 )
           - 1.156 * s3h * cos ( d2 )
       + ( - 0.633 * s2h
           + 0.723 * sh ) * sin ( d2 + d2 + dec );
      type1 = 'X';
      break;

/*--------------------------------------------------------------------
**
** HFD: AAT horseshoe flexure NS
**
** Change NS equal to coefficient HFD multiplied by a canonical
** expression in HA and Dec.
*/
   case 54:
      h2 = ha + ha;
      h3 = h2 + ha;
      h5 = h3 + h2;
      h12 = h5 + h5 + h2;
      d2 = dec + dec;
      d3 = d2 + dec;
      d6 = d3 + d3;
      ch = cos ( ha );
      ef1 = - 2.529 * sin ( h2 )
           + 3.964 * sin ( h3 )
           + 3.008 * cos ( h3 )
           - 2.428 * sin ( h3 + ha )
           - 0.921 * cos ( h5 )
           - 0.233 * cos ( h12 )
           + 0.163 * sin ( d6 )
           - 0.147 * cos ( d6 + d3 + d2 )
           - 0.244 * cos ( h12 + h3 ) * ze
       + ( - 4.330 * ch
           - 2.958 * cos ( h2 )
           - 1.805 * sin ( ha ) ) * cos ( dec )
           + 1.666 * ch * cos ( d2 );
      type1 = 'D';
      break;

/*--------------------------------------------------------------------
**
** FLOP: constant error in zenith distance
**
** Change in ZD equal to coefficient FLOP.
**
** In altazimuth spherical coordinates:
**    dAz = 0
**    dEl = -FLOP
*/
   case 55:
      ef1 = 1.0;
      type1 = 'Z';
      break;

/*--------------------------------------------------------------------
**
** CRX: Altazimuth coude displacement (NS)
**
** Change LR on the sky equal to coefficient CRX multiplied
** by cos(Az+El), and change in El equal to the same coefficient
** multiplied by -sin(Az+El).
**
** In altazimuth spherical coordinates (n.b. S=0, E=90 !):
**    dAz = +CRX*cos(AZ+EL)/cos(EL)
**    dEl = -CRX*sin(AZ+EL)
*/
   case 56:
      w = az + el;
      ef1 = cos ( w );
      type1 = 'S';
      ef2 = - sin ( w );
      type2 = 'E';
      break;

/*--------------------------------------------------------------------
**
** CRY: Altazimuth coude displacement (EW)
**
** Change LR on the sky equal to coefficient CRY multiplied
** by sin(Az+El), and change in El equal to the same coefficient
** multiplied by cos(Az+El).
**
** In altazimuth spherical coordinates (n.b. S=0, E=90 !):
**    dAz = CRY*sin(AZ+EL)/COS(EL)
**    dEl = CRY*cos(AZ+EL)
*/
   case 57:
      w = az + el;
      ef1 = sin ( w );
      type1 = 'S';
      ef2 = cos ( w );
      type2 = 'E';
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

/* Success. */
   return 0;
}
