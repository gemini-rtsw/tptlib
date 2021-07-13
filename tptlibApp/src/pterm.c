int tptPterms ( int,
                double, double[3], double, double, int, double[3][3],
                int, int, double[3], double[3], double[8] );
int tptPterml ( int,
                double, double[3], double, double, int, double[3][3],
                int, int, double[3], double[3], double[8] );
int tptPtermg ( char[9], char[9],
                double, double[3], double, double, int, double[3][3],
                int, int, double[3], double[3], double[8] );

int tptPterm ( int nterm, double cval, double auxil[], int jbp, double vw[],
               double sinphi, double cosphi, int mount, double ae2nm[3][3],
               int nterml, int ntermx, int nterms,
               char coeffn[][9], char coform[][9], int method,
               double vc[3], double pmodel[7] )
/*
**  - - - - - - - - -
**   t p t P t e r m
**  - - - - - - - - -
**
**  Evaluate the effect of one term in the pointing model at one position.
**
**  Given (arguments):
**     nterm    int         coefficient number (1..., Note 1)
**     cval     double      coefficient value (radians)
**     auxil    double[]    auxiliary readings if any
**     jbp      int         "below pole" flag
**     vw       double[]    working [-h,d] vector for calculations
**     sinphi   double      sine of telescope latitude (Note 2)
**     cosphi   double      cosine of telescope latitude (Note 2)
**     mount    int         mount type (Note 3)
**     ae2nm    double[][]  [Az,El] to nominal mount matrix (Note 3)
**     nterml   int         number of terms implemented locally
**     ntermx   int         number of terms implemented explicitly (loc+std)
**     nterms   int         number of terms available currently (loc+std+gen)
**     coeffn   char[][]    coefficient names
**     coform   char[][]    tokenized form of generic term names
**     method   int         1 = add correction, else = subtract correction
**     vc       double[]    current [-h,d] vector (Note 4)
**     pmodel   double[]    basic pointing model (Note 5)
**
**  Returned (argument):
**     vc       double[]    updated
**     pmodel   double[]    updated
**
**  Returned (function value):
**              int        status: 0 = OK (Note 7)
**
**  Called:  tptPterml, tptPterms, tptPtermg
**
**  Notes:
**
**  1) It is up to the caller to make sure that nterm is within the
**     range 1 thru nterms.  If it is not, the Cartesian Az/El
**     supplied is not changed.
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
**  4) The "vc" [-h,d] vector is the direction to which the correction
**     is being applied.
**
**  5) The array "pmodel" is a list of seven coefficients which
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
**  6) The corrections are computed in the routines tptPterml (locally
**     implemented terms), tptPterms (standard terms) and tptPtermg
**     (generic terms).  These routines, like the present one, avoid use
**     of #includes, to facilitate their use in other contexts, for
**     example as part of a telescope control system.  The function's
**     usefulness in this regard is supplemented by the correction
**     being applied in two forms, as corrected coordinates and as
**     an updated "basic model".  The basic model is an 8-coefficient
**     generic model which summarizes, for the current telescope
**     position, the many and varied terms of the operational model.
**
**  7) The function value is the status returned by tptPterml,
**     tptPterms or tptPtermg, whichever is called.
**
**  Last revision:   15 September 2000
**
**  Copyright 2000 P.T.Wallace.  All rights reserved.
*/

{
/* Decide with which routine to apply the correction. */
   if ( nterm >= 1 && nterm <= nterms ) {

      if ( nterm <= nterml ) {

      /* Apply a LOCAL term. */
         return tptPterml ( nterm,
                            cval, auxil,
                            sinphi, cosphi, mount, ae2nm, jbp, method,
                            vw, vc, pmodel );

      } else if ( nterm <= ntermx ) {

      /* Apply a STANDARD term. */
         return tptPterms ( nterm,
                            cval, auxil,
                            sinphi, cosphi, mount, ae2nm, jbp, method,
                            vw, vc, pmodel );

      } else {

      /* Apply a GENERIC term. */
         return tptPtermg ( coeffn[nterm-1], coform[nterm-1],
                            cval, auxil,
                            sinphi, cosphi, mount, ae2nm, jbp, method,
                            vw, vc, pmodel );
      }
   }

/* No such term. */
   return -1;
}
