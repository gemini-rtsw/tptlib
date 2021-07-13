void tptTrmloc ( int*, char[][9] );
void tptTrmstd ( int*, char[][9] );

int tptMinit ( int maxtrm, int ntroom, int model[], int *nterml,
               int *ntermx, int *nterms, char coeffn[][9] )

/*
**  - - - - - - - - -
**   t p t M i n i t
**  - - - - - - - - -
**
**  Initialize the pointing model.
**
**  Given (arguments):
**     maxtrm     int        maximum number of terms in model
**     ntroom     int        maximum number of terms allowed for
**
**  Returned (arguments):
**     model      int[]      term numbers for current model (0 = end)
**     nterml     int        number of local terms
**     ntermx     int        number of terms implemented explicitly (loc+std)
**     nterms     int        number of terms available currently (loc+std+gen)
**     coeffn     char[][]   coefficient names (local, standard, generic)
**
**  Returned (function value):
**                int        status: 0 = OK
**                                  -1 = insufficient room for all terms
**
**  Called:  tptTrmloc, tptTrmstd
**
**  Last revision:   22 August 1999
**
**  Copyright 1999 P.T.Wallace.  All rights reserved.
*/

{

   int n, navail;



/* Fill the model with endmarks. */
   for ( n = 0; n < maxtrm; n++ ) {
      model[n] = 0;
   }

/* Coefficient name table: preset to empty. */
   for ( n = 0; n < ntroom; n++ ) {
      coeffn[n][0] = (char) '\0';
   }
   navail = ntroom;

/* LOCAL terms. */
   tptTrmloc ( &navail, coeffn );
   if ( navail < 0 ) return -1;

/* Number of local terms. */
   *nterml = ntroom - navail;

/* STANDARD terms. */
   tptTrmstd ( &navail, coeffn + *nterml );
   if ( navail < 0 ) return -1;

/* Number of non-generic terms. */
   *ntermx = ntroom - navail;

/* Total number of terms (n.b. no GENERIC terms yet - dynamically loaded). */
   *nterms = *ntermx;

   return 0;
}
