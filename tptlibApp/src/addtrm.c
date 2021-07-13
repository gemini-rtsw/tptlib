#include "tptlib.h"

int tcsAddtrm ( int maxtrm, int ntroom,
                char cname[9], double cvalue, int* nterms,
                char coeffn[][9], char coform[][9],
                int model[], double coeffv[] )

/*
**  - - - - - - - - - -
**   t c s A d d t r m
**  - - - - - - - - - -
**
**  Add a pointing term to the model if it is not already present,
**  and set its value.
**
**  Given:
**     maxtrm  int             maximum number of terms in model
**     ntroom  int             maximum number of terms in repertoire
**     cname   char[9]         name of term
**     cvalue  double          coefficient value (radians)
**     nterms  int*            # terms available currently (loc+std+gen)
**     coeffn  char[ntroom][9] names of coefficients
**
**  Returned:
**     nterms  int*            updated
**     coeffn  char[ntroom][9] updated
**     coform  char[ntroom][9] format of generic terms added to coeffn
**     model   double[maxtrm]  term numbers for model (0 = end)
**     coeffv  double[maxtrm]  coefficient values
**
**  Returned (function value):
**             int      -2 = no room left in model
**                      -1 = unrecognized term
**                       0 = no room in repertoire for new generic term
**                      >0 = term number (1...)
**
**  Called:  tptTerm
**
**  Last revision:   19 July 2000
**
**  Copyright 2000 P.T.Wallace.  All rights reserved.
*/

{
   int n, i;


/* Make sure the term is in the repertoire and get its number. */
   if ( ( n = tptTerm ( cname, ntroom, nterms, coeffn, coform ) ) <= 0 )
      return n;

/* It is (or is now);  locate the existing, or first unused, model slot. */
   for ( i = 0; i < maxtrm; i++ ) {

   /* If unused slot, add the new term to the model. */
      if ( ! model[i] ) model[i] = n;

   /* If this is the right term, set the coefficient value and exit. */
      if ( n == model[i] ) {
         coeffv[i] = cvalue;
         return n;
      }
   }

/* The model is full. */
   return -2;

}
