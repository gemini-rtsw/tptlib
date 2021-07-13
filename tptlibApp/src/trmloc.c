#include <string.h>
void tptTrmloc ( int *navail, char coeffn [ ] [ 9 ] )
/*
**  - - - - - - - - - -
**   t p t T r m l o c
**  - - - - - - - - - -
**
**  Insert into the coefficient name table the names of any LOCALLY
**  IMPLEMENTED pointing terms, and perform any other local
**  initialization actions.  The terms are defined mathematically
**  in the tptPterml routine.
**
**  The present routine should be edited to meet local requirements.
**
**  The STANDARD set of terms is added in the tptTrmstd routine.
**
**  Given:
**     navail   int*             number of vacant elements in coeffn array
**     coeffn   char[NTROOM][9]  term table (coefficient names)
**
**  Returned:
**     navail   int*             updated, less than zero if not enough room
**     coeffn   char[NTROOM][9]  updated
**
**  After calling this routine, navail can have the following values:
**     > 0 = number of vacant elements left in the coeffn array
**       0 = all terms were added, exactly filling the coeffn array
**     < 0 = coeffn array filled up before all terms could be added
**
**  Last revision:   21 August 1999
**
**  Copyright 1999 P.T.Wallace.  All rights reserved.
*/
{
#if 0

/* Note that all the code is disabled in this demo version. */

   int n, m;
 
/* Coefficient names */
   static char coname [ ] [ 9 ] = {
      "BEND    ",
      "SLIP    ",
      "JUDDER  ",
      "SHAKE   ",
      "RATTLE  ",
      "ROLL    ",
      "CRUNCH  ",
      "SLUMP   ",
      "SLOUCH  ",
      "WANDER  ",
      "WOBBLE  ",
      "DRIFT   "
   };
   int ncoeff = sizeof ( coname ) / sizeof ( coname [ 0 ] );
 
 
/* Add coefficient names to the table. */
   m = ( ncoeff < *navail ) ? ncoeff : *navail;
   for ( n = 0; n < m; n++ ) {
      strcpy ( coeffn [ n ], coname [ n ] );
   }
 
/* Set navail to indicate spare/full/overflow. */
   *navail -= ncoeff;

#endif
}
