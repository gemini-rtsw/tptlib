#include <string.h>
void tptTrmstd ( int *navail, char coeffn [ ] [ 9 ] )
/*
**  - - - - - - - - - -
**   t p t T r m s t d
**  - - - - - - - - - -
**
**  Insert into the coefficient name table the names of the STANDARD
**  pointing terms.  The terms are defined mathematically in the
**  tptPterms routine.
**
**  Any LOCALLY IMPLEMENTED terms are added in the tptTrmloc routine.
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
**  Last revision:   18 May 2000
**
**  Copyright 2000 P.T.Wallace.  All rights reserved.
*/
{
   int n, m;
 
/* Coefficient names */
   static char coname [ ] [ 9 ] = {
      "IH",
      "ID",
      "ZH",
      "NP",
      "CH",
      "HF",
      "ZE",
      "ME",
      "MA",
      "TF",
      "TX",
      "FO",
      "IE",
      "IA",
      "HCES",
      "HCEC",
      "DCES",
      "DCEC",
      "DNP",
      "X2HC",
      "HGES",
      "HGEC",
      "DGES",
      "DGEC",
      "TFP",
      "CA",
      "CD4A",
      "CD4B",
      "CD5A",
      "CD5B",
      "AN",
      "AW",
      "NPAE",
      "ACES",
      "ACEC",
      "ECES",
      "ECEC",
      "NRX",
      "NRY",
      "AUX1H",
      "AUX1X",
      "AUX1D",
      "AUX1A",
      "AUX1S",
      "AUX1E",
      "AUX2H",
      "AUX2X",
      "AUX2D",
      "AUX2A",
      "AUX2S",
      "AUX2E",
      "DAF",
      "HFX",
      "HFD",
      "FLOP",
      "CRX",
      "CRY"
   };
   int ncoeff = sizeof ( coname ) / sizeof ( coname [ 0 ] );
 

/* Add coefficient names to the table. */
   m = ( ncoeff < *navail ) ? ncoeff : *navail;
   for ( n = 0; n < m; n++ ) {
      strcpy ( coeffn [ n ], coname [ n ] );
   }
 
/* Set navail to indicate spare/full/overflow. */
   *navail -= ncoeff;

}
