#include <stdio.h>
#include <string.h>

int tptScomp ( char*, char* );
char* tptTrmgen ( char* );

int tptTerm ( char *name, int ntroom,
              int *nterms, char coeffn[][9], char coform[][9] )
/*
**  - - - - - - - -
**   t p t T e r m
**  - - - - - - - -
**
**  Translate a term name into the corresponding term number.  If
**  the term is a generic one that has not previously been referred
**  to, it is added to the current repertoire ready for if and when
**  it is needed.
**
**  Given (arguments):
**     name    char[]    name of term
**     ntroom  int       maximum number of terms in repertoire
**     nterms  int*      number of terms available currently (loc+std+gen)
**     coeffn  char[][]  names of coefficients
**
**  Returned (arguments):
**     nterms  int*      updated
**     coeffn  char[][]  updated
**     coform  char[][]  format of generic term added to coeffn
**
**  Returned (function value):
**             int      -1 = unrecognized term
**                       0 = no room for new generic term
**                      >0 = term number (1...)
**
**  Called:  tptScomp, tptTrmgen
**
**  Last revision:   29 August 1999
**
**  Copyright 1999 P.T.Wallace.  All rights reserved.
*/

{
   int n;
   char* f;



/* Search table of coefficient names for match. */
   for ( n = 0; n < *nterms; n++ ) {
      if ( ! tptScomp ( name, coeffn [ n ] ) ) return ++n;
   }

/* Not matched: see whether the term has a valid generic name. */
   if ( ( f = tptTrmgen ( name ) ) == NULL ) return -1;

/* New generic term recognized: exit if no room for it. */
   if ( *nterms >= ntroom - 1 ) return 0;

/* Enough room: add it to the repertoire, then exit. */
   n = (*nterms)++;
   strcpy ( coeffn [ n ], name );
   strcpy ( coform [ n ], f );
   return *nterms;

}
