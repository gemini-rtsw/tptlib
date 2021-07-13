#include <string.h>

int tptLength ( char* );

int tptScomp ( char *string1, char *string2 )
/*
**  - - - - - - - - -
**   t p t S c o m p
**  - - - - - - - - -
**
**  Compare two strings ignoring trailing whitespace.
**
**  Given:
**     *string1     char     first string
**     *string2     char     second string
**
**  Returned (function value):
**                  int      0 = strings are the same apart
**                               from any trailing whitespace
**
**  Called:  tptLength
**
**  Last revision:   10 August 1996
**
**  Copyright 1996 P.T.Wallace.  All rights reserved.
*/
{
   int l1, l2;

   l1 = tptLength ( string1 );
   l2 = tptLength ( string2 );
   return strncmp ( string1, string2, ( l1 > l2 ) ? l1 : l2 );
}
