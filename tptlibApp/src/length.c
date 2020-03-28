#include <ctype.h>
int tptLength ( char *string )
/*
**  - - - - - - - - - -
**   t p t L e n g t h
**  - - - - - - - - - -
**
**  Length of a string excluding any trailing whitespace.
**
**  Given:
**     *string      char     the string
**
**  Returned (function value):
**                  int      length excluding trailing spaces
**
**  The minimum length returned is 0.
**
**  Last revision:   13 January 1996
**
**  Copyright 1996 P.T.Wallace.  All rights reserved.
*/
{
   int l, i, c;

   l = 0;
   for ( i = 0;
         ( c = (int) string [ i++ ] );
         l = isspace ( c ) ? l : i );
   return l;
}
