#include <stdio.h>

char* tptTrmgen ( char *name )
/*
**  - - - - - - - - - -
**   t p t T r m g e n
**  - - - - - - - - - -
**
**  Check a generic term name for validity and return a pointer to a
**  string which identifies the format.
**
**  Given (argument):
**     name    char[]   name of term
**
**  Returned (function value):
**             char*    NULL = not a valid generic term name (Note 1)
**                      else = pointer to format string (Note 2)
**
**  Notes
**
**  1  There are three types of generic term, called P, H and A
**     respectively.
**
**     P terms (polynomials) apply adjustments to a "result" coordinate
**     which are proportional to the product of two "input" coordinates
**     each raised to an integer power greater than or equal to zero.
**
**     H terms (harmonics) apply adjustments to a "result" coordinate
**     which are proportional to the product of sines or cosines of two
**     "input" coordinates.
**
**     A terms (auxiliaries) apply adjustments to a "result" coordinate
**     which are proportional to an auxiliary reading.
**
**     A valid generic term name is one of the following:
**
**        Prc[i[c[i]]]
**        Hrfc[i][fc[i]]
**        Hrfci[i[i]]
**        Hrfcfcii
**        Hrfciifc
**        Anr
**
**     Where:
**
**        P = the character P
**        H = the character H
**        A = the character A
**        r = a single-letter code specifying the "result" coordinate:
**               H    result is in hour angle
**               X    result is East-West on the sky
**               D    result is in declination
**               U    result moves polar axis up/down
**               L    result moves polar axis left/right
**               P    result changes HA/Dec nonperpendicularity
**               A    result is in azimuth
**               S    result is left-right on the sky
**               Z    result is in zenith distance (sky)
**               E    result is in elevation (mount)
**               N    result moves azimuth axis north/south
**               W    result moves azimuth axis east/west
**               V    result changes Az/El nonperpendicularity
**        c = a single-letter code specifying an "input" coordinate:
**               H    hour angle
**               D    declination
**               A    azimuth
**               Z    zenith distance
**               E    elevation
**               Q    parallactic angle
**        i = a decimal digit 0-9
**        f = function: 'S' for sine or 'C' for cosine
**        n = a decimal digit 1-3
**
**     A missing i, ci or fci defaults to unity.
**
**     Examples of valid generic term names:
**        PDH       = HA adjustment to Dec
**        PZZ2      = ZD**2 adjustment to ZD
**        PXH3D2    = (HA**3)*(DEC**2) adjustment east-west
**        HHSH      = sin(HA) adjustment to HA
**        HACA12    = cos(12*AZ) adjustment to Az
**        HESA3     = sin(3*AZ) adjustment to El
**        HDSH2CD7  = sin(2*HA)*cos(7*DEC) adjustment to Dec
**        A1E       = auxiliary reading #1 adjustment to mount elevation
**
**  2  The format string pointed to by this function can be any of those
**     given in the format array, below.  These are used by the tptTerm
**     routine to reduce the amount of string parsing that goes on during
**     intensive processing.
**
**  3  This routine is self-contained and can be used outside the
**     TPOINT package.
**
**  Last revision:   17 August 1999
**
**  Copyright 1999 P.T.Wallace.  All rights reserved.
*/

#define TRUE 1
#define FALSE 0

{
   int nform, i, c1, c2, ok, ncom;
   char* f;

/* Table of valid formats */
   char format [ ] [ 9 ] = {
      "PRCICI",
      "PRCIC",
      "PRCI",
      "PRC",
      "HRFCIII",
      "HRFCII",
      "HRFCI",
      "HRFC",
      "HRFCIFCI",
      "HRFCIFC",
      "HRFCFCI",
      "HRFCFC",
      "HRFCIIFC",
      "HRFCFCII",
      "ANR"
   };

/* Number of valid formats */
   int nforms = sizeof ( format ) / sizeof ( format [ 0 ] );

/* Table of valid format/name character options */
   int com [ ] [ 2 ] = {
      { 'P','P' },
      { 'H','H' },
      { 'A','A' },
      { 'R','H' },
      { 'R','X' },
      { 'R','D' },
      { 'R','U' },
      { 'R','L' },
      { 'R','P' },
      { 'R','A' },
      { 'R','S' },
      { 'R','Z' },
      { 'R','E' },
      { 'R','N' },
      { 'R','W' },
      { 'R','V' },
      { 'C','H' },
      { 'C','D' },
      { 'C','A' },
      { 'C','Z' },
      { 'C','E' },
      { 'C','Q' },
      { 'I','0' },
      { 'I','1' },
      { 'I','2' },
      { 'I','3' },
      { 'I','4' },
      { 'I','5' },
      { 'I','6' },
      { 'I','7' },
      { 'I','8' },
      { 'I','9' },
      { 'F','S' },
      { 'F','C' },
      { 'N','1' },
      { 'N','2' },
      { 'N','3' },
      { ' ',' ' }
   };

/* Number of valid format/name character pairs */
   int ncoms = sizeof ( com ) / sizeof ( com [ 0 ] );



/* For each allowed format ... */
   for ( nform = 0; nform < nforms; nform++ ) {

   /* Current format. */
      f = format [ nform ];

   /* For each character position in format/name... */
      for ( i = 0; i < 9; i++ ) {

      /* Format and name characters. */
         c1 = (int) f [ i ];
         if ( c1 == ' ' ) c1 = '\0';
         c2 = (int) name [ i ];
         if ( c2 == ' ' ) c2 = '\0';

      /* If both are null or space we have a match. */
         if ( c1 == '\0' && c2 == '\0' ) return f;

      /* Not matched yet:  set flag. */
         ok = FALSE;

      /* For each option in com table... */
         for ( ncom = 0; ncom < ncoms; ncom++ ) {

         /* Is this one character consistent with both name & format? */
            if ( c1 == com [ ncom ] [ 0 ] &&
                 c2 == com [ ncom ] [ 1 ] ) {

            /* Yes, name and format consistent so far. */
               ok = TRUE;
               break;
            }

         /* Next com option. */
         }

      /* Give up if we didn't match this character. */
         if ( ! ok ) break;

      /* Otherwise move on to next character in format/name. */
      }

   /* Try next format. */
   }

/* Not a valid generic term name: exit. */
   return NULL;
}
