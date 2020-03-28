#ifndef TPTIC
#define TPTIC

#ifdef __cplusplus
extern "C" {
#endif

/*
**  - - - - - - - - -
**   t p t l i b . h
**  - - - - - - - - -
**
**  TPT macros and function prototypes.
**
**  Last revision:   27 October 2000
**
**  Copyright P.T.Wallace.  All rights reserved.
*/

/*
**  ------
**  Macros
**  ------
*/

/* Mount codes */
#define EQUAT 1               /* equatorial */
#define ALTAZ 2               /* altazimuth */
#define GIMBAL 3              /* generalized gimbal */

/*
** -------------------
** Function prototypes
** -------------------
*/
int tcsAddtrm ( int, int, char[9], double, int*, char[][9], char[][9], 
                int[], double[] );
int tcsPup ( int, int[], double[], int, int, int, char[][9], char[][9],
             int, double[3][3], double, double, int, double[3],
             double, double*, double*, double*,
             double*, double*, double[3][3] );

/* Directly-referenced TPOINT functions. */
int tptMinit ( int, int, int[], int*, int*, int*, char[][9] );
int tptTerm ( char*, int, int*, char[][9], char[][9] );
int tptPterm ( int, double, double[], int, double[3], double, double,
               int, double[3][3], int, int, int, char[][9], char[][9],
               int, double[3], double[8] );

#ifdef __cplusplus
}
#endif

#endif
