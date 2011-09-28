/*
 * Switchable support for double and extended precision.
 *
 * We define a type "myfloat" which will be either double or long
 * double.
 *
 * Use of double precision constants (like 0.5) in extended
 * precision-based codes can result in limited accuracy.  In extended
 * precision, we need 0.5L.  For this purpose we define FPX where X is
 * 0, 1, 2, 3, 4, 5 and "half".
 *
 */

#ifndef MULTIPREC_H
#define MULTIPREC_H


/* Get appropriate math.h lib */
#ifdef EXTENDEDPRECISION
#include <tgmath.h>
#else
#include <math.h>
/* (tgmath.h would work if available.  Could probably do something
   Windows specific here.) */
#endif


/* Define my_float */
#ifdef EXTENDEDPRECISION
typedef long double myfloat;
#else
typedef double myfloat;
#endif


/* TODO */
/* Using this instead of macro: works with different FP. */
/*const myfloat CONST_PI = 3.14159265358979323846L;*/


/* Floating point constants, add more if you need them */
#ifdef EXTENDEDPRECISION
#define FP0 0.0L
#define FP1 1.0L
#define FP2 2.0L
#define FP3 3.0L
#define FP4 4.0L
#define FP5 5.0L
#define FPhalf 0.5L
#else
#define FP0 0.0
#define FP1 1.0
#define FP2 2.0
#define FP3 3.0
#define FP4 4.0
#define FP5 5.0
#define FPhalf 0.5
#endif


#endif  /* MULTIPREC_H */
