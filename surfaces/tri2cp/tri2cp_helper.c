/********************************************************************
 * TRI2CP_HELPER.c
 * Takes triangles (as faces and vertices) and returns closest points.
 * Written as a mex interface for matlab/octave.
 *
 * Authors: Colin Macdonald with routines by Steve Ruuth
 *
 * To compile, type this in matlab:
 *
 * mex helper_tri2cp.c
 *
 * (or "mex -O CFLAGS='\$CFLAGS -Wall' helper_tri2cp.c")
 *
 * Maintainence: recompiling this on many different OSes is a hassle:
 * its better to remove intelligence from this and put it in the M
 * file instead.
 *
 * Limitations:
 *
 * Can use extended precision internally: may not be completely tested
 * yet.  Results are copied back into doubles to return to matlab.
 *
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
/*#include <search.h>*/
#include <stddef.h>    /* hash needs this for offsetof */
#include <string.h>
#include <errno.h>
/* don't include math.h here, see multiprec.h */

#include <limits.h>

#include <mex.h>

#include "uthash.h"

/* Use extended precision or not */
/*#define EXTENDEDPRECISION*/
#include "multiprec.h"

/* What type of boundingSpheres */
#define boundingSpheres boundingSpheresOptimal
/*#define boundingSpheres boundingSpheresCentroid */

/* Also form 3D matrices: uses a lot of memory, just for debugging.
   RELPT must be specified carefully to use this. */
/*#define ALSOTHREEDMATRICES
  #include "debug_mat3d.h"*/



/***************************************************************************/



/* use 1st one for calls that set errno, 2nd otherwise */
/*#define myerr_nix(s) { perror((s)); exit(EXIT_FAILURE); }
  #define myerr(s) { fprintf(stderr, "%s\n", (s)); exit(EXIT_FAILURE); }*/
#define myerr_nix(s) { mexErrMsgTxt(s); exit(EXIT_FAILURE); }
#define myerr(s) { mexErrMsgTxt(s); exit(EXIT_FAILURE); }


typedef struct {
  myfloat x, y, z;
  myfloat confidence;
  myfloat intensity;
} struct_vertex;

typedef struct {
  long v1, v2, v3;
} struct_face;

/* Struct used for the hashtable.  Don't need the (x,y,z) but it might
   be useful for debugging. */
typedef struct {
  unsigned long long key;
  int i, j, k;  /* the hash table key */
  myfloat cpx, cpy, cpz;
  myfloat dd;
  myfloat x, y, z;
  long which_face;
  UT_hash_handle hh;
} struct_cp;



/* Global variables */
myfloat RELPTX, RELPTY, RELPTZ;
myfloat DX;
myfloat BANDWIDTH;
int DEBUG_LEVEL;  /* messages with smaller level will display */
unsigned long number_vertices, number_faces;

/* TODO: probably no need to be global */
/* TODO: REMOVE? */
unsigned long ExpectedGridSize, HashTableSize;

/* new hash */
struct_cp *newhashlist = NULL;

/*struct_cp **gridPtList;
  char **keyList;*/

unsigned long numgridpts;

struct_vertex *vertex;
struct_face   *face;
myfloat *global_Sphere_Cx;
myfloat *global_Sphere_Cy;
myfloat *global_Sphere_Cz;
myfloat *global_Sphere_w;

struct {int i; int j; int k; } *lookup_key;
unsigned int keylen;


/*
 * A debugging print function
 */
#define BUFMAXSIZE 1000
void dbg_printf(int level, const char *fmt, ...)
{
  char buf[BUFMAXSIZE];
  va_list ap;

  if (level <= DEBUG_LEVEL) {
    va_start(ap,fmt);
    vsnprintf(buf, BUFMAXSIZE-1, fmt, ap);
    va_end(ap);
    /*mexPrintf("  %s: %s", mexFunctionName(), buf);*/
    mexPrintf("  %s: %s", "MEXDEBUG", buf);
  }
}

/*
#define hashAddGridPt hashAddGridPt_MF
#define hashFind hashFind_MF
*/
#define hashAddGridPt hashAddGridPt_LONG
#define hashFind hashFind_LONG


/* Needs to be 0.5*10^k  */
#define MAXIND (50000UL)
#define KEYMULTFAC (2*MAXIND)
#define KEYSHIFT (MAXIND-1)


/*
 *
 */
void hashInit()
{
  double A, B;

  /* empty hash list */
  newhashlist = NULL;

  /* for the multifield key, need these global vars */
  keylen = offsetof(struct_cp, k) + sizeof(int) - offsetof(struct_cp, i);
  lookup_key = malloc(sizeof(*lookup_key));
  memset(lookup_key, 0, sizeof(*lookup_key));  /* zero fill */

  /* Some checking for bounds.  Can't check if an integer is larger
     than LONG_MAX so use floating points. */
  /*
  dbg_printf(10, "KEYMULTFAC^3 = %llu, LONG_MAX = %llu,%llu,%llu,%llu\n", \
	     KEYMULTFAC*KEYMULTFAC*KEYMULTFAC, ULONG_MAX, ULLONG_MAX, LONG_MAX, LLONG_MAX);
  dbg_printf(10, "KEYMULTFAC = %d, KEYMULTFAC^3 = %ld, KEYMULTFAC^3 = %llu\n", \
	     KEYMULTFAC, KEYMULTFAC*KEYMULTFAC*KEYMULTFAC, KEYMULTFAC*KEYMULTFAC*KEYMULTFAC);
  */
  A = (double) KEYMULTFAC;
  A = A*A*A;
  B = (double) ULLONG_MAX;
  dbg_printf(10, "(in FP) KEYMULTFAC^3 = %g < ULLONG_MAX = %g\n", A, B);
  dbg_printf(99, "sizeof(unsigned int)=%d\n", sizeof(unsigned int));
  dbg_printf(99, "sizeof(unsigned long)=%d\n", sizeof(unsigned long));
  dbg_printf(99, "sizeof(unsigned long long)=%d\n", sizeof(unsigned long long));
  if (0.95*A > B) {
    myerr("not enough space in unsigned long long int");
  }
}

/*
 *
 */
void hashDeleteAll()
{
  struct_cp *s, *tmp;

  HASH_ITER(hh, newhashlist, s, tmp) {
    HASH_DEL(newhashlist, s);  /* delete; s advances to next */
    free(s);
  }

  free(lookup_key);
}


/*
 * Add a gridpoint to the hash table.
 * Integer key version.
 */
void hashAddGridPt_LONG(int i, int j, int k, \
			myfloat dd, myfloat cpx, myfloat cpy, myfloat cpz, \
			myfloat x, myfloat y, myfloat z, long which_face)
{
  struct_cp *gridpt;

  gridpt = malloc(sizeof(struct_cp));

  /* this gives a unique (long) integer key for i,j,k \in [-4999 4999]
     which we can then hash.  This seems to be faster (maybe 4%)
     compared to using the "multiple contiguous fields" of struct_cp
     in uthash.  Same formula must be used below in hashFind() */
  if ( (abs(i) > MAXIND) | (abs(j) > MAXIND) | (abs(k) > MAXIND) ) {
    dbg_printf(0, "MAXIND=%d, (i,j,k)=(%d,%d,%d)\n", MAXIND, i, j, k);
    myerr("indices are too big");
  }
  gridpt->key = i*(KEYMULTFAC*KEYMULTFAC) + (j+KEYSHIFT)*KEYMULTFAC + (k+KEYSHIFT);

  gridpt->i = i;
  gridpt->j = j;
  gridpt->k = k;
  gridpt->cpx = cpx;
  gridpt->cpy = cpy;
  gridpt->cpz = cpz;
  gridpt->dd = dd;
  gridpt->x = x;
  gridpt->y = y;
  gridpt->z = z;
  gridpt->which_face = which_face;

  HASH_ADD( hh, newhashlist, key, sizeof(unsigned long long), gridpt );
  /* HASH_ADD_INT( newhashlist, key, gridpt ); */
}


/*
 * Return a pointer to the grid point or NULL if cannot be found.
 * Integer key version.
 */
struct_cp *hashFind_LONG(int i, int j, int k)
{
  struct_cp *s;
  unsigned long long findkey;
  findkey = i*(KEYMULTFAC*KEYMULTFAC) + (j+KEYSHIFT)*KEYMULTFAC + (k+KEYSHIFT);
  HASH_FIND( hh, newhashlist, &findkey, sizeof(unsigned long long), s );
  /* HASH_FIND_INT( newhashlist, &findkey, s ); */
  return(s);
}


/*
 * Add a gridpoint to the hash table.
 * Multifield key version.
 */
void hashAddGridPt_MF(int i, int j, int k, \
		   myfloat dd, myfloat cpx, myfloat cpy, myfloat cpz, \
		   myfloat x, myfloat y, myfloat z)
{
  struct_cp *gridpt;

  gridpt = malloc(sizeof(struct_cp));
  memset(gridpt, 0, sizeof(struct_cp));  /* zero fill */

  gridpt->i = i;
  gridpt->j = j;
  gridpt->k = k;
  gridpt->cpx = cpx;
  gridpt->cpy = cpy;
  gridpt->cpz = cpz;
  gridpt->dd = dd;
  gridpt->x = x;
  gridpt->y = y;
  gridpt->z = z;

  /* add our structure to the hash table */
  HASH_ADD( hh, newhashlist, i, keylen, gridpt);
}



/*
 * Return a pointer to the grid point or NULL if cannot be found.
 * Multifield key version
 */
struct_cp *hashFind_MF(int i, int j, int k)
{
  struct_cp *s;
  lookup_key->i = i;
  lookup_key->j = j;
  lookup_key->k = k;
  HASH_FIND( hh, newhashlist, &lookup_key->i, keylen, s );
  return(s);
}



/*
 * functions to malloc and free space for the ply data
 */
void mallocPlyData(long numv, long numf)
{
  if ((vertex = (struct_vertex *)			\
       malloc(numv * sizeof(struct_vertex))) == NULL)
    myerr_nix("Error mallocing ply data")

  if ((face = (struct_face *) malloc(numf * sizeof(struct_face))) == NULL)
    myerr_nix("Error mallocing ply data")

  if ((global_Sphere_Cx = (myfloat *) malloc(numf * sizeof(myfloat))) == NULL)
    myerr_nix("Error mallocing ply data")
  if ((global_Sphere_Cy = (myfloat *) malloc(numf * sizeof(myfloat))) == NULL)
    myerr_nix("Error mallocing ply data")
  if ((global_Sphere_Cz = (myfloat *) malloc(numf * sizeof(myfloat))) == NULL)
    myerr_nix("Error mallocing ply data")

  if ((global_Sphere_w = (myfloat *) malloc(numf * sizeof(myfloat))) == NULL)
    myerr_nix("Error mallocing ply data")
}

void freePlyData()
{
  free(vertex);
  free(face);
  free(global_Sphere_Cx);
  free(global_Sphere_Cy);
  free(global_Sphere_Cz);
  free(global_Sphere_w);
}




/*
 * Copy Faces and Vertices from the matlab to the global variables.
 * (probably this copy is not necessary because we don't change this
 * data but this way keeps matlab arrays maximally separate from
 * native ones)
 */
void initShapeFromMatlabArray(double *Faces, long numF, double *Vertices, long numV)
{
  unsigned long i;

  mallocPlyData(numV, numF);

  /*mexPrintf("%ld %ld %ld %ld\n", i, vertex[i].x, vertex[i].y, vertex[i].z);*/
  for (i=0; i<numV; i++) {
    vertex[i].x = Vertices[0*numV+i];
    vertex[i].y = Vertices[1*numV+i];
    vertex[i].z = Vertices[2*numV+i];
    /*mexPrintf("%ld %g %g %g\n", i, vertex[i].x, vertex[i].y, vertex[i].z);
      if (vertex[i].z < -FP1 | vertex[i].z > FP1) { printf("HELP!\n"); error(-1); }*/
  }
  for (i=0; i<number_faces; i++) {
    /* -1 here because Faces is from Matlab */
    face[i].v1 = (long) Faces[0*numF+i] - 1;
    face[i].v2 = (long) Faces[1*numF+i] - 1;
    face[i].v3 = (long) Faces[2*numF+i] - 1;
    /*mexPrintf("%ld %ld %ld %ld\n", i, face[i].v1, face[i].v2, face[i].v3);*/
  }
}



/*
 * int away from zero: like ceil, floor but always moves away from
 * zero.  Special case: iafz(+-0.0) = 0
 * TODO: should that be x==0L in extended precision? (FP0?)
 */
int iafz(myfloat x)
{
  if (x==0) {
    return 0;
  } else if (x > 0) {
    return (int) ceil(x);
  } else {
    return (int) floor(x);
  }
  /* wrong b/c of signed zero */
  /*return( (int) (x + signbit(x)) ); */

  /* this approach cannot be trusted for large x */
  /*int signum;
    signum = (x > 0) - (x < 0);
    return( (int) (x + signum) ); */
}



/*
 * Functions for max/min floating point.
 * (these were inline in C99 but not sure how to do that in ansi C)
 */
myfloat float_max(myfloat a, myfloat b)
{
  return (a > b) ? a : b;
}
myfloat float_min(myfloat a, myfloat b)
{
  return (a < b) ? a : b;
}



/*
 * find optimal bounding spheres for the triangles
 */
void boundingSpheresOptimal()
{
  /* Wikipedia: The useful minimum bounding circle of three points is
     defined either by the circumcircle (where three points are on the
     minimum bounding circle) or by the two points of the longest side
     of the triangle (where the two points define a diameter of the
     circle). */
  unsigned long i;
  myfloat ang[3];
  myfloat m[3];
  myfloat a[3], b[3], c[3];
  /* TODO: change to arrays */
  myfloat xba,yba,zba, xca,yca,zca;
  myfloat R;
  myfloat minang;
  unsigned short minangi;

  for (i=0; i<number_faces; i++) {  /* for each triangle */
    a[0] = vertex[face[i].v1].x;
    a[1] = vertex[face[i].v1].y;
    a[2] = vertex[face[i].v1].z;
    b[0] = vertex[face[i].v2].x;
    b[1] = vertex[face[i].v2].y;
    b[2] = vertex[face[i].v2].z;
    c[0] = vertex[face[i].v3].x;
    c[1] = vertex[face[i].v3].y;
    c[2] = vertex[face[i].v3].z;

    xba = b[0] - a[0];
    yba = b[1] - a[1];
    zba = b[2] - a[2];
    xca = c[0] - a[0];
    yca = c[1] - a[1];
    zca = c[2] - a[2];

    /* anga = (C-A)'*(B-A)
       angb = (A-B)'*(C-B)
       angc = (A-C)'*(B-C) */
    ang[0] = xca*xba + yca*yba + zca*zba;
    ang[1] = (-xba) * (c[0] - b[0]) +  \
             (-yba) * (c[1] - b[1]) +  \
             (-zba) * (c[2] - b[2]);
    ang[2] = (-xca) * (b[0] - c[0]) +  \
             (-yca) * (b[1] - c[1]) +  \
             (-zca) * (b[2] - c[2]);

    /* ang[0] = (c[0] - a[0]) * (b[0] - a[0]) + \ */
    /*          (c[1] - a[1]) * (b[1] - a[1]) + \ */
    /*          (c[2] - a[2]) * (b[2] - a[2]); */
    /* ang[1] = (a[0] - b[0]) * (c[0] - b[0]) + \ */
    /*          (a[1] - b[1]) * (c[1] - b[1]) + \ */
    /*          (a[2] - b[2]) * (c[2] - b[2]); */
    /* ang[2] = (a[0] - c[0]) * (b[0] - c[0]) + \ */
    /*          (a[1] - c[1]) * (b[1] - c[1]) + \ */
    /*          (a[2] - c[2]) * (b[2] - c[2]); */

    minang = ang[0];
    minangi = 0;
    if (ang[1] < minang) {
      minang = ang[1];
      minangi = 1;
    }
    if (ang[2] < minang) {
      minang = ang[2];
      minangi = 2;
    }

    if (minang < 0) {
      /* obtuse triangle, use longest side */
      if (minangi == 0) {
	/* m = (B+C)/2 */
	/* TODO: better?: m = (B-C)/2 + C */
	m[0] = (b[0] + c[0]) / FP2;
	m[1] = (b[1] + c[1]) / FP2;
	m[2] = (b[2] + c[2]) / FP2;
	R = ( (b[0] - c[0]) * (b[0] - c[0]) +
	      (b[1] - c[1]) * (b[1] - c[1]) +
	      (b[2] - c[2]) * (b[2] - c[2]) ) / FP4;
      } else if (minangi == 1) {
	/* m = (A+C)/2 */
	m[0] = (a[0] + c[0]) / FP2;
	m[1] = (a[1] + c[1]) / FP2;
	m[2] = (a[2] + c[2]) / FP2;
	R = ( xca*xca + yca*yca + zca*zca ) / FP4;
      } else if (minangi == 2) {
	/* m = (A+B)/2 */
	m[0] = (a[0] + b[0]) / FP2;
	m[1] = (a[1] + b[1]) / FP2;
	m[2] = (a[2] + b[2]) / FP2;
	R = ( xba*xba + yba*yba + zba*zba ) / FP4;
      } else {
	myerr("no case");
      }
      R = sqrt(R);
    } else {
      /* acute triangle, use circumcenter

	 uses code from
	 http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html

          |                                                           |
          | |c-a|^2 [(b-a)x(c-a)]x(b-a) + |b-a|^2 (c-a)x[(b-a)x(c-a)] |
          |                                                           |
      r = -------------------------------------------------------------,
                               2 | (b-a)x(c-a) |^2

              |c-a|^2 [(b-a)x(c-a)]x(b-a) + |b-a|^2 (c-a)x[(b-a)x(c-a)]
      m = a + ---------------------------------------------------------.
                                 2 | (b-a)x(c-a) |^2
      */

      myfloat balength, calength;
      myfloat xcrossbc, ycrossbc, zcrossbc;
      myfloat denominator;

      /* Squares of lengths of the edges incident to `a'. */
      balength = xba * xba + yba * yba + zba * zba;
      calength = xca * xca + yca * yca + zca * zca;

      /* Cross product of these edges. */
      /* Take your chances with floating-point roundoff.  (code on
	 webpage could also use a more robust library for this). */
      xcrossbc = yba * zca - yca * zba;
      ycrossbc = zba * xca - zca * xba;
      zcrossbc = xba * yca - xca * yba;

      /* Calculate the denominator of the formulae. */
      /* (FPhalf is 0.5) */
      denominator = FPhalf / (xcrossbc * xcrossbc + ycrossbc * ycrossbc +
			      zcrossbc * zcrossbc);

      /* Calculate offset (from `a') of circumcenter. */
      m[0] = ((balength * yca - calength * yba) * zcrossbc -
	      (balength * zca - calength * zba) * ycrossbc) * denominator;
      m[1] = ((balength * zca - calength * zba) * xcrossbc -
	      (balength * xca - calength * xba) * zcrossbc) * denominator;
      m[2] = ((balength * xca - calength * xba) * ycrossbc -
	      (balength * yca - calength * yba) * xcrossbc) * denominator;
      /*R = sqrt(xcirca*xcirca + ycirca*ycirca + zcirca*zcirca) */
      R = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
      m[0] += a[0];
      m[1] += a[1];
      m[2] += a[2];
    }
    global_Sphere_Cx[i] = m[0];
    global_Sphere_Cy[i] = m[1];
    global_Sphere_Cz[i] = m[2];
    R += BANDWIDTH;
    global_Sphere_w[i] = R;
    /*printf("A=[%.20g,%.20g,%.20g]'\n", a[0],a[1],a[2]);
      printf("B=[%.20g,%.20g,%.20g]'\n", b[0],b[1],b[2]);
      printf("C=[%.20g,%.20g,%.20g]'\n", c[0],c[1],c[2]);
      printf("Cen=[%.20g,%.20g,%.20g]'\n", global_Sphere_Cx[i], global_Sphere_Cy[i], global_Sphere_Cz[i]);
      printf("rad=%.20g\n", global_Sphere_w[i]);*/
  } /* end loop over triangles */
}



/*
 * a helper function for boundingSpheresCentroid() below
 */
myfloat distance(unsigned long i, unsigned long j)
{
  myfloat sqrd;

  sqrd = (vertex[i].x-global_Sphere_Cx[j])*(vertex[i].x-global_Sphere_Cx[j])
       + (vertex[i].y-global_Sphere_Cy[j])*(vertex[i].y-global_Sphere_Cy[j])
       + (vertex[i].z-global_Sphere_Cz[j])*(vertex[i].z-global_Sphere_Cz[j]);
  return(sqrt(sqrd));
}

/*
myfloat distance2(myfloat x, myfloat y, myfloat z, long j)
{
  myfloat dd;

  dd = (x-global_Sphere_Cx[j])*(x-global_Sphere_Cx[j])
     + (y-global_Sphere_Cy[j])*(y-global_Sphere_Cy[j])
     + (z-global_Sphere_Cz[j])*(z-global_Sphere_Cz[j]);
  return(dd);
}
*/



/*
 * compute bounding spheres using the centroid, easy algorithm but not
 * a very tight bound.
 */
void boundingSpheresCentroid()
{
  unsigned long i;
  myfloat R;

  for (i=0; i<number_faces; i++) {
    global_Sphere_Cx[i] = (vertex[face[i].v1].x+vertex[face[i].v2].x+vertex[face[i].v3].x)/FP3;
    global_Sphere_Cy[i] = (vertex[face[i].v1].y+vertex[face[i].v2].y+vertex[face[i].v3].y)/FP3;
    global_Sphere_Cz[i] = (vertex[face[i].v1].z+vertex[face[i].v2].z+vertex[face[i].v3].z)/FP3;
    R = float_max(distance(face[i].v1,i), float_max(distance(face[i].v2,i),distance(face[i].v3,i)));
    R += BANDWIDTH;
    global_Sphere_w[i]  = R;
  }
}



/*
 * Project the point (c1,c2,c3) onto the line segment specified by
 * (p1,p2,p3) and (q1,q2,q3).
 * (code by Steve Ruuth)
 */
void ProjectOnSegment(myfloat *c1, myfloat *c2, myfloat *c3, \
		      myfloat p1, myfloat p2, myfloat p3, \
		      myfloat q1, myfloat q2, myfloat q3)
{
  myfloat lambda, lambda_star;
  myfloat cmp1, cmp2, cmp3;
  myfloat qmp1, qmp2, qmp3;

  cmp1 = *c1-p1;
  cmp2 = *c2-p2;
  cmp3 = *c3-p3;
  qmp1 =  q1-p1;
  qmp2 =  q2-p2;
  qmp3 =  q3-p3;

  lambda = (cmp1*qmp1+cmp2*qmp2+cmp3*qmp3)/(qmp1*qmp1+qmp2*qmp2+qmp3*qmp3);
  lambda_star = float_max(FP0, float_min(lambda, FP1));

  *c1 = p1+lambda_star*qmp1;
  *c2 = p2+lambda_star*qmp2;
  *c3 = p3+lambda_star*qmp3;
}



/*
 * Closest point and distance from a point (a1,a2,a3) to a triangle
 * indexed by face_index.  Returns the *squared* distance and the
 * closest point in (c1,c2,c3).  Uses global vars `face' and
 * `vertex'.  (Based on code by Steve Ruuth)
 * TODO: this code may not be robust to degenerate triangles (line
 * segments and points).  More testing required.
 */
myfloat FindClosestPointToOneTri(myfloat a1, myfloat a2, myfloat a3, \
				 unsigned long face_index, \
				 myfloat *c1, myfloat *c2, myfloat *c3)
{
  unsigned long index_p, index_q, index_r;
  myfloat a11, a12, a22, b1, b2;
  myfloat i11, i12, i22, factor;
  myfloat q1, q2, q3;
  myfloat r1, r2, r3;
  myfloat lambda, mu;
  myfloat dd;

  /* obtain the indices to the three vertices */
  index_p = face[face_index].v1;
  index_q = face[face_index].v2;
  index_r = face[face_index].v3;

  /* translate so the p is at the origin */
  a1 -= vertex[index_p].x;
  a2 -= vertex[index_p].y;
  a3 -= vertex[index_p].z;
  q1 = vertex[index_q].x-vertex[index_p].x;
  q2 = vertex[index_q].y-vertex[index_p].y;
  q3 = vertex[index_q].z-vertex[index_p].z;
  r1 = vertex[index_r].x-vertex[index_p].x;
  r2 = vertex[index_r].y-vertex[index_p].y;
  r3 = vertex[index_r].z-vertex[index_p].z;

  /* evaluate the various matrix entries */
  a11 = q1*q1+q2*q2+q3*q3;
  a12 = q1*r1+q2*r2+q3*r3;
  a22 = r1*r1+r2*r2+r3*r3;
  b1  = a1*q1+a2*q2+a3*q3;
  b2  = a1*r1+a2*r2+a3*r3;

  /* find the inverse matrix and solve for lambda and mu */
  factor = FP1/(a11*a22-a12*a12);
  i11 = a22*factor;
  i12 =-a12*factor;
  i22 = a11*factor;
  lambda = i11*b1+i12*b2;
  mu     = i12*b1+i22*b2;
  *c1 = lambda*q1+mu*r1;
  *c2 = lambda*q2+mu*r2;
  *c3 = lambda*q3+mu*r3;

  if ((lambda<0) && (mu<0) && (lambda+mu<=1)) {
    *c1 = *c2 = *c3 = FP0;
  } else if ((lambda>=0) && (mu<0) && (lambda+mu<=1)) {
    ProjectOnSegment(c1,c2,c3,FP0,FP0,FP0,q1,q2,q3);
  } else if ((lambda>=0) && (mu<0) && (lambda+mu>1)) {
    *c1 = q1;
    *c2 = q2;
    *c3 = q3;
  } else if ((lambda>=0) && (mu>=0) && (lambda+mu>1)) {
    ProjectOnSegment(c1,c2,c3,q1,q2,q3,r1,r2,r3);
  } else if ((lambda<0) && (mu>=0) && (lambda+mu>1)) {
    *c1 = r1;
    *c2 = r2;
    *c3 = r3;
  } else if ((lambda<0) && (mu>=0) && (lambda+mu<=1)) {
    ProjectOnSegment(c1,c2,c3,r1,r2,r3,FP0,FP0,FP0);
  } else if ((lambda>=0) && (mu>=0) && (lambda+mu<=1)) {
    /* TODO: do nothing, what case is this? */
  } else {
    dbg_printf(0, "Error: non-enumerated case, can this happen?\n");
#ifdef EXTENDEDPRECISION
    dbg_printf(0, "lambda mu %Lg %Lg\n", lambda, mu);
#else
    dbg_printf(0, "lambda mu %g %g\n", lambda, mu);
    dbg_printf(0, "factor %g\n", factor);
    dbg_printf(0, "a11,a22,a12=%g,%g,%g\n", a11, a22, a12);
    dbg_printf(0, "det?=%g\n", a11*a22 - a12*a12);
#endif
    myerr("Error in CP to tri, unanticipated case");
  }

  /* Calculate distance */
  /* Note: dd is dist squared! */
  /* TODO: HORRIBLE HACK 2010-07-28, this was to deal with a ply file
     with degenerate triangles. */
  /*if (isinf(factor)) {*/
  /* TODO, just gets worse and worse, where is isinf on windows? */
  if ((factor > 1e20) || (factor < -1e20)) {
    dd = 10000.0;
    myerr("'factor' is infinite: panic!");
  } else {
    dd  = (a1-*c1)*(a1-*c1)+(a2-*c2)*(a2-*c2)+(a3-*c3)*(a3-*c3);
  }

  /*  DEBUGGING
printf("lambda mu %g %g\n",lambda,mu);
printf("q %g %g %g\n",q1,q2,q3);
printf("r %g %g %g\n",r1,r2,r3);
	printf("%g %g %g %g\n",dd, *c1, *c2, *c3);
printf("c %g %g %g\n",*c1,*c2,*c3);
printf("a11 a12 a22 %g %g %g\n",a11,a12,a22);
printf("c %g %g %g\n",*c1,*c2,*c3);

printf("p %g %g %g\n",vertex[index_p].x,vertex[index_p].y,vertex[index_p].z);
printf("a %g %g %g\n",a1,a2,a3);
printf("%g\n",dd);
for (i=0;i<1000000;i++)
{
	myfloat w1,w2;
	w1=((myfloat)(rand()%640001))/640000.;
	w2=((myfloat)(rand()%640001))/640000.;
        *c1 = w1*q1+(1-w1)*w2*r1;
        *c2 = w1*q2+(1-w1)*w2*r2;
        *c3 = w1*q3+(1-w1)*w2*r3;
        dd  = (a1-*c1)*(a1-*c1)+(a2-*c2)*(a2-*c2)+(a3-*c3)*(a3-*c3);
	printf("%g %g %g %g\n",dd, *c1, *c2, *c3);
}
  */

  /* Shift everything back */
  *c1+= vertex[index_p].x;
  *c2+= vertex[index_p].y;
  *c3+= vertex[index_p].z;
  return (dd);  /* return square distance */
}



/*
 * A global search: for a given grid point, search all triangles one
 * at a time.  This is not used as part of the algorithm, just for
 * debugging
 */
myfloat FindClosestPointsGlobally(myfloat x, myfloat y, myfloat z, \
				  myfloat *c1, myfloat *c2, myfloat *c3)
{
  unsigned long i;
  myfloat dd, dd_min;
  myfloat t1, t2, t3;

  dd_min = FindClosestPointToOneTri(x,y,z, 0, &t1, &t2, &t3);
  *c1 = t1;
  *c2 = t2;
  *c3 = t3;
  for (i=1; i<number_faces; i++) {
    dd =  FindClosestPointToOneTri(x,y,z, i, &t1, &t2, &t3);
    if (dd<dd_min) {
      dd_min = dd;
      *c1 = t1;
      *c2 = t2;
      *c3 = t3;
    }
  }
  return (dd_min);
}


/* Ruuth's algorithm as described in [MR2008,MR2009]
 *
 *  preprocessing: find the bounding spheres around each triangle.
 *  Specifically compute the centroid of the triangle.  The bounding
 *  sphere is then the maximum distance from the centroid to a vertex
 *  plus the bandwidth.  Points outside of this bounding sphere cannot
 *  be in the bandwidth of the surface.
 *
 *
 *  iterate over each triangle,
 *     find a (short) list of nearby grid points
 *     for each, compute the closest point to this triangle
 *     store the smallest
 * (based on code by Steve Ruuth, hash table stuff added by cbm)
 */
void FindClosestPointsFromTriangulation()
{
  unsigned long n;
  int i,j,k;
  int iL,iU;
  int jL,jU;
  int kL,kU;
  myfloat ih;
  myfloat dd;
  myfloat x,y,z;
  myfloat cpx, cpy, cpz;
  struct_cp *gridpt;
  unsigned long cc = 0;
  unsigned long total_count = 0;

  /* introduce dummy distance to the surface.  Outside band means it
     might as well be infinite */
#ifdef ALSOTHREEDMATRICES
  for (i=0; i<NPOINTS; i++)
    for (j=0; j<NPOINTS; j++)
      for (k=0; k<NPOINTS; k++)
	CPdd[i][j][k] = 987654321.;
#endif



  /* 1/h: inverse of h */
  ih = 1/DX;
  for (n=0; n<number_faces; n++) {
    /* a bounding cube in grid space around each bounding sphere */
    iL = iafz( (global_Sphere_Cx[n]-global_Sphere_w[n] - RELPTX) * ih );
    iU = iafz( (global_Sphere_Cx[n]+global_Sphere_w[n] - RELPTX) * ih );
    jL = iafz( (global_Sphere_Cy[n]-global_Sphere_w[n] - RELPTY) * ih );
    jU = iafz( (global_Sphere_Cy[n]+global_Sphere_w[n] - RELPTY) * ih );
    kL = iafz( (global_Sphere_Cz[n]-global_Sphere_w[n] - RELPTZ) * ih );
    kU = iafz( (global_Sphere_Cz[n]+global_Sphere_w[n] - RELPTZ) * ih );

    for (i=iL; i<=iU; i++) {
      x = ((myfloat)i)*DX + RELPTX;
      for (j=jL; j<=jU; j++) {
	y = ((myfloat)j)*DX + RELPTY;
	for (k=kL; k<=kU; k++) {
	  z = ((myfloat)k)*DX + RELPTZ;
	  total_count++;
	  if ( ( (x-global_Sphere_Cx[n])*(x-global_Sphere_Cx[n]) +
		 (y-global_Sphere_Cy[n])*(y-global_Sphere_Cy[n]) +
		 (z-global_Sphere_Cz[n])*(z-global_Sphere_Cz[n]) )
	       <= global_Sphere_w[n]*global_Sphere_w[n] ) {
	    dd = FindClosestPointToOneTri(x,y,z, n, &cpx, &cpy, &cpz);
	    cc++;

#ifdef ALSOTHREEDMATRICES
	    if (dd<CPdd[i][j][k]) {
	      CPdd[i][j][k] = dd;
	      CPx[i][j][k]  = cpx;
	      CPy[i][j][k]  = cpy;
	      CPz[i][j][k]  = cpz;
	    }
#endif

	    /* Search hashtable */
	    gridpt = hashFind(i, j, k);

	    if (gridpt == NULL) {
	      /* not found, add new one */
	      hashAddGridPt(i,j,k, dd, cpx,cpy,cpz, x,y,z, n);
	      numgridpts++;
	    } else {
	      /* already there, check if closer and update */
	      if (dd < (gridpt->dd)) {
		gridpt->dd = dd;
		gridpt->cpx = cpx;
		gridpt->cpy = cpy;
		gridpt->cpz = cpz;
		gridpt->which_face = n;
		/* no need to update i,j,k,x,y,z because they haven't changed */
		if ( (gridpt->i != i) | (gridpt->j != j) | (gridpt->k != k) | \
		     (gridpt->x != x) | (gridpt->y != y) | (gridpt->z != z) ) {
		  dbg_printf(10, "i,j,k=%d,%d,%d\n", i,j,k);
		  myerr("panic");
		}
	      }
	    }

	  } /* end inside sphere */
	}
      }
    }
  }


  dbg_printf(10, "total loops: %ld\n", total_count);
  dbg_printf(10, "made %ld calls to Find_CP_to_one_tri()\n", cc);
  dbg_printf(10, "found %ld grid points\n", numgridpts);

  /* TODO: double check only, remove? */
  /*
  unsigned long numcheck;
  numcheck = HASH_COUNT(newhashlist);
  if (numgridpts != numcheck) {
    dbg_printf(0, "not the same: %ld, %ld\n", numgridpts, numcheck);
    myerr("sanity check on number of grid points failed");
  }
  */
}




/*
 * output a grid to a file, for debugging only
 */
void outputGridToFile(const char *fname, int withPruning)
{
  unsigned long d = 0;
  FILE *fd;
  struct_cp *s;

  myerr("Deprecated");

  if ((fd = fopen(fname, "w")) == NULL) myerr_nix("Error writing griddata");

#ifdef EXTENDEDPRECISION
#define PSTR "%d %d %d %.21Le %.21Le %.21Le %.21Le %.21Le %.21Le %.21Le\n"
#else
#define PSTR "%d %d %d %.18g %.18g %.18g %.18g %.18g %.18g %.18g\n"
#endif

  d = 0;
  for(s = newhashlist; s != NULL; s = s->hh.next) {
    if ( (!withPruning) || (s->dd <= BANDWIDTH*BANDWIDTH) ) {
      fprintf(fd, PSTR, s->i,s->j,s->k, sqrt(s->dd),	\
	      s->cpx,s->cpy,s->cpz, s->x,s->y,s->z);
      d++;
    }
  }
  if ( fclose(fd) != 0 )  myerr_nix("Error closing file");
  dbg_printf(10, "Saved %ld gridpts to file\n", d);
}






/*
 * this was main() before adding mex
 */
void mainRoutine(void)
{

#ifdef ALSOTHREEDMATRICES
  /* This makes lots of assumptions about the grid.  Most importantly
     assumes relpt is the lower right-hand corner of the grid and that
     the model is centered at the origin. */
  myfloat temp;
  /*temp = min(RELPTX,RELPTY,RELPTZ);*/
  temp = (RELPTX<=RELPTY)?RELPTX:RELPTY; temp = (temp<=RELPTZ)?temp:RELPTZ;
  NPOINTS = (int) ( -(2*temp)/DX );
  dbg_printf(10, "3D matrices: assuming NPOINTS=%d\n", NPOINTS);
  dbg_printf(10, "3D matrices: mallocing...\n");
  CPx = mallocMatrix();
  CPy = mallocMatrix();
  CPz = mallocMatrix();
  CPdd = mallocMatrix();
#endif

  numgridpts = 0;


  dbg_printf(10, "running boundingSpheres()...\n");
  boundingSpheres();
  dbg_printf(10, "Finding Grid...\n");
  FindClosestPointsFromTriangulation();

  /* myfloat maxdiff = -1e42; */
  /* for (i = 0; i < numgridpts; i++) { */
  /*   myfloat dd1; */
  /*   myfloat dd2; */
  /*   int ii,jj,kk; */
  /* dd1 = gridPtList[i]->dd; */
  /* ii = gridPtList[i]->i; */
  /* jj = gridPtList[i]->j; */
  /* kk = gridPtList[i]->k; */
  /* dd2 = CPdd[ii][jj][kk]; */
  /*   //printf("(%d,%d,%d) = \t%g, diff = %g\n", ii,jj,kk, dd1, dd2-dd1); */
  /*   if (abs(dd1-dd2) >= maxdiff) { */
  /*     maxdiff = abs(dd1-dd2); */
  /*   } */
  /* } */
  /* printf("DEBUG max diff: %g\n", maxdiff); */

  /*dbg_printf(10, "Outputting to file...\n");*/
  /*outputGridToFile("Griddata", 1);*/
#ifdef ALSOTHREEDMATRICES
  /*outputMatrixToFile(CPx, 0, "CPdatax");
    outputMatrixToFile(CPy, 0, "CPdatay");
    outputMatrixToFile(CPz, 0, "CPdataz");
    outputMatrixToFile(CPdd, 0, "CPdatadd");*/
#endif

}



/*
 *
 */
void cleanup(void) {

#ifdef ALSOTHREEDMATRICES
  freeMatrix(CPx);
  freeMatrix(CPy);
  freeMatrix(CPz);
  freeMatrix(CPdd);
#endif

  freePlyData();

  dbg_printf(10, "destroying hashtable...\n");
  hashDeleteAll();
}





/*
 * The function matlab calls
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *indx, *inrelpt, *inbw, *inF, *inV, *inPrune, *inDebugLevel;
  double *mxIJK, *mxDD, *mxCP, *mxXYZ, *mxWHICHF;

  int i;
  const int numInputs = 7;
  struct_cp *s;
  unsigned long d;
  int withPruning;
  long bandsz;


  /* input checking */
  if (nrhs != numInputs) {
    mexErrMsgTxt("wrong number of arguments");
  }

  for (i=0; i < numInputs; i++) {
    if ( mxIsComplex(prhs[i]) || mxIsClass(prhs[i],"sparse") || mxIsChar(prhs[i]) ) {
      mexErrMsgTxt("Input must be real, full, and nonstring");
    }
  }

  /* TODO: maybe do some dimension checking here too */
  /*dims = mxGetDimensions(prhs[3]);
    numdims = mxGetNumberOfDimensions(prhs[3]);
    dimy = (int)dims[0]; dimx = (int)dims[1];
    mexPrintf("DEBUG: [%d,%d]\n", dimy, dimx);*/


  /* extract data from the inputs */
  indx = mxGetPr(prhs[0]);
  inrelpt = mxGetPr(prhs[1]);
  inbw = mxGetPr(prhs[2]);
  inF = mxGetPr(prhs[3]);
  inV = mxGetPr(prhs[4]);
  inPrune = mxGetPr(prhs[5]);
  inDebugLevel = mxGetPr(prhs[6]);

  number_faces = (unsigned long) mxGetM(prhs[3]);
  number_vertices = (unsigned long) mxGetM(prhs[4]);

  DEBUG_LEVEL = (int) inDebugLevel[0];
  dbg_printf(99, "setting debug level to %d: (display messages with level lower than %d)\n", DEBUG_LEVEL, DEBUG_LEVEL);

  DX = indx[0];
  dbg_printf(10, "DX=%g\n", DX);
  RELPTX = inrelpt[0];
  RELPTY = inrelpt[1];
  RELPTZ = inrelpt[2];
  dbg_printf(10, "RELPT=(%g,%g,%g)\n", RELPTX, RELPTY, RELPTZ);
  BANDWIDTH = inbw[0];
  dbg_printf(10, "BANDWIDTH=%g (%g*dx)\n", BANDWIDTH, BANDWIDTH/DX);

  if (inPrune[0] >= 0) {
    withPruning = 1;
  } else {
    withPruning = 0;
  }
  dbg_printf(10, "withPruning=%d\n", withPruning);

  dbg_printf(10, "number_faces=%ld, number_vertices=%ld\n", number_faces, number_vertices);



  /* todo: move to main()? */
  hashInit();

  initShapeFromMatlabArray(inF, number_faces, inV, number_vertices);

  mainRoutine();


  /* This first pass through the results it just to find how many
     there are.  TODO: a bad idea: we rely on processing these
     exactly the same twice. */
  d = 0;
  for(s=newhashlist; s != NULL; s=s->hh.next) {
    if ( (!withPruning) || (s->dd <= BANDWIDTH*BANDWIDTH) ) {
      d++;
    }
  }
  bandsz = d;
  dbg_printf(2, "after pruning, counted %ld gridpts\n", bandsz);

  /* create output arrays */
  plhs[0] = mxCreateDoubleMatrix(bandsz, 3, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(bandsz, 1, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(bandsz, 3, mxREAL);
  /* TODO: IJK and which_face should be int? */
  plhs[3] = mxCreateDoubleMatrix(bandsz, 3, mxREAL);
  plhs[4] = mxCreateDoubleMatrix(bandsz, 1, mxREAL);
  /* associate pointers */
  mxIJK = mxGetPr(plhs[0]);
  mxDD  = mxGetPr(plhs[1]);
  mxCP  = mxGetPr(plhs[2]);
  mxXYZ = mxGetPr(plhs[3]);
  mxWHICHF = mxGetPr(plhs[4]);


  d = 0;
  for(s = newhashlist; s != NULL; s = s->hh.next) {
    if ( (!withPruning) || (s->dd <= BANDWIDTH*BANDWIDTH) ) {
      mxIJK[0*bandsz+d] = s->i + 1;  /* plus 1 here for matlab indexing */
      mxIJK[1*bandsz+d] = s->j + 1;
      mxIJK[2*bandsz+d] = s->k + 1;
      mxDD[d] = s->dd;  /* squared distance */
      mxCP[0*bandsz+d] = s->cpx;
      mxCP[1*bandsz+d] = s->cpy;
      mxCP[2*bandsz+d] = s->cpz;
      mxXYZ[0*bandsz+d] = s->x;
      mxXYZ[1*bandsz+d] = s->y;
      mxXYZ[2*bandsz+d] = s->z;
      mxWHICHF[d] = s->which_face + 1;  /* +1 for matlab indexing */

      if ( (d < 5) || (d >= bandsz - 5) ) {
	dbg_printf(20,							\
		   "%5d: %d %d %d   %9.6g %9.6g %9.6g   %9.6f   %7.4g %7.4g %7.4g, face=%d\n", \
		   d, s->i,s->j,s->k, s->dd, s->cpx,s->cpy,s->cpz, s->x,s->y,s->z, s->which_face);
      }

      d++;
    }
  }
  dbg_printf(10, "Passing %ld gridpts back to Matlab\n", d);

  cleanup();

  return;
}
