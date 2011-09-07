/********************************************************************
 * HELPER_TRI2CP.c
 * Takes triangles (as faces and vertices) and returns closest points.
 * Written as a mex interface for matlab/octave.
 *
 * Authors: Colin Macdonald with routines by Steve Ruuth
 *
 * To compile, type this in matlab:
 * 
 * mex -O CFLAGS='\$CFLAGS -Wall -std\=c99' helper_tri2cp.c
 *
 * Maintainence: recompiling this on many different OSes is a hassle:
 * its better to remove intelligence from this and put it in the M
 * file instead.
 *
 * Limitations:
 *
 * If use ANSI C, could compile easier in matlab:
 *   "mex helper_tri2cp.c"
 *
 * command line also did extended precision: may not be completely
 * working yet.  Anyway, might not be a good idea with mex anyway.
 *
 * Might be some bug where different RELPTs take quite different
 * amounts of time: this shouldn't happen.  problem with hashtable,
 * try a different one?  Maybe the google one?
 * See also the funny snprintf(tupstr, TUPSTRSZ, "(%d,%d,%d)") results
 *
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
/* tgmath.h works for all precisions, but tgmath.h not easily
   available on cygwin. */
/*#include <tgmath.h>*/
#include <math.h>
#include <search.h>
#include <string.h>
#include <errno.h>
#include <mex.h>




/* What type of boundingSpheres: might as well use the optimal one */
/*#define boundingSpheres boundingSpheresCentroid*/
#define boundingSpheres boundingSpheresOptimal


/* Use extended precision or not */
/* TODO: values like 0.5 appear in code but may need L for extended,
   example 0.0, 1.0, 2.0, 4.0 */
/*#define EXTENDEDPRECISION*/
#ifdef EXTENDEDPRECISION
typedef long double myfloat;
#else
typedef double myfloat;
#endif




/* Also form the 3D matrices: uses a lot of RAM / disk space. */
/* TODO: probably doesn't work with the new RELPT code */
/*#define ALSOTHREEDMATRICES*/
/*#include "mat3d.h"*/

/* 23 is enough space for 3D, 5 digit entries in tuple */
#define TUPSTRSZ 23


/* Using this instead of macro: works with different FP. */
/*const myfloat CONST_PI = 3.14159265358979323846L;*/


/* use 1st one for calls that set errno, 2nd otherwise */
/*
#define myerr_nix(s) { perror((s)); exit(EXIT_FAILURE); }
#define myerr(s) { fprintf(stderr, "%s\n", (s)); exit(EXIT_FAILURE); }
*/
#define myerr_nix(s) { mexErrMsgTxt(s); exit(EXIT_FAILURE); }
#define myerr(s) { mexErrMsgTxt(s); exit(EXIT_FAILURE); }




struct struct_vertex {
  myfloat x, y, z;
  myfloat confidence;
  myfloat intensity;
};

struct struct_face {
  long v1, v2, v3;
};

/* Struct used for the hashtable.  Don't need the (x,y,z) but it might
   be useful for debugging.  Do need (i,j,k) b/c search.h hashtable
   doesn't seem to give a way to walk the table. */
struct struct_cp {
  myfloat dd;
  myfloat cpx, cpy, cpz;
  int i, j, k;
  myfloat x, y, z;
};




/* Global variables */
myfloat RELPTX, RELPTY, RELPTZ;
myfloat DX;
myfloat BANDWIDTH;
int DEBUG_LEVEL = 0;  /* messages with smaller level will display */
long number_vertices, number_faces;

/* TODO: probably no need to be global */
long ExpectedGridSize, HashTableSize;

struct struct_cp **gridPtList;
char **keyList;
long numgridpts = 0;

struct struct_vertex *vertex;
struct struct_face   *face;
myfloat *global_Sphere_Cx;
myfloat *global_Sphere_Cy;
myfloat *global_Sphere_Cz;
myfloat *global_Sphere_w;




#define BUFMAXSIZE 1000
void dbg_printf(int level, const char *fmt,...)
{
  char buf[BUFMAXSIZE];
  va_list ap;

  if (level <= DEBUG_LEVEL) {
    va_start(ap,fmt);
    vsnprintf(buf,BUFMAXSIZE,fmt,ap);
    va_end(ap);
    /*mexPrintf("  %s: %s", mexFunctionName(), buf);*/
    mexPrintf("  %s: %s", "DEBUG", buf);
  }
}



void mallocPlyData(long numv, long numf)
{
  if ((vertex = (struct struct_vertex *)			\
       malloc(numv * sizeof(struct struct_vertex))) == NULL)
    myerr_nix("Error mallocing ply data")

  if ((face = (struct struct_face *) malloc(numf * sizeof(struct struct_face))) == NULL)
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
 * read data from a ply file
 */
#define CONST_XSHIFT	0.0
#define CONST_YSHIFT	0.0
#define CONST_ZSHIFT	0.0
#define XSCALE    1.0
#define YSCALE    1.0
#define ZSCALE    1.0
void initShapeFromFile(const char* fname)
{
  long i;
  FILE *fp;
  long dummy;
  char dstr[256];
  char dstr2[256];
  char dstr3[256];

  if ((fp = fopen(fname, "r")) == NULL) myerr_nix("Error: can't open ply file");

  /* could auto detect this? */
  if (1==0) {
    /* nonstandard simply ply file */
    fscanf(fp, "%ld %ld",&number_vertices, &number_faces);
  } else {
    /* standard ply file (only one comment line allowed) */
    fgets(dstr, 255, fp);
    dbg_printf(5, "Debug: read input: %s", dstr);
    if (strcmp("ply\n", dstr) != 0) {
      dbg_printf(5, "not a ply file (or not one I can read!\n");
      dbg_printf(5, "input line was: %s", dstr);
    }
    while(1) {
      fgets(dstr, 255, fp);
      dbg_printf(5, "Debug: read input: %s", dstr);
      if (strncmp("element vertex", dstr, 14) == 0) {
	dbg_printf(5, "Debug: ... found the 'element vertex'\n");
	break;
      } else {
	//dbg_printf(5, "Debug: ... disgarding\n");
      }
    }
    sscanf(dstr, "%s %s %ld\n", dstr2, dstr3, &number_vertices);
    dbg_printf(5, "Debug: vertices: %ld\n", number_vertices);

    while(1) {
      fgets(dstr, 255, fp);
      dbg_printf(5, "Debug: reading input: %s", dstr);
      if (strncmp("element face", dstr, 12) == 0) {
	dbg_printf(5, "Debug: ... found the 'element face'\n");
	break;
      } else {
	//dbg_printf(5, "Debug: ... disgarding\n");
      }
    }
    sscanf(dstr, "%s %s %ld\n", dstr2, dstr3, &number_faces);
    dbg_printf(5, "Debug: faces: %ld\n", number_faces);

    //fgets(dstr, 255, fp);
    //fgets(dstr, 255, fp);
    //fscanf(fp, "%s %s %ld\n", dstr, dstr2, &number_vertices);
    //fgets(dstr, 255, fp);
    //fgets(dstr, 255, fp);
    //fgets(dstr, 255, fp);
    //fscanf(fp, "%s %s %ld\n", dstr, dstr2, &number_faces);

    /* TODO: hardcoded for two more lines:
       property list uchar int vertex_indices
       end_header
    */
    fgets(dstr, 255, fp);
    fgets(dstr, 255, fp);
  }

  dbg_printf(5, "Debug: vertices, faces: %ld %ld\n", number_vertices, number_faces);

  mallocPlyData(number_vertices, number_faces);

  for (i=0; i<number_vertices; i++) {
#ifdef EXTENDEDPRECISION
    fscanf(fp, "%Lf %Lf %Lf", &vertex[i].x, &vertex[i].y, &vertex[i].z);
#else
    fscanf(fp, "%lf %lf %lf", &vertex[i].x, &vertex[i].y, &vertex[i].z);
#endif
    vertex[i].x = (vertex[i].x - CONST_XSHIFT) / XSCALE;
    vertex[i].y = (vertex[i].y - CONST_YSHIFT) / YSCALE;
    vertex[i].z = (vertex[i].z - CONST_ZSHIFT) / ZSCALE;
    //if (vertex[i].z < -1.0 | vertex[i].z > 1.0) {
    //  myerr("HELP!");
    //}
    /*
      fscanf(fp,"%lf %lf %lf %lf %lf", &vertex[i].x, &vertex[i].y, &vertex[i].z, &vertex[i].confidence, &vertex[i].intensity);
      fscanf(fp,"%lf %lf %lf", &vertex[i].x, &vertex[i].y, &vertex[i].z);
      printf("%g %g %g %d %d %d\n", vertex[i].x, vertex[i].y, vertex[i].z, (vertex[i].x>0), 0, (vertex[i].x<0));
    */
  }
  for (i=0; i<number_faces; i++) {
    fscanf(fp,"%ld %ld %ld %ld",&dummy, &face[i].v1, &face[i].v2, &face[i].v3);
    //printf("%ld %ld %ld %ld\n",  dummy,  face[i].v1,  face[i].v2,  face[i].v3);
  }
}



void initShapeFromMatlabArray(double *Faces, long numF, double *Vertices, long numV)
{
  long i;

  mallocPlyData(numV, numF);

  //mexPrintf("%ld %ld %ld %ld\n", i, vertex[i].x, vertex[i].y, vertex[i].z);
  for (i=0; i<numV; i++) {
    vertex[i].x = Vertices[0*numV+i];
    vertex[i].y = Vertices[1*numV+i];
    vertex[i].z = Vertices[2*numV+i];
    //mexPrintf("%ld %g %g %g\n", i, vertex[i].x, vertex[i].y, vertex[i].z);
    //if (vertex[i].z < -1.0 | vertex[i].z > 1.0) {
    //  printf("HELP!\n");
    //  error(-1);
    //}
  }
  for (i=0; i<number_faces; i++) {
    /* -1 here because Faces is from Matlab */
    face[i].v1 = Faces[0*numF+i] - 1;
    face[i].v2 = Faces[1*numF+i] - 1;
    face[i].v3 = Faces[2*numF+i] - 1;
    //mexPrintf("%ld %ld %ld %ld\n", i, face[i].v1, face[i].v2, face[i].v3);
  }
}



/*
 * int away from zero: like ceil, floor but always moves away from
 * zero.  Special case: iafz(+-0.0) = 0
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
  /* wrong b/c of signed zero: */
  //return( (int) (x + signbit(x)) ); */

  /* this approach cannot be trusted for large x */
  //int signum;
  //signum = (x > 0) - (x < 0);
  //return( (int) (x + signum) );
}



inline myfloat double_max(myfloat a, myfloat b)
{
  if (a>=b) return(a);
  return(b);
}


inline myfloat double_min(myfloat a, myfloat b)
{
  if (a<=b) return(a);
  return(b);
}



myfloat distance2(myfloat x, myfloat y, myfloat z, long j)
{
	myfloat dd;

	dd = (x-global_Sphere_Cx[j])*(x-global_Sphere_Cx[j])
	   + (y-global_Sphere_Cy[j])*(y-global_Sphere_Cy[j])
	   + (z-global_Sphere_Cz[j])*(z-global_Sphere_Cz[j]);
	return(dd);
}

myfloat distance(long i, long j)
{
	myfloat dd;

	dd = (vertex[i].x-global_Sphere_Cx[j])*(vertex[i].x-global_Sphere_Cx[j])
	   + (vertex[i].y-global_Sphere_Cy[j])*(vertex[i].y-global_Sphere_Cy[j])
	   + (vertex[i].z-global_Sphere_Cz[j])*(vertex[i].z-global_Sphere_Cz[j]);
	return(sqrt(dd));
}




void boundingSpheresOptimal()
{
  /* Wikipedia: The useful minimum bounding circle of three points is
     defined either by the circumcircle (where three points are on the
     minimum bounding circle) or by the two points of the longest side
     of the triangle (where the two points define a diameter of the
     circle). */
  long i;
  myfloat ang[3];
  myfloat m[3];
  myfloat a[3], b[3], c[3];
  // TODO: change to arrays
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

    //anga = (C-A)'*(B-A);
    //angb = (A-B)'*(C-B);
    //angc = (A-C)'*(B-C);
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
      //dbg_printf(99, "obtuse triangle, use longest side\n");
      if (minangi == 0) {
	//m = (B+C)/2.0;
	//TODO: better?: m = (B-C)/2.0 + C
	m[0] = (b[0] + c[0]) / 2.0;
	m[1] = (b[1] + c[1]) / 2.0;
	m[2] = (b[2] + c[2]) / 2.0;
	R = ( (b[0] - c[0]) * (b[0] - c[0]) +
	      (b[1] - c[1]) * (b[1] - c[1]) +
	      (b[2] - c[2]) * (b[2] - c[2]) ) / 4.0;
      } else if (minangi == 1) {
	//m = (A+C)/2.0;
	m[0] = (a[0] + c[0]) / 2.0;
	m[1] = (a[1] + c[1]) / 2.0;
	m[2] = (a[2] + c[2]) / 2.0;
	R = ( xca*xca + yca*yca + zca*zca ) / 4.0;
      } else if (minangi == 2) {
	//m = (A+B)/2.0;
	m[0] = (a[0] + b[0]) / 2.0;
	m[1] = (a[1] + b[1]) / 2.0;
	m[2] = (a[2] + b[2]) / 2.0;
	R = ( xba*xba + yba*yba + zba*zba ) / 4.0;
      } else {
	myerr("no case");
      }
      R = sqrt(R);
    } else {
      //dbg_printf(99, "acute triangle, use circumcenter\n");
      /* uses code from
	 http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html */

      /*
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
      /* TODO: 0.5 here is double, need to fix for extended */
      denominator = 0.5 / (xcrossbc * xcrossbc + ycrossbc * ycrossbc +
			   zcrossbc * zcrossbc);

      /* Calculate offset (from `a') of circumcenter. */
      m[0] = ((balength * yca - calength * yba) * zcrossbc -
	      (balength * zca - calength * zba) * ycrossbc) * denominator;
      m[1] = ((balength * zca - calength * zba) * xcrossbc -
	      (balength * xca - calength * xba) * zcrossbc) * denominator;
      m[2] = ((balength * xca - calength * xba) * ycrossbc -
	      (balength * yca - calength * yba) * xcrossbc) * denominator;
      //R = sqrt(xcirca*xcirca + ycirca*ycirca + zcirca*zcirca);
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
    //printf("A=[%.20g,%.20g,%.20g]'\n", a[0],a[1],a[2]);
    //printf("B=[%.20g,%.20g,%.20g]'\n", b[0],b[1],b[2]);
    //printf("C=[%.20g,%.20g,%.20g]'\n", c[0],c[1],c[2]);
    //printf("Cen=[%.20g,%.20g,%.20g]'\n", global_Sphere_Cx[i], global_Sphere_Cy[i], global_Sphere_Cz[i]);
    //printf("rad=%.20g\n", global_Sphere_w[i]);
  } /* end loop over triangles */
}



/* compute bounding spheres using the centroid, easy algorithm but not
   very tight. */
void boundingSpheresCentroid()
{
  long i;
  myfloat R;

  for (i=0; i<number_faces; i++) {
    global_Sphere_Cx[i] = (vertex[face[i].v1].x+vertex[face[i].v2].x+vertex[face[i].v3].x)/3.0;
    global_Sphere_Cy[i] = (vertex[face[i].v1].y+vertex[face[i].v2].y+vertex[face[i].v3].y)/3.0;
    global_Sphere_Cz[i] = (vertex[face[i].v1].z+vertex[face[i].v2].z+vertex[face[i].v3].z)/3.0;
    R = double_max(distance(face[i].v1,i),double_max(distance(face[i].v2,i),distance(face[i].v3,i)));
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
  lambda_star = double_max(0.0,double_min(lambda,1.0));

  *c1 = p1+lambda_star*qmp1;
  *c2 = p2+lambda_star*qmp2;
  *c3 = p3+lambda_star*qmp3;
}


/* Closest point and distance from a point (a1,a2,a3) to a triangle
 * indexed by face_index.  Returns the distance and the closest point
 * is in (c1,c2,c3).  Uses global vars `face' and `vertex'.
 * (code by Steve Ruuth)
 */
myfloat FindClosestPointToOneTri(myfloat a1, myfloat a2, myfloat a3, \
				 long face_index, \
				 myfloat *c1, myfloat *c2, myfloat *c3)
{
  long index_p, index_q, index_r;
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
  factor = 1.0/(a11*a22-a12*a12);
  i11 = a22*factor;
  i12 =-a12*factor;
  i22 = a11*factor;
  lambda = i11*b1+i12*b2;
  mu     = i12*b1+i22*b2;
  *c1 = lambda*q1+mu*r1;
  *c2 = lambda*q2+mu*r2;
  *c3 = lambda*q3+mu*r3;

  if ((lambda<0) && (mu<0) && (lambda+mu<=1)) {
    *c1 = *c2 = *c3 = 0.0;
  } else if ((lambda>=0) && (mu<0) && (lambda+mu<=1)) {
    ProjectOnSegment(c1,c2,c3,0.0,0.0,0.0,q1,q2,q3);
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
    ProjectOnSegment(c1,c2,c3,r1,r2,r3,0.0,0.0,0.0);
  } else if ((lambda>=0) && (mu>=0) && (lambda+mu<=1)) {
    /* do nothing */
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
  // HORRIBLE HACK 2010-07-28
  if (isinf(factor)) {
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
  return (dd);
}



/* a real global search */
myfloat FindClosestPointGlobally2(myfloat a1, myfloat a2, myfloat a3, \
				  myfloat *c1, myfloat *c2, myfloat *c3)
{
  long i;
  myfloat dd, dd_min;
  myfloat t1, t2, t3;

  dd_min = FindClosestPointToOneTri(a1,a2,a3, 0, &t1, &t2, &t3);
  *c1 = t1;
  *c2 = t2;
  *c3 = t3;
  for (i=1; i<number_faces; i++) {
    dd =  FindClosestPointToOneTri(a1,a2,a3, i, &t1, &t2, &t3);
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
void FindClosestPointGlobally()
{
  long n;
  int i,j,k;
  int iL,iU;
  int jL,jU;
  int kL,kU;
  myfloat ih;
  myfloat dd;
  myfloat x,y,z;
  myfloat cpx, cpy, cpz;
  long cc = 0;
  ENTRY e, *ep;
  struct struct_cp *gridpt;
  char *tupstr, *tupstr2;

  /* colin.key = strdup("(10,-20,300)"); */
  /* printf("key: %s\n", colin.key); */
  /* ep = hsearch(colin, FIND); */
  /* printf("%s -> %s:%d\n", colin.key, */
  /* 	 ep ? ep->key : "NULL", ep ? (int)(ep->data) : 0); */


  if ((tupstr = (char *) malloc(TUPSTRSZ)) == NULL)
    myerr_nix("Error mallocing tuple string");


  /* introduce dummy distance to the surface.  Outside band means it mgiht as well be infinite */
#ifdef ALSOTHREEDMATRICES
  for (i=0; i<NPOINTS; i++)
    for (j=0; j<NPOINTS; j++)
      for (k=0; k<NPOINTS; k++)
	CPdd[i][j][k] = 987654321.;
#endif

  long total_count = 0;

  /* 1/h or inverse of h */
  //ih = (NPOINTS-1.)/(DOMAIN_B-DOMAIN_A);
  ih = 1/DX;
  for (n=0; n<number_faces; n++) {
    //i0 =  ceil(((global_Sphere_Cx[n]-global_Sphere_w[n])-DOMAIN_A)*ih);
    //in =  (((global_Sphere_Cx[n]+global_Sphere_w[n])-DOMAIN_A)*ih);
    //if (i0<0) i0=0;
    //if (in>=NPOINTS) in=NPOINTS-1;
    //j0 =  ceil(((global_Sphere_Cy[n]-global_Sphere_w[n])-DOMAIN_A)*ih);
    //jn =  (((global_Sphere_Cy[n]+global_Sphere_w[n])-DOMAIN_A)*ih);
    //if (j0<0) j0=0;
    //if (jn>=NPOINTS) jn=NPOINTS-1;
    //k0 =  ceil(((global_Sphere_Cz[n]-global_Sphere_w[n])-DOMAIN_A)*ih);
    //kn =  (((global_Sphere_Cz[n]+global_Sphere_w[n])-DOMAIN_A)*ih);
    //if (k0<0) k0=0;
    //if (kn>=NPOINTS) kn=NPOINTS-1;
    /* arbitrary ref point support: */
    iL = iafz( (global_Sphere_Cx[n]-global_Sphere_w[n] - RELPTX) * ih );
    iU = iafz( (global_Sphere_Cx[n]+global_Sphere_w[n] - RELPTX) * ih );
    jL = iafz( (global_Sphere_Cy[n]-global_Sphere_w[n] - RELPTY) * ih );
    jU = iafz( (global_Sphere_Cy[n]+global_Sphere_w[n] - RELPTY) * ih );
    kL = iafz( (global_Sphere_Cz[n]-global_Sphere_w[n] - RELPTZ) * ih );
    kU = iafz( (global_Sphere_Cz[n]+global_Sphere_w[n] - RELPTZ) * ih );

    for (i=iL; i<=iU; i++) {
      //x = ((myfloat)i)*(DOMAIN_B-DOMAIN_A)/(NPOINTS-1)+DOMAIN_A;
      x = ((myfloat)i)*DX + RELPTX;
      for (j=jL; j<=jU; j++) {
	//y = ((myfloat)j)*(DOMAIN_B-DOMAIN_A)/(NPOINTS-1)+DOMAIN_A;
	y = ((myfloat)j)*DX + RELPTY;
	for (k=kL; k<=kU; k++) {
	  //z = ((myfloat)k)*(DOMAIN_B-DOMAIN_A)/(NPOINTS-1)+DOMAIN_A;
	  z = ((myfloat)k)*DX + RELPTZ;
	  total_count++;
	  if ( ( (x-global_Sphere_Cx[n])*(x-global_Sphere_Cx[n]) +
		 (y-global_Sphere_Cy[n])*(y-global_Sphere_Cy[n]) +
		 (z-global_Sphere_Cz[n])*(z-global_Sphere_Cz[n]) )
	       <= global_Sphere_w[n]*global_Sphere_w[n] ) {
	    dd =  FindClosestPointToOneTri(x,y,z, n, &cpx, &cpy, &cpz);
	    cc++;
	    //printf("face: %d, point: (%d,%d,%d), x: (%g,%g,%g)\n", n, i,j,k, x,y,z);
	    //printf("%d %d %d %d %g %g %g %g\n", n, i,j,k, dd, cpx,cpy,cpz);

#ifdef ALSOTHREEDMATRICES
	    if (dd<CPdd[i][j][k]) {
	      CPdd[i][j][k] = dd;
	      CPx[i][j][k]  = cpx;
	      CPy[i][j][k]  = cpy;
	      CPz[i][j][k]  = cpz;
	    }
#endif

	    // TODO: strangely, this  strongly effects the resulting time
	    //if ( snprintf(tupstr, TUPSTRSZ, "(%d,%d,%d)", i+640,j+640,k+640) >= TUPSTRSZ )
	    if ( snprintf(tupstr, TUPSTRSZ, "(%d,%d,%d)", i,j,k) >= TUPSTRSZ )
	      myerr("Couldn't write tuple string, increase TUPSTRSZ?");
	    //dbg_printf(5, "DEBUG: tupstr=\"%s\"\n", tupstr);

	    e.key = tupstr;
	    e.data = NULL;  // doesn't matter, only need key

	    ep = hsearch(e, FIND);
	    if (ep == NULL) {
	      /* we didn't find it so enter a new one */
	      if ((gridpt = (struct struct_cp *)		\
		   malloc(sizeof(struct struct_cp))) == NULL)
		myerr_nix("Error allocating gridpt");

	      /* we need a duplicate string to put in the table b/c
		 tupstr gets reused next time through the loop. */
	      if ((tupstr2 = strdup(tupstr)) == NULL)
		myerr_nix("Error duplicating tuple string");

	      gridpt->dd = dd;
	      gridpt->cpx = cpx;
	      gridpt->cpy = cpy;
	      gridpt->cpz = cpz;
	      gridpt->i = i;
	      gridpt->j = j;
	      gridpt->k = k;
	      gridpt->x = x;
	      gridpt->y = y;
	      gridpt->z = z;
	      //dbg_printf(5, "DEBUG: gridPtCtr=%d\n", gridPtCtr);
	      gridPtList[numgridpts] = gridpt;
	      keyList[numgridpts] = tupstr2;
	      numgridpts++;
	      e.key = tupstr2;
	      e.data = (void *) gridpt;
	      /* add it e to the hashtable */
	      ep = hsearch(e, ENTER);
	      if (ep == NULL) myerr_nix("Error, maybe hashtable is full\? Msg:");
	    } else {
	      /* it was already in the table, check if new pt is
		 closer */
	      gridpt = (struct struct_cp*) ep->data;
	      if (dd < (gridpt->dd)) {
		gridpt->dd = dd;
		gridpt->cpx = cpx;
		gridpt->cpy = cpy;
		gridpt->cpz = cpz;
		gridpt->i = i;
		gridpt->j = j;
		gridpt->k = k;
		gridpt->x = x;
		gridpt->y = y;
		gridpt->z = z;
	      }
	    }
	  } // end inside sphere
	}
      }
    }
  }
  dbg_printf(10, "file writer: total loops: %ld\n", total_count);
  dbg_printf(10, "file writer: made %ld calls to Find_CP_to_tri()\n", cc);
  dbg_printf(10, "file writer: found %ld grid points\n", numgridpts);

  free(tupstr);
}






/*
 * output a grid to a file
 */
void outputGridToFile(const char *fname, int withPruning)
{
  long c, d = 0;
  FILE *fd;

  if ((fd = fopen(fname, "w")) == NULL) myerr_nix("Error writing griddata");

  for (c = 0; c < numgridpts; c++) {
    if ((!withPruning) || ((gridPtList[c]->dd <= BANDWIDTH*BANDWIDTH))) {
      // ! (1 and 1 = 1) = 0
      // ! (1 and 0 = 0) = 1
      // ! (0 and d = 0) = 1
      d++;
#ifdef EXTENDEDPRECISION
#define PSTR "%d %d %d %.21Le %.21Le %.21Le %.21Le %.21Le %.21Le %.21Le\n"
#else
#define PSTR "%d %d %d %.18g %.18g %.18g %.18g %.18g %.18g %.18g\n"
#endif
      fprintf(fd, PSTR,						\
	      gridPtList[c]->i,					\
	      gridPtList[c]->j,					\
	      gridPtList[c]->k,					\
	      sqrt(gridPtList[c]->dd),				\
	      gridPtList[c]->cpx,				\
	      gridPtList[c]->cpy,				\
	      gridPtList[c]->cpz,				\
	      gridPtList[c]->x,					\
	      gridPtList[c]->y,					\
	      gridPtList[c]->z);
    }
  }
  if ( fclose(fd) != 0 )  myerr_nix("Error closing file");
  dbg_printf(10, "Saved %ld gridpts to file\n", d);
}






/*
 * this was main() before adding mex
 */
void findGrid(void)
{
  long i;
  long expectedgridsz, hashtablesz;

  /*
  printf("%g, %d\n", -0.0, iafz(-0.0));
  printf("%g, %d\n", 0.0, iafz(0.0));
  printf("%g, %d\n", -0.1, iafz(-0.1));
  printf("%g, %d\n", 0.1, iafz(0.1));
  printf("%g, %d\n", -1.1, iafz(-1.1));
  printf("%g, %d\n", 1.1, iafz(1.1));
  */

#ifdef ALSOTHREEDMATRICES
  dbg_printf(10, "Debug: mallocing 3D matrices...\n");
  CPx = mallocMatrix();
  CPy = mallocMatrix();
  CPz = mallocMatrix();
  CPdd = mallocMatrix();
#endif

  expectedgridsz = ExpectedGridSize;
  hashtablesz = HashTableSize;
  dbg_printf(2, "expected grid size: %ld\n", expectedgridsz);
  dbg_printf(2, "max hash tabel size: %ld\n", hashtablesz);


  dbg_printf(12, "Debug: mallocing lists... \n");
  gridPtList = (struct struct_cp **) malloc(expectedgridsz*sizeof(struct struct_cp *));
  if (gridPtList == NULL) myerr_nix("Error: malloc() grid list failed, msg");
  for (i = 0; i < expectedgridsz; i++)
    gridPtList[i] = (struct struct_cp *) NULL;

  if ((keyList = (char **) malloc(expectedgridsz*sizeof(char *))) == NULL)
    myerr_nix("Error: malloc() key list failed, msg");
  for (i = 0; i < expectedgridsz; i++)
    keyList[i] = (char *) NULL;

  dbg_printf(4, "done mallocing\n");

  /* Create a new hash table.  To be portable we only can make one of
     these.  On GNU, could use hcreate_r(). */
  if (hcreate(hashtablesz) == 0) myerr_nix("Error creating hash table");


  dbg_printf(10, "running boundingSpheres()...\n");
  boundingSpheres();
  dbg_printf(10, "Finding Grid...\n");
  FindClosestPointGlobally();

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

  //dbg_printf(10, "Outputting to file...\n");
  //outputGridToFile("Griddata", 1);
#ifdef ALSOTHREEDMATRICES
  outputMatrixToFile(CPx, 0, "CPdatax");
  outputMatrixToFile(CPy, 0, "CPdatay");
  outputMatrixToFile(CPz, 0, "CPdataz");
  outputMatrixToFile(CPdd, 0, "CPdatadd");
#endif

}



void cleanup(void) {
  long i;

#ifdef ALSOTHREEDMATRICES
  freeMatrix(CPx);
  freeMatrix(CPy);
  freeMatrix(CPz);
  freeMatrix(CPdd);
#endif

  freePlyData();
  hdestroy();

  for (i = 0; i < numgridpts; i++) {
    free(gridPtList[i]);
    free(keyList[i]);
  }
  free(gridPtList);
  free(keyList);
  //return(0);
}





/*
 * The function matlab calls
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  //mxArray *a_in_m, *b_in_m, *c_in_m, *d_in_m, *e_in_m;
  double *indx, *inrelpt, *inbw, *inF, *inV, *inMaxHashSz, *inDebugLevel;
  //const mwSize *dims;
  double *mxIJK, *mxDD, *mxCP, *mxXYZ;
  //int dimx, dimy, numdims;
  int i;
  const int numInputs = 7;

  if (nrhs != numInputs) {
    mexErrMsgTxt("wrong number of arguments");
  }

  for (i=0; i < numInputs; i++) {
    if ( mxIsComplex(prhs[i]) || mxIsClass(prhs[i],"sparse") || mxIsChar(prhs[i]) ) {
      mexErrMsgTxt("Input must be real, full, and nonstring");
    }
  }
  //associate inputs
  /*
  a_in_m = mxDuplicateArray(prhs[0]);
  b_in_m = mxDuplicateArray(prhs[1]);
  c_in_m = mxDuplicateArray(prhs[2]);
  d_in_m = mxDuplicateArray(prhs[3]);
  e_in_m = mxDuplicateArray(prhs[4]);

  in0 = mxGetPr(a_in_m);
  in1 = mxGetPr(b_in_m);
  in2 = mxGetPr(c_in_m);
  inF = mxGetPr(d_in_m);
  inV = mxGetPr(e_in_m);
  */

  indx = mxGetPr(prhs[0]);
  inrelpt = mxGetPr(prhs[1]);
  inbw = mxGetPr(prhs[2]);
  inF = mxGetPr(prhs[3]);
  inV = mxGetPr(prhs[4]);
  inMaxHashSz = mxGetPr(prhs[5]);
  inDebugLevel = mxGetPr(prhs[6]);

  number_faces = mxGetM(prhs[3]);
  number_vertices = mxGetM(prhs[4]);

  DEBUG_LEVEL = (int) inDebugLevel[0];
  dbg_printf(99, "setting debug level to %d: (display messages with level lower than %d)\n", DEBUG_LEVEL, DEBUG_LEVEL);

  // TODO: error checking here on size of this array
  ExpectedGridSize = (int) inMaxHashSz[0];
  HashTableSize = (int) inMaxHashSz[1];
  dbg_printf(2, "expected number of gridpoints: %ld\n", ExpectedGridSize);
  dbg_printf(2, "hashtable size: %ld\n", HashTableSize);

  DX = indx[0];
  dbg_printf(10, "DX=%g\n", DX);
  RELPTX = inrelpt[0];
  RELPTY = inrelpt[1];
  RELPTZ = inrelpt[2];
  dbg_printf(10, "RELPT=(%g,%g,%g)\n", RELPTX, RELPTY, RELPTZ);
  BANDWIDTH = inbw[0];
  dbg_printf(10, "BANDWIDTH=%g (%g*dx)\n", BANDWIDTH, BANDWIDTH/DX);


  dbg_printf(10, "number_faces=%ld, number_vertices=%ld\n", number_faces, number_vertices);
  // TODO: extra trim boolean

  //figure out dimensions
  //dims = mxGetDimensions(prhs[3]);
  //numdims = mxGetNumberOfDimensions(prhs[3]);
  //dimy = (int)dims[0]; dimx = (int)dims[1];
  //mexPrintf("DEBUG: [%d,%d]\n", dimy, dimx);

  initShapeFromMatlabArray(inF, number_faces, inV, number_vertices);
  //initShapeFromFile("annies_pig.ply");

  findGrid();

  long c, d;
  int withPruning = 1;
  long bandsz;

  // TODO: bad idea here: we rely on processing these the same twice
  d = 0;
  for (c = 0; c < numgridpts; c++) {
    if ((!withPruning) || ((gridPtList[c]->dd <= BANDWIDTH*BANDWIDTH))) {
      // ! (1 and 1 = 1) = 0
      // ! (1 and 0 = 0) = 1
      // ! (0 and d = 0) = 1
      d++;
    }
  }
  bandsz = d;
  dbg_printf(2, "after pruning, counted %ld gridpts\n", bandsz);

  //dimx = 2;
  //dimy = 4;

  //associate outputs
  plhs[0] = mxCreateDoubleMatrix(bandsz, 3, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(bandsz, 1, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(bandsz, 3, mxREAL);
  plhs[3] = mxCreateDoubleMatrix(bandsz, 3, mxREAL);
  //associate pointers
  mxIJK = mxGetPr(plhs[0]);
  mxDD  = mxGetPr(plhs[1]);
  mxCP  = mxGetPr(plhs[2]);
  mxXYZ = mxGetPr(plhs[3]);

  //d_out_m = plhs[1] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);
  //a = mxGetPr(a_in_m);
  //b = mxGetPr(b_in_m);
  //C = mxGetPr(c_out_m);
  //D = mxGetPr(d_out_m);

  d = 0;
  for (c = 0; c < numgridpts; c++) {
    if ((!withPruning) || ((gridPtList[c]->dd <= BANDWIDTH*BANDWIDTH))) {
      if ( (d < 5) || (d >= bandsz - 5) ) {
	dbg_printf(20, "%5d: %d %d %d   %9.6g %9.6g %9.6g   %9.6f   %7.4g %7.4g %7.4g\n", \
		  d,							\
		  gridPtList[c]->i,					\
		  gridPtList[c]->j,					\
		  gridPtList[c]->k,					\
		  sqrt(gridPtList[c]->dd),				\
		  gridPtList[c]->cpx,					\
		  gridPtList[c]->cpy,					\
		  gridPtList[c]->cpz,					\
		  gridPtList[c]->x,					\
		  gridPtList[c]->y,					\
		  gridPtList[c]->z);
      }
      // TODO: I've added one here for option base 1 in matlab
      mxIJK[0*bandsz+d] = gridPtList[c]->i+1;
      mxIJK[1*bandsz+d] = gridPtList[c]->j+1;
      mxIJK[2*bandsz+d] = gridPtList[c]->k+1;
      // TODO: make a decision about this sqrt!
      mxDD[d] = sqrt(gridPtList[c]->dd);
      mxCP[0*bandsz+d] = gridPtList[c]->cpx;
      mxCP[1*bandsz+d] = gridPtList[c]->cpy;
      mxCP[2*bandsz+d] = gridPtList[c]->cpz;
      mxXYZ[0*bandsz+d] = gridPtList[c]->x;
      mxXYZ[1*bandsz+d] = gridPtList[c]->y;
      mxXYZ[2*bandsz+d] = gridPtList[c]->z;
      d++;
    }
  }
  dbg_printf(10, "Passing %ld gridpts back to Matlab\n", d);

  cleanup();

  return;
}
