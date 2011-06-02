/*********************************************************************
README:

** Compiling:
  gcc -lm cp_tri_hash.c -o cp_tri_hash
then run "./cp_tri_hash"

** Changing input file: search for "annies_pig", then recompile

** relevant things to change:

   NPOINTS: dx = 4/NPOINTS

   RELPTX,Y,Z: grid points are indexed relative to this point

** Output: a file named "Griddata"

  each row contains a grid point

  see function "outputGridToFile", to see what each column is (TODO:
  should document it better)

** Limitations:

  many!  input .ply file must have a very specific structure (must be
  exactly the same as annies_pig.ply

  Hardcoded for cubic interpolation (should be easy to fix though).

  Assumes model is roughly centered around origin (just to estimate
  size of hashtable)

  extended precision: may not be completely working yet.

  Might be some bug where different RELPTs take quite different
  amounts of time: this shouldn't happen.

  TODO: could use a sparse matrix class instead of a hashtable?
  Essentially Martin's idea.

*********************************************************************/

/* tgmath.h not easily available on cygwin, can use math.h */
#include <tgmath.h>
//#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <search.h>
#include <string.h>
#include <errno.h>


/************************************************************************
 *  GLOBAL CONSTANTS
 ***********************************************************************/
#define NPOINTS 161
/* Use 4/N here for historial reasons (used to use [-2,2]^3) */



/* Use extended precision or not */
/* TODO: values like 0.5 appear in code but may need L for extended,
   example 0.0, 1.0, 2.0, 4.0 */
//#define EXTENDEDPRECISION
#ifdef EXTENDEDPRECISION
typedef long double myfloat;
//#define myfloat long double
//#define DOMAIN_B ( 2.0L)
#define RELPTX (-2.0L)
#define RELPTY RELPTX
#define RELPTZ RELPTX
#define DX (4.0L / (NPOINTS-1))
#else
typedef double myfloat;
//#define myfloat double
//#define DOMAIN_B ( 2.0)
#define RELPTX (-2.0)
#define RELPTY RELPTX
#define RELPTZ RELPTX
#define DX (4.0 / (NPOINTS-1))
#endif
//#define DOMAIN_A	(-DOMAIN_B)


/* Also form the 3D matrices: uses a lot of RAM / disk space. */
/* TODO: probably doesn't work with the new RELPT code */
//#define ALSOTHREEDMATRICES

/* We estimate how many grid points we might need based on a sphere,
   then multiply by this much.  Maybe increase if your surface is
   really convoluted! */
#define EXPECTEDGRIDSZ_SAFETY_FACTOR 6


/* Degree of polynomial interpolation.  TODO: not used, needs a
   cleanup, don't hardcode BW below, etc. */
#define DEGREEP 3
#define POLYN DEGREEP
/* TODO: sort out bandwidth, maybe pass in as parameter */
//#define POINTS_DIRECTION	((POLYN+1.)/2.)
//#define BANDWIDTH (((DOMAIN_B-DOMAIN_A)*sqrt((POINTS_DIRECTION+1)*(POINTS_DIRECTION+1)+2*(POINTS_DIRECTION)*(POINTS_DIRECTION))+0.001)/(NPOINTS-1))
//#define DX ((DOMAIN_B-DOMAIN_A)/(NPOINTS-1))
#define BANDWIDTH (4.1235*DX)
//#define BANDWIDTH (4.1235179361802218878L*DX)


/* 23 is enough space for 3D, 5 digit entries in tuple */
#define TUPSTRSZ 23


/* Using this instead of macro: works with different FP. */
const myfloat CONST_PI = 3.14159265358979323846L;
//#define CONST_PI ((myfloat) 3.14159265358979323846L)

/* Where to direct debugging messages. */
#define dbout stderr

/* use 1st one for calls that set errno, 2nd otherwise */
#define myerr_nix(s) { perror((s)); exit(EXIT_FAILURE); }
#define myerr(s) { fprintf(stderr, "%s\n", (s)); exit(EXIT_FAILURE); }


#define matrix myfloat ***


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
matrix CPx;
matrix CPy;
matrix CPz;
matrix CPdd;

long number_vertices, number_faces;
struct struct_vertex *vertex;
struct struct_face   *face;
myfloat *global_Sphere_Cx;
myfloat *global_Sphere_Cy;
myfloat *global_Sphere_Cz;
myfloat *global_Sphere_w;

struct struct_cp **gridPtList;
char **keyList;
long numgridpts = 0;




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
// hand
#define CONST_XSHIFT	0.0
#define CONST_YSHIFT	0.0
#define CONST_ZSHIFT	-850.46
#define XSCALE    100.0
#define YSCALE    100.0
#define ZSCALE    100.0
*/

/*
// bumpy_sphere
#define CONST_XSHIFT	0.0
#define CONST_YSHIFT	0.0
#define CONST_ZSHIFT	0.0
#define XSCALE    10.0
#define YSCALE    10.0
#define ZSCALE    10.0
*/

// no changes to ply data
#define CONST_XSHIFT	0.0
#define CONST_YSHIFT	0.0
#define CONST_ZSHIFT	0.0
#define XSCALE    1.0
#define YSCALE    1.0
#define ZSCALE    1.0



/*
 * read data from a ply file
 */
void initShape(const char* fname)
{
  long i;
  FILE *fp;
  long dummy;
  char dstr[256];
  char dstr2[256];
  char dstr3[256];

  if ((fp = fopen(fname, "r")) == NULL) myerr_nix("Error: can't open ply file");

  /* TODO: auto detect this? */
  if (1==0) {
    /* nonstandard simply ply file */
    fscanf(fp, "%ld %ld",&number_vertices, &number_faces);
  } else {
    /* standard ply file (only one comment line allowed) */
    fgets(dstr, 255, fp);
    fprintf(dbout, "Debug: read input: %s", dstr);
    if (strcmp("ply\n", dstr) != 0) {
      fprintf(stderr, "not a ply file (or not one I can read!\n");
      fprintf(stderr, "input line was: %s", dstr);
    }
    while(1) {
      fgets(dstr, 255, fp);
      fprintf(dbout, "Debug: read input: %s", dstr);
      if (strncmp("element vertex", dstr, 14) == 0) {
	fprintf(dbout, "Debug: ... found the 'element vertex'\n");
	break;
      } else {
	//fprintf(dbout, "Debug: ... disgarding\n");
      }
    }
    sscanf(dstr, "%s %s %ld\n", dstr2, dstr3, &number_vertices);
    fprintf(dbout, "Debug: vertices: %ld\n", number_vertices);

    while(1) {
      fgets(dstr, 255, fp);
      fprintf(dbout, "Debug: reading input: %s", dstr);
      if (strncmp("element face", dstr, 12) == 0) {
	fprintf(dbout, "Debug: ... found the 'element face'\n");
	break;
      } else {
	//fprintf(dbout, "Debug: ... disgarding\n");
      }
    }
    sscanf(dstr, "%s %s %ld\n", dstr2, dstr3, &number_faces);
    fprintf(dbout, "Debug: faces: %ld\n", number_faces);

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

  fprintf(dbout, "Debug: vertices, faces: %ld %ld\n", number_vertices, number_faces);

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
    //  printf("HELP!\n");
    //  error(-1);
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



myfloat distance2(x,y,z, j)
myfloat x,y,z;
long j;
{
	myfloat dd;

	dd = (x-global_Sphere_Cx[j])*(x-global_Sphere_Cx[j])
	   + (y-global_Sphere_Cy[j])*(y-global_Sphere_Cy[j])
	   + (z-global_Sphere_Cz[j])*(z-global_Sphere_Cz[j]);
	return(dd);
}

myfloat distance(i, j)
long i,j;
{
	myfloat dd;

	dd = (vertex[i].x-global_Sphere_Cx[j])*(vertex[i].x-global_Sphere_Cx[j])
	   + (vertex[i].y-global_Sphere_Cy[j])*(vertex[i].y-global_Sphere_Cy[j])
	   + (vertex[i].z-global_Sphere_Cz[j])*(vertex[i].z-global_Sphere_Cz[j]);
	return(sqrt(dd));
}



//#define boundingSpheres boundingSpheresCentroid
#define boundingSpheres boundingSpheresOptimal

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
      //printf("obtuse triangle, use longest side\n");
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
      //printf("acute triangle, use circumcenter\n");
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



/* compute bounding spheres using the centroid, easy algorithm but no
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
void ProjectOnSegment(c1, c2, c3, p1, p2, p3, q1, q2, q3)
myfloat *c1, *c2, *c3;
myfloat p1, p2, p3;
myfloat q1, q2, q3;
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
myfloat FindClosestPointToOneTri(a1,a2,a3, face_index, c1, c2, c3)
myfloat a1, a2, a3;
myfloat *c1, *c2, *c3;
long face_index;
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
    printf("Error: non-enumerated case, can this happen?\n");
#ifdef EXTENDEDPRECISION
    printf("lambda mu %Lg %Lg\n", lambda, mu);
#else
    printf("lambda mu %g %g\n", lambda, mu);
    printf("factor %g\n", factor);
    printf("a11,a22,a12=%g,%g,%g\n",a11,a22,a12);
    printf("det?=%g\n",a11*a22-a12*a12);
#endif
    //if (factor == 0) {
    //  
    //}
    //error("Error in CP to tri, unanticipated case");
    printf("Error in CP to tri, unanticipated case\n");
  }

  /* Calculate distance */
  // TODO: found my bug, no sqrt here, dd is dist squared
  // HORRIBLE HACK 2010-07-28
  if (isinf(factor)) {
    dd = 10000.0;
    printf("OMG: I cannot believe you just did that!\n");
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
myfloat FindClosestPointGlobally2(a1,a2,a3, c1, c2, c3)
myfloat a1, a2, a3;
myfloat *c1, *c2, *c3;
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
 *  preprosing: find the bounding spheres around each triangle.
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

	    // TODO: what is 640?
	    if ( snprintf(tupstr, TUPSTRSZ, "(%d,%d,%d)", i+640,j+640,k+640) >= TUPSTRSZ )
	      myerr("Couldn't write tuple string, increase TUPSTRSZ?");
	    //fprintf(dbout, "DEBUG: tupstr=\"%s\"\n", tupstr);

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
	      //fprintf(dbout, "DEBUG: gridPtCtr=%d\n", gridPtCtr);
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
	      gridpt = ep->data;
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
  printf("*** total loops: %ld\n", total_count);
  printf("*** made %ld calls to Find_CP_to_tri()\n", cc);
  printf("*** found %ld grid points\n", numgridpts);

  free(tupstr);
}





/*
 * output a matrix w to a file
 */
void outputMatrixToFile(matrix w, myfloat t, const char *fname)
{
  int i, j, k;
  FILE *fd;

  if ((fd = fopen(fname, "w")) == NULL)
    myerr_nix("Error writing matrix to file");

  for (i=0; i<NPOINTS; i++) {
    for (j=0; j<NPOINTS; j++) {
      for (k=0; k<NPOINTS; k++) {
#ifdef EXTENDEDPRECISION
	fprintf(fd, "%.21Le\n", w[i][j][k]);
#else
	fprintf(fd, "%.18e\n", w[i][j][k]);
#endif
      }
    }
  }

  if ( fclose(fd) != 0 )  myerr_nix("Error closing file");
}



/*
 * output a matrix w to a file
 */
void outputGridToFile(const char *fname, int withPruning)
{
  long c, d = 0;
  FILE *fd;
  //int i,j,k;

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
  printf("Saved %ld gridpts to file\n", d);
}



/*
 * output a matrix w to a file
 */
void outputBoundingSpheres(const char *fname)
{
  long i;
  FILE *fd;

  if ((fd = fopen(fname, "w")) == NULL) myerr_nix("Error writing griddata");
  for (i=0; i<number_faces; i++) {
    // TODO: output
  }
  if ( fclose(fd) != 0 )  myerr_nix("Error closing file");

}



/*
 * allocate a new 3D NPOINT x NPOINT x NPOINT matrix
 */
#ifdef ALSOTHREEDMATRICES
matrix mallocMatrix()
{
  matrix w;
  int i, j;

  w = (myfloat ***) malloc((NPOINTS+1) * sizeof(myfloat **));
  if (w == NULL) myerr_nix("Error: malloc() failed");
  for(i=0; i<NPOINTS; i++) {
    w[i] = (myfloat **) malloc((NPOINTS+1) * sizeof(myfloat *));
    if (w[i] == NULL) myerr_nix("error: malloc() failed");
    for(j=0; j<NPOINTS; j++) {
      w[i][j] = (myfloat *) malloc((NPOINTS+1) * sizeof(myfloat));
      if (w[i][j] == NULL) myerr_nix("error: malloc() failed");
    }
  }
  return(w);
}



/*
 * free a 3D matrix
 */
void freeMatrix(matrix w)
{
  int i, j;

  for(i=0; i<NPOINTS; i++) {
    for(j=0; j<NPOINTS; j++) {
      free(w[i][j]);
    }
    free(w[i]);
  }
  free(w);
}
#endif


int main(void)
{
  long i;
  long expectedgridsz, hashtablesz;

  printf("%g, %d\n", -0.0, iafz(-0.0));
  printf("%g, %d\n", 0.0, iafz(0.0));
  printf("%g, %d\n", -0.1, iafz(-0.1));
  printf("%g, %d\n", 0.1, iafz(0.1));
  printf("%g, %d\n", -1.1, iafz(-1.1));
  printf("%g, %d\n", 1.1, iafz(1.1));


#ifdef ALSOTHREEDMATRICES
  fprintf(dbout, "Debug: mallocing 3D matrices... ");
  CPx = mallocMatrix();
  CPy = mallocMatrix();
  CPz = mallocMatrix();
  CPdd = mallocMatrix();
  fprintf(dbout, "...done\n");
#endif


  //printf("initShape()...\n");
  //initShape("annies_pig.ply");
  //initShape("pig_loop2.ply");
  //initShape("beacon.ply");
  //initShape("test2.ply");
  //initShape("bunny_input.ply");
  initShape("bumpy_torus_scaled.ply");

  /* find max,min vertices */
  myfloat min1, max1;
  min1 = 1e42;
  max1 = -1e42;
  for (i=0; i<number_vertices; i++) {
    if (vertex[i].x < min1) min1 = vertex[i].x;
    if (vertex[i].y < min1) min1 = vertex[i].y;
    if (vertex[i].z < min1) min1 = vertex[i].z;
    if (vertex[i].x > max1) max1 = vertex[i].x;
    if (vertex[i].y > max1) max1 = vertex[i].y;
    if (vertex[i].z > max1) max1 = vertex[i].z;
  }
  fprintf(dbout, "TODO: hash-table size estimator currently assumes model is roughly centered around origin\n");
  myfloat RAD;
  RAD = fmax(fabs(min1),fabs(max1));
#ifdef EXTENDEDPRECISION
  printf("dx=%Lg, BW=%Lg\n", DX, BANDWIDTH);
  printf("min1, max1=(%Lg,%Lg), RAD=%Lg\n", min1,max1, RAD);
  printf("pi=%.21Lg\n", CONST_PI);
#else
  printf("dx=%g, BW=%g\n", DX, BANDWIDTH);
  printf("min1, max1=(%g,%g), RAD=%g\n", min1,max1, RAD);
  printf("pi=%.18g\n", CONST_PI);
#endif


  //expectedgridsz = ceil((BANDWIDTH*DX)*(2*CONST_PI/DX));
  expectedgridsz = EXPECTEDGRIDSZ_SAFETY_FACTOR * \
    ceil((BANDWIDTH/DX)*(4*CONST_PI*(RAD)*(RAD)/(DX*DX)));
  hashtablesz = ceil(1.25*expectedgridsz);
  printf("(over)-estimated number of gridpoints: %ld\n", expectedgridsz);
  printf("hashtable size: %ld\n", hashtablesz);


  fprintf(dbout, "Debug: mallocing lists... ");
  gridPtList = (struct struct_cp **) malloc(expectedgridsz*sizeof(struct struct_cp *));
  if (gridPtList == NULL) myerr_nix("Error: malloc() grid list failed, msg");
  for (i = 0; i < expectedgridsz; i++)
    gridPtList[i] = (struct struct_cp *) NULL;

  if ((keyList = (char **) malloc(expectedgridsz*sizeof(char *))) == NULL)
    myerr_nix("Error: malloc() key list failed, msg");
  for (i = 0; i < expectedgridsz; i++)
    keyList[i] = (char *) NULL;

  fprintf(dbout, "...done\n");

  /* Create a new hash table.  To be portable we only can make one of
     these.  On GNU, could use hcreate_r(). */
  if (hcreate(hashtablesz) == 0) myerr_nix("Error creating hash table");


  printf("Debug: running boundingSpheres()...\n");
  boundingSpheres();
  printf("Debug: Finding Grid...\n");
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

  printf("Debug: Outputting to file...\n");
  outputGridToFile("Griddata", 1);

#ifdef ALSOTHREEDMATRICES
  outputMatrixToFile(CPx, 0, "CPdatax");
  outputMatrixToFile(CPy, 0, "CPdatay");
  outputMatrixToFile(CPz, 0, "CPdataz");
  outputMatrixToFile(CPdd, 0, "CPdatadd");
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

  return(0);
}

