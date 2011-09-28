/*
 * helper_tri2cp.c can optionally build some dense 3D matrices to help
 * with debugging.  This provides some routines for that purpose.
 */

#ifndef MAT3D_H
#define MAT3D_H

#define matrix myfloat ***
matrix CPx;
matrix CPy;
matrix CPz;
matrix CPdd;

int NPOINTS;

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
 * allocate a new 3D NPOINT x NPOINT x NPOINT matrix
 */
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
