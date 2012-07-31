function E = interp3_matrix_band(x, y, z, xi, yi, zi, p, band, use_ndgrid)
%INTERP3_MATRIX_BAND  Return a interpolation matrix over a band
%   E = INTERP3_MATRIX(X,Y,Z,XI,YI,ZI,P,BAND)
%   Build a matrix which interpolates grid data on a grid defined
%   by the product of the lists X, Y and Z onto the
%   points specified by the lists XI, YI, ZI.  Interpolation is
%   done using degree P barycentric Lagrange interpolation.
%
%   BAND is a list of linear indices into a (possibly fictious) 3D
%   array of points constructed with meshgrid.
%
%   Works by calling INTERP3_MATRIX and then keeping only those
%   columns that correspond to BAND.
%
%   Does no error checking up the equispaced nature of x,y,z

  if (nargin < 9)
    use_ndgrid = false;
  end

  warning('this code is deprecated, call interp3_matrix directly');

  E = interp3_matrix_oldloop(x, y, z, xi, yi, zi, p, use_ndgrid);

  % sanity check: the columns outside of band should all be zero
  Nx = length(x);
  Ny = length(y);
  Nz = length(z);
  Eout = E(:,setdiff(1:(Nx*Ny*Nz),band));
  if (nnz(Eout) > 0)
    nnz(Eout)
    warning('Lost some non-zero coefficients (from outside the innerband)');
  end

  E = E(:,band);
