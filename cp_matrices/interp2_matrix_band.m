function E = interp2_matrix_band(x, y, xi, yi, p, band, use_ndgrid)
%INTERP2_MATRIX_BAND  Return a interpolation matrix over a band
%   E = INTERP2_MATRIX(X,Y,XI,YI,P,BAND)
%   Build a matrix which interpolates grid data on a grid defined
%   by the product of the lists X and Y onto the
%   points specified by the lists XI and YI.  Interpolation is
%   done using degree P barycentric Lagrange interpolation.
%
%   BAND is a list of linear indices into a (possibly fictious) 2D
%   array of points constructed with meshgrid.
%
%   MESHGRID is important: this code makes assumptions about
%   orderings of the grid that are true for MESHGRID
%
%   Works by calling INTERP2_MATRIX and then keeping only those
%   columns that correspond to BAND.
%
%   Does no error checking up the equispaced nature of x,y,z

  if (nargin < 7)
    use_ndgrid = false;
  end

  warning('this code is deprecated, call interp2_matrix directly');

  E = interp2_matrix_oldloop(x, y, xi, yi, p, use_ndgrid);

  % sanity check: the columns outside of band should all be zero
  Nx = length(x);
  Ny = length(y);
  Eout = E(:,setdiff(1:(Nx*Ny),band));
  if (nnz(Eout) > 0)
    nnz(Eout)
    warning('Lost some non-zero coefficients (from outside the innerband)');
  end

  % remove columns that are outside of band
  E = E(:,band);
