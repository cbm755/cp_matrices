function bw = rm_bandwidth(dim, p, fdrad, safety)
%RM_BANDWIDTH  Ruuth--Merriman bandwidth formula
%  The width of a computational domain for the Ruuth--Merriman
%  closest point method.  If you include grid points within 'bw*dx'
%  of the surface, the discrete operators will fit.
%
%  bw = rm_bandwidth(dim, p)
%     This gives the bandwidth for degree 'p' interpolation with
%     embedding dimension 'dim'.  Assuming dist is (possible
%     signed) distance of a superset of the band, we can do:
%  band = find(abs(dist) <= bw*dx);
%
%  bw = rm_bandwidth(dim, p, fd_rad, safety)
%     'fd_rad' specifies the "radius" of the finite difference scheme,
%     defaults to 1 (e.g., for the {-4,1,1,1,1} stencil for the
%     Laplacian or other first-order differences).
%     'safety' is a multiplicative safety factor which slightly
%     widens the band, mainly so that rounding error doesn't place
%     outside the band.  Defaults to 1.0001.
%
%  TODO: double-check what 'fd_sten' should be for other schemes
%  (e.g., with diagonal components).
%
%  Reference: [Ruuth & Merriman 2008, pg 1951]

  if (nargin < 4)
    safety = 1.0001;
  end
  if (nargin < 3)
    fdrad = 1;
  end
  if (nargin < 2)
    p = 3;
  end

  bw = safety*sqrt( (dim-1)*((p+1)/2)^2 + ((fdrad + (p+1)/2)^2) );

