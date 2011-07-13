function [Ix,D1xx,D1xc,D1xb,D1xf] = diff_matrices1d(N, dx)
%DIFF_MATRICES1D  Build 1D operators with periodic BCs
%
%   todo: support other BCs?

  Ix = speye(N,N);

  e = ones(N,1);
  D1xx = spdiags([e  -2*e  e], [-1 0 1], N, N);
  D1xx(1,N) = 1;
  D1xx(N,1) = 1;
  D1xx = D1xx/dx^2;

  D1xc = spdiags([-e  e], [-1 1], N, N);
  D1xc(1,N) = -1;
  D1xc(N,1) = 1;
  D1xc = D1xc/(2*dx);

  D1xb = spdiags([-e  e], [-1 0], N, N);
  D1xb(1,N) = -1;
  %D1xb(N,1) = 1;
  D1xb = D1xb/dx;

  D1xf = spdiags([-e  e], [0 1], N, N);
  %D1xf(1,N) = -1;
  D1xf(N,1) = 1;
  D1xf = D1xf/dx;
