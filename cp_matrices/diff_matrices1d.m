function [Ix,D1xx,D1xc,D1xb,D1xf] = diff_matrices1d(N, dx, BC)
%DIFF_MATRICES1D  Build 1D difference operators (default periodic BCs)
%   [Ix,D1xx,D1xc,D1xb,D1xf] = diff_matrices1d(N, dx)
%      N is length and dx is the grid spacing, with periodic BCs.
%
%   [Ix,D1xx,D1xc,D1xb,D1xf] = diff_matrices1d(N, dx, 'n')
%      Homogeneous Neumann BCs, using mirrored ghost values.
%
%   [Ix,D1xx,D1xc,D1xb,D1xf] = diff_matrices1d(N, dx, 'd')
%      Homoegeneous Dirichlet BCs, using mirrored ghost values.
%
%   TODO: there are many ways of imposing BC: these may not be choices
%   you want...  The Neumann conditions for forward/backward
%   differences is the one that's consistent with half-point
%   evaluations.

  if nargin < 3
    BC = 'p'
  end

  switch BC
    case 'p'  % periodic BCs
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
      D1xb = D1xb/dx;

      D1xf = spdiags([-e  e], [0 1], N, N);
      D1xf(N,1) = 1;
      D1xf = D1xf/dx;
    case 'n'  % neumann BCs
      Ix = speye(N,N);
      e = ones(N,1);
      D1xx = spdiags([e  -2*e  e], [-1 0 1], N, N);
      D1xx(1,2) = 2;
      D1xx(N,N-1) = 2;
      D1xx = D1xx/dx^2;

      D1xc = spdiags([-e  e], [-1 1], N, N);
      D1xc(1,2) = 0;
      D1xc(N,N-1) = 0;
      D1xc = D1xc/(2*dx);

      D1xb = spdiags([-e  e], [-1 0], N, N);
      D1xb(1,2) = -1;
      %D1xb(1,1) = 0;  D1xb(N,N) = 0;  D1xb(N,N-1) = 0;
      D1xb = D1xb/dx;

      D1xf = spdiags([-e  e], [0 1], N, N);
      D1xf(N,N-1) = 1;
      %D1xf(N,N) = 0;  D1xf(1,1) = 0;  D1xf(1,2) = 0;
      D1xf = D1xf/dx;

    case 'd'  % dirichlet BCs
      Ix = speye(N,N);
      e = ones(N,1);
      % or maybe user would want one-sided 2nd-difference...
      D1xx = spdiags([e  -2*e  e], [-1 0 1], N, N);
      %D1xx(1,1) = 0;
      D1xx(1,2) = 0;
      %D1xx(N,N) = 0;
      D1xx(N,N-1) = 0;
      D1xx = D1xx/dx^2;

      D1xc = spdiags([-e  e], [-1 1], N, N);
      D1xc(1,2) = 2;
      D1xc(N,N-1) = -2;
      D1xc = D1xc/(2*dx);

      D1xb = spdiags([-e  e], [-1 0], N, N);
      D1xb(1,2) = 1;
      D1xb = D1xb/dx;

      D1xf = spdiags([-e  e], [0 1], N, N);
      D1xf(N,N-1) = -1;
      D1xf = D1xf/dx;

  end

