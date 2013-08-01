function [cpx, dist] = cpSphereInHighDim(x, R, cen)
%CPSPHEREINHIGHDIM   A sphere embedded in a R^d where d>=3

  if ~iscell(x)
    error('expected cell array for input x');
  end

  dim = length(x);

  if (nargin < 2)
    R = 1;
  end
  if (nargin < 3)
    cen = zeros(1,dim);
  end

  cpx = cell(1,dim);

  [th, phi, tilde] = cart2sph(x{1}-cen(1), x{2}-cen(2), x{3}-cen(3));
  [cpx{1},cpx{2},cpx{3}] = sph2cart(th, phi, R);

  for d = 1:3
    cpx{d} = cpx{d} + cen(d);
  end
  % extra dimension are just constants
  for d = 4:dim
    cpx{d} = cen(d) * ones(size(cpx{1}));
  end

  dist = (x{1}-cpx{1}).^2;
  for d = 2:dim
    dist = dist + (x{d}-cpx{d}).^2;
  end
  dist = sqrt(dist);

