function [cpx, dist] = cpCircleInHighDim(x, R, cen)
%CPCIRCLEINHIDIM  A circle embedded in a higher dimensional space
%   closest point to a circle in the x,y plane embedded in (at
%   least R^2).
%
%   [cpx, dist] = cpCircleInHighDim(x, radius, center)

  if ~iscell(x)
    error('expected cell array for input x');
  end

  dim = length(x);

  [th, r] = cart2pol(x{1}-cen(1), x{2}-cen(2));
  [cpx{1}, cpx{2}] = pol2cart(th, R);

  for d = 1:2
    cpx{d} = cpx{d} + cen(d);
  end
  % extra dimension are just constants
  for d = 3:dim
    cpx{d} = cen(d) * ones(size(cpx{1}));
  end

  dist = (x{1}-cpx{1}).^2;
  for d = 2:dim
    dist = dist + (x{d}-cpx{d}).^2;
  end
  dist = sqrt(dist);


  %[th, r] = cart2pol(xx{1}, xx{2});
  %[cpx, cpy] = pol2cart(th, 1);
  %cpz = ZC*ones(size(cpx));
  %cpw = WC*ones(size(cpx));
  %cpX = {cpx cpy cpz cpw};
  %dist = (xx{1}-cpX{1}).^2;
  %for d = 2:dim
  %  dist = dist + (xx{d}-cpX{d}).^2;
  %end
  %dist = sqrt(dist);

