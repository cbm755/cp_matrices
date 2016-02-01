function [cpx, dist] = cpCircleInHighDim(x, R, cen)
%CPCIRCLEINHIDIM  A circle embedded in a higher dimensional space
%   closest point to a circle in the x,y plane embedded in (at
%   least R^2).
%
%   [cp, dist] = cpCircleInHighDim(x, radius, center)
%
%   Example:
%   >> [cp, dist] = cpCircleInHighDim({3 0 0}, 1, [0 0 0])
%      cp =
%      {
%        [1,1] = 1
%        [1,2] = 0
%        [1,3] = 0
%      }
%      dist = 2


  if ~iscell(x)
    error('expected cell array for input x');
  end

  dim = length(x);

  cpx = cell(1,dim);

  [th, tilde] = cart2pol(x{1}-cen(1), x{2}-cen(2));
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

end


%!test
%! [cp, dist] = cpCircleInHighDim({2, 0, 1}, 1, [0; 0; 1]);
%! assert(cp, {1, 0, 1});
%! assert(dist, 1);

%!test
%! [cp, dist] = cpCircleInHighDim({-2, 0, 1}, 10, [0; 0; 1]);
%! assert(cp, {-10, 0, 1}, 10*eps);
%! assert(dist, 8);

%!test
%! [cp, dist] = cpCircleInHighDim({0, 0, 2}, 1, [0; 0; 1]);
%! assert(dist, sqrt(2), eps);

%!test
%! [cp, dist] = cpCircleInHighDim({2, 2, 4}, 1, [1; 2; 3]);
%! assert(cp, {2, 2, 3});
%! assert(dist, 1);

%!test
%! % random points go to correct plane
%! x = 10*rand(3, 3, 3) - 5;
%! y = 10*rand(3, 3, 3) - 5;
%! z = 10*rand(3, 3, 3) - 5;
%! [cp, dist] = cpCircleInHighDim({x, y, z}, 1.3, [1.1; 2.2; 3.3]);
%! assert(cp{3}, 3.3*ones(size(x)));
%! assert((cp{1}-1.1).^2 + (cp{2}-2.2).^2, 1.3^2*ones(size(x)), 5*eps);

%!test
%! % compare to specific R^4 implementation
%! % first, random data
%! cen = rand(4, 1);
%! R = 1.2*rand;
%! for i=1:4
%!   x{i} = 2*rand(5, 5, 5) - 1;
%! end
%! cen = rand(4, 1);
%! R = 1.2*rand;
%! [CP, DIST] = cpCircleInHighDim(x, R, cen);
%! % now, the other implementation
%! [th, r] = cart2pol(x{1} - cen(1), x{2} - cen(2));
%! [cpx, cpy] = pol2cart(th, R);
%! cpx = cpx + cen(1);
%! cpy = cpy + cen(2);
%! cpz = cen(3)*ones(size(cpx));
%! cpw = cen(4)*ones(size(cpx));
%! cp = {cpx cpy cpz cpw};
%! dist = (x{1}-cp{1}).^2;
%! for d = 2:4
%!   dist = dist + (x{d}-cp{d}).^2;
%! end
%! dist = sqrt(dist);
%! assert(cp, CP);
%! assert(dist, DIST);
