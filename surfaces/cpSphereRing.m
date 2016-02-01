function [cpx,cpy,cpz,dist,bdy] = cpSphereRing(x,y,z,zlim,R,cen)
%CPSPHERERING   Closest point fcn for a ring, subset of sphere
%   This can be used to make Rings, Fishbowls and other things
%   that are subsets of a sphere
%
%   Code is vectorized: any size/shape for x should work.


  % defaults
  if (nargin < 5)
    R = 1;
  end
  if (nargin < 6)
    cen = [0,0,0];
  end
  if (nargin < 4)
    % default is a fishbowl
    zlim = [-2*R  0.5];
  end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);
  z = z - cen(3);
  zlim = zlim - cen(3);

  [cpx, cpy, cpz] = cpSphere(x, y, z, R);

  % should not include the bdy itself (hence <=)
  bdy1 = (cpz < zlim(1));
  bdy2 = (cpz > zlim(2));

  % Next we find the CPs to two circles, each centered at (0,0, zlim(i))
  % with appropriate radius.
  radlim = sqrt(R^2 - zlim.^2);

  [CP, dist1] = cpCircleInHighDim({x, y, z}, radlim(1), [0; 0; zlim(1)]);
  [cpx1, cpy1, cpz1] = deal(CP{:});

  [CP, dist2] = cpCircleInHighDim({x, y, z}, radlim(2), [0; 0; zlim(2)]);
  [cpx2, cpy2, cpz2] = deal(CP{:});

  % for those points whos closest point was outside the zlim range, we replace
  % them with the closest points on one of the two circles.
  cpx(bdy1) = cpx1(bdy1);
  cpy(bdy1) = cpy1(bdy1);
  cpz(bdy1) = cpz1(bdy1);

  cpx(bdy2) = cpx2(bdy2);
  cpy(bdy2) = cpy2(bdy2);
  cpz(bdy2) = cpz2(bdy2);

  bdy = (bdy1 | bdy2);

  dist = sqrt( (x-cpx).^2 + (y-cpy).^2 + (z-cpz).^2 );

  % shift back to center
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
  cpz = cpz + cen(3);

end


%!test
%! [cpx, cpy, cpz, dist, bdy] = cpSphereRing(0, 0, 10, [1.4 inf], 2);
%! assert(cpx, 0, 2*eps);
%! assert(cpy, 0, 2*eps);
%! assert(cpz, 2, 2*eps);
%! assert(dist, 8, 2*eps);
%! assert(bdy, false);

%!test
%! [cpx, cpy, cpz, dist, bdy] = cpSphereRing(1, 0, -10, [1.4 inf], 2);
%! R = sqrt(2^2 - 1.4^2);
%! assert(cpx, R, 5*eps);
%! assert(cpy, 0, 2*eps);
%! assert(cpz, 1.4, 2*eps);
%! assert(bdy, true);

%!test
%! [cpx, cpy, cpz, dist, bdy] = cpSphereRing(0, -1, -10, [1.4 inf], 2);
%! R = sqrt(2^2 - 1.4^2);
%! assert(cpx, 0, 2*eps);
%! assert(cpy, -R, 2*eps);
%! assert(cpz, 1.4, 2*eps);
%! assert(bdy, true);

%!test
%! % random points are in the ring
%! x = 10*rand(10, 10, 10) - 5;
%! y = 10*rand(10, 10, 10) - 5;
%! z = 10*rand(10, 10, 10) - 5;
%! [cpx, cpy, cpz, dist, bdy] = cpSphereRing(x, y, z, [1.4 1.5], 2, [1.1; 1.2; 1.3]);
%! assert(cpz >= 1.4 & cpz <= 1.5);
%! R1 = 4 - (1.5-1.3)^2;
%! R2 = 4 - (1.5-1.4)^2;
%! rad = (cpx-1.1).^2 + (cpy-1.2).^2;
%! assert(rad >= R1 - 10*eps)
%! assert(rad <= R2 + 10*eps)

%!test
%! % reduces to circle: random points
%! x = 10*rand(10, 10, 10) - 5;
%! y = 10*rand(10, 10, 10) - 5;
%! z = 10*rand(10, 10, 10) - 5;
%! [cpx, cpy, cpz, dist, bdy] = cpSphereRing(x, y, z, [1.5 1.5], 2, [1.1; 1.2; 1.3]);
%! R = sqrt(4 - (1.5-1.3)^2);
%! [CP, dist1] = cpCircleInHighDim({x y z}, R, [1.1; 1.2; 1.5]);
%! assert(dist, dist1, 10*eps);
%! assert(cpx, CP{1}, 10*eps);
%! assert(cpy, CP{2}, 10*eps);
%! assert(cpz, CP{3}, 10*eps);
