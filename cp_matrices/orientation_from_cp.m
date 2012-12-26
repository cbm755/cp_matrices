function [in2,sdist] = orientation_from_cp(xx,yy,zz,cpx,cpy,cpz,dist,dx,seedin,seedout,verbose)
%ORIENTATION_FROM_CP  Orient a cp-representation from seed points
%   Given a closest point representation of a closed surface
%   embedded in 3D, classify each grid points as inside/outside
%   based on a few seed points.
%
%   [inside,sdist] = orientation_from_cp(xx,yy,zz, cpx,cpy,cpz, ...
%                        dist,dx, seedin, seedout, verbose)
%
%   'seedin' and 'seedout' should be lists of indices into the arrays
%   xx, yy, zz, dist.  They must be points which are known to be
%   inside/outside.
%
%   This code starts at the points in 'seedin' and flood fills toward
%   the surface.  Then it does the same from 'seedout'.
%
%   'inside' will be the same shape as xx and contains ones for each
%   points classified as inside.  'sdist' is signed distance.
%
%     verbose = 0: quiet (default if omitted)
%     verbose = 1: progress
%     verbose = 2: progress + plots (must have a meshgrid)
%
%   TODO: not robust: will fail, possibly badly, for open surfaces or
%   if there are holes.  Will probably need some tweaks for point
%   clouds.

  if (nargin < 11)
    verbose = 0;
  end

  inside = orientation_fill(xx,yy,zz,dist,dx,seedin,verbose);
  outside = orientation_fill(xx,yy,zz,dist,dx,seedout,verbose);

  % these tests only valid if we already know signed distance
  %test1 = all(dist(logical(inside)) < 0)
  %test2 = all(dist(logical(outside)) > 0)

  ni = nnz(inside);
  no = nnz(outside);
  total = length(xx(:));

  if ~(ni + no <= total)
    error('open or must have a leak');
  end

  noclass = find(~inside & ~outside);
  nnc = length(noclass);

  if ~(ni + no + nnc == total)
    error('something doesn''t add up!')
  end

  [in2,out2] = orientation_stage2(xx,yy,zz, ...
                                  cpx,cpy,cpz, dist, dx, ...
                                  inside, outside, verbose);

  noclass2 = find(~in2 & ~out2);
  nnc2 = length(noclass2);

  if ~(nnc2 == 0)
    error('some points unclassifed after stage 2');
  end

  sdist = -1*in2.*abs(dist) + out2.*abs(dist);

