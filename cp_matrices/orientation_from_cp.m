function [in2,sdist,unknown] = orientation_from_cp(xx,yy,zz,cpx,cpy,cpz,dist,dx,E,seedin,seedout,verbose)
%ORIENTATION_FROM_CP  Orient a cp-representation from seed points
%   Given a closest point representation of a closed surface
%   embedded in 3D, classify each grid points as inside/outside
%   based on a few seed points.
%
%   WARNING: This should not be considered reliable!  "stage 2" of
%   the algorithm is heurestic and a bit ad hoc.
%
%   [inside,sdist,unknown] = orientation_from_cp(xx,yy,zz, ...
%               cpx,cpy,cpz, dist,dx, E, seedin, seedout, verbose)
%
%   'seedin' and 'seedout' should be lists of indices into the arrays
%   xx, yy, zz, dist.  They must be points which are known to be
%   inside/outside.
%
%   'E' is the extension matrix, used in stage 2 of the algorithm.
%
%   This code starts at the points in 'seedin' and flood fills toward
%   the surface.  Then it does the same from 'seedout'.
%
%   'inside' will be the same shape as xx and contains ones for each
%   points classified as inside.  'sdist' is signed distance.
%
%     verbose = 0: very quiet
%     verbose = 1: quiet (default if omitted)
%     verbose = 2: progress
%     verbose = 3: progress + plots (must have a meshgrid)
%
%   'unknown' contains a list of points which the algorithm was
%   unable to classify.
%
%   Not robust: will fail, possibly badly, for open surfaces or if
%   there are holes.  Will probably need some tweaks for point clouds.
%   Certainly has trouble near regions of high-curvature (relative to
%   the grid spacing).  Using a wider band can help.
%
%   TODO: drops to keyboard eventually for unknown points: this
%   should be a user option.
%
%   TODO: currently 3D only although should be easy to port to 2D.

  if (nargin < 12)
    verbose = 1;
  end

  if (verbose >= 1)
    disp('Orientation: starting fill from inside');
  end
  inside = orientation_fill(xx,yy,zz,dist,dx,seedin,verbose);
  if (verbose >= 1)
    disp('Orientation: starting fill from outside');
  end
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

  if verbose >= 1
    fprintf('Orientation: flood fill classified %d of %d: %d remain\n', ...
            ni+no, total, nnc);
  end

  if ~(ni + no + nnc == total)
    error('something doesn''t add up!')
  end

  [in2,out2,unknown] = orientation_stage2(xx,yy,zz, ...
                            cpx,cpy,cpz, dist, dx, E, ...
                            inside, outside, verbose);

  noclass2 = find(~in2 & ~out2);
  nnc2 = length(noclass2);

  if ~(nnc2 == 0)
    error('some points unclassifed after stage 2');
  end

  sdist = -1*in2.*abs(dist) + out2.*abs(dist);

