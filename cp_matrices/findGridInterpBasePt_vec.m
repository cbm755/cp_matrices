function [I,X] = findGridInterpBasePt_vec(x, p, relpt, dx)
%FINDGRIDINTERPBASEPT_VEC   Find the "base point" for an interpolation
% This is best explained in the diagram below.
%
% [I,X] = findGridInterpBasePt(x, p, relpt, dx)
%   x: is the interpolation point, must lie inside the === signs
%      below.  (x should be Nx3 matrix).
%   p: degree interpolation (N-1 point interp), a scalar.
%   relpt: a reference point, corresponding to (1,1) in your grid (a
%          1x3 row vector).
%   dx: grid spacing (a 1x3 row vector).
%   X: is the "basepoint", B below, the lower left corner of an interpolation
%      stencil hypercube
%   I: 2D/3D/etc index of X, measured relative to 'relpt'
%
% picture:
% p=0: ==B==
% p=1:   B====x
% p=2:   B  ==x==  x
% p=3:   B    x====x    x
% p=4:   B    x  ==x==  x    x
% p=5:   B    x    x====x    x    x
% etc
%
% Notes:
% * This index I is base 1 (i.e., B=(1,1,...) is the bottom left of
%   the whole grid)
% * Seems to work for negative I too...
% * TODO: could worry about using doubles for indexing.  This could
%   matter when indices become as large as 1e16 (note int64(1e15+1) vs
%   int64(1e16+1)).  sub2ind could get that large if we had a cube
%   with something like 10000 points per side in 3D (or less in 4D).
%   This is not as ridiculous and it sounds: number of grid points
%   actually used in the Closest Point Method scales like 2D (for a 2D
%   surface).  For now I'll add a warning to the 3D code.

  % TODO: could simplify, if isscalar(dx)

  if (mod(p,2) == 0)  % even
    I = round( ( x - repmat(relpt,size(x,1),1) ) ./ repmat(dx,size(x,1),1) ) + 1  -  p/2;
  else % odd
    I = floor( ( x - repmat(relpt,size(x,1),1) ) ./ repmat(dx,size(x,1),1) ) + 1  -  (p-1)/2;
  end
  % I-1 here because I is a matlab index
  X = repmat(relpt,size(x,1),1) + (I-1).*repmat(dx,size(x,1),1);

  % in Matlab R14, isinteger is really slow, seems ok in 2010b
  %if ~(isinteger(p))
  %  p = int64(p)
  %end
  %if (mod(p,2) == 0)  % even
  %  I = int64(round( (x-relpt) ./ dx )) + 1  -  p/2;
  %else % odd
  %  I = int64(floor( (x-relpt) ./ dx )) + 1  -  (p-1)/2;
  %end
  % I-1 here because I is a matlab index
  %X = relpt + double(I-1).*dx;
end
