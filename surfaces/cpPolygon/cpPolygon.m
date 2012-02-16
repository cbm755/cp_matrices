function [cpx, cpy, sdist] = cpPolygon(x, y, poly)
%CPPOLYGON  Closest Point function for a polygon
%   The polygon does not need to be convex.
%
%   [cpx, cpy] = cpPolygon(x, y, poly)
%   [cpx, cpy, sdist] = cpPolygon(x, y, poly)
%
%   You can generate 'poly' interactively by calling
%
%      poly = makePolyInteractive()
%
%   If you leave off poly it will use a default one.
%
%   Note: returns signed distance (with negative inside).
%   TODO: currently sign of signed distance seems to depend on
%   orientation of the bdy: should fix that!
%
%   Uses code by Tom Maerz.

  % defaults
  if (nargin < 3)
    poly = makePolyDefault();
  end

  [cpx,cpy] = helper_CPopPoly(x, y, poly);
  sdist = helper_sDistPoly(x, y, poly);
  sdist = -sdist;
