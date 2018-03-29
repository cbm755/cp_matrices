function [cpx, cpy, dist, bdy] = cpEllipticalDisc(x, y, varargin)
%cpEllipticalDisc  Closest Point function for a elliptical-shaped disc.
%
%   >> [cpx, cpy, dist, bdy] = cpEllipticalDisc([0 2 0], [0 0 2], 2, 1)
%   cpx = 0  2  0
%   cpy = 0  0  1
%   dist = 0
%   bdy = 0
%


  [cpx, cpy, dist] = cpEllipse(x, y, varargin{:});

  I = dist < 0;

  cpx(I) = x(I);
  cpy(I) = y(I);
  dist(I) = 0;

  bdy = ~I;
