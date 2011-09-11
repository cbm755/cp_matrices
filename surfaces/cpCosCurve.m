function [cpx, cpy, dist, bdy] = cpCosCurve(x, y)
%CPCOSCURVE  Closest Point function for cosine curve
%
%   Depends on the spline toolbox
%   Uses cpSpline2DClosed

  % interactive tool to make a new spline
  %[xy, spcv] = getcurve

  %scale = 2.46581390599295;  % 2*pi length
  scale = 1.23290695299646;  % pi  length
  pts = [...
      1.0      0.1 ;
      0.83     0.4 ;
      0.4      0.4 ;
      0.0     -0.1 ;
     -0.4     -0.3 ;
     -1.0     -0.1 ];

  pts = scale*pts;

  myspline = cscvn(pts');

  [cpx, cpy, dist, bdy] = cpSpline2D(x, y, myspline, 0);
