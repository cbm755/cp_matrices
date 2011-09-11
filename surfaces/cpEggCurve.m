function [cpx, cpy, dist] = cpEggCurve(x, y)
%CPEGGCURVE  Closest Point function for the egg-shaped curve
%
%   Depends on the spline toolbox
%   Uses cpSpline2DClosed

  % interactive tool to make a new spline
  %[xy, spcv] = getcurve

  scale = 1.6917222947492;
  origPts = [...
      0          0.7 ;
      0.4        0.1 ;
      0.2       -0.7 ;
     -0.4       -0.5 ;
     -0.4        0.2 ;
      0          0.7 ];

  pts = scale*origPts;

  myspline = cscvn(pts');

  [cpx, cpy, dist] = cpSpline2D(x, y, myspline, 1);
