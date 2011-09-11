function [cpx, cpy, dist] = cpBeanCurve(x, y)
%CPBEANCURVE  Closest Point function for a bean-shaped curve
%
%   Depends on the spline toolbox
%   Uses cpSpline2DClosed


  % TODO: The spline data here is duplicated in the param file, should
  % use a common helper function?


  % interactive tool to make a new spline:
  %[xy, spcv] = getcurve

  scale = 1.0;
  pts = [...
         -0.4          0.42
         -0.6          0.4
         -0.8          0.25
         -0.9            0
         -0.8         -0.3
         -0.6         -0.43
         -0.3         -0.5
          0.1         -0.5
          0.4         -0.44
          0.7         -0.25
          0.8            0
          0.75          0.24
          0.6          0.35
          0.4          0.35
          0.17         0.25
         -0.05         0.25
         -0.22          0.35
         -0.4          0.42];

  % make it roughly centered at origin
  pts = pts + repmat([0.05 0.04], [size(pts,1), 1]);
  pts = scale*pts;

  myspline = cscvn(pts');

  DEBUG = 0;
  [cpx, cpy, dist] = cpSpline2D(x, y, myspline, 1, DEBUG);