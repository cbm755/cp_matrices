function [x,y,tt] = paramBeanCurve(N)
%PARAMBEANCURVE   A parameterization of the bean curve


  % TODO: The spline data here is duplicated in the cp file, should use
  % a common helper function?

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

  %myspline.isPeriodic = 1;

  sp = myspline;
  tt = linspace(sp.breaks(1), sp.breaks(end), N);

  A = ppval(myspline, tt);
  x = transpose(A(1,:));
  y = transpose(A(2,:));
