function [x,y,tt] = paramEggCurve(N)
%PARAMEGGCURVE   A parameterization of the egg-shaped curve


  % TODO: The spline data here is duplicated in the cp file, should use
  % a common helper function?

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

  %myspline.isPeriodic = 1;

  sp = myspline;
  tt = linspace(sp.breaks(1), sp.breaks(end), N);

  A = ppval(myspline, tt);
  x = transpose(A(1,:));
  y = transpose(A(2,:));
