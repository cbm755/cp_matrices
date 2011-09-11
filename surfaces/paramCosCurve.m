function [x,y,tt] = paramCosCurve(N)
%PARAMCOSCURVE   A parameterization of the cosine curve

  % TODO: The spline data here is duplicated in the cp file, should use
  % a common helper function?

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

  sp = myspline;
  tt = linspace(sp.breaks(1), sp.breaks(end), N);

  A = ppval(myspline, tt);
  x = transpose(A(1,:));
  y = transpose(A(2,:));
