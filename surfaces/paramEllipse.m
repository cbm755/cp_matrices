function [x,y,th] = paramEllipse(n, aa, bb, cen)
%PARAMELLIPSE  A parameterization of an ellipse
%   [x,y] = paramEllipse(N,A,B) returns an N-point mesh for an ellipse
%   centered at the origin with major axis A and minor axis B.
%
%   [x,y] = paramEllipse(N,A,B,CEN) returns a mesh for an ellipse
%   centered at CEN.
%
%   [x,y,th] = paramEllipse(...) returns the value of the parameter
%   as well.

  % defaults
  if (nargin < 3)
    if (nargin == 2)
      error('must provide both or either of a,b');
    end
    aa = 1.5;
    bb = 0.75;
  end
  if (nargin < 4)
    cen = [0 0];
  end

  %th = ( 0:(2*pi/n):(2*pi-2*pi/n) )';
  th = ( 0:(2*pi/n):(2*pi) )';
  x = aa*cos(th) + cen(1);
  y = bb*sin(th) + cen(2);
