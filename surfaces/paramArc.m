function [x, y, th] = paramArc(n, R, cen, angle1, angle2)
%PARAMARC  Parameterization of a circular arc.
%   [x, y, th] = paramArc(n, R, cen, angle1, angle2)
%   'n' specifies the approximate number of points to be returned.
%   'R', 'cen', 'angle1' and 'angle2' all take default values if
%   omitted.


  % defaults
  if (nargin < 2)
    R = 1;
  end
  if (nargin < 3)
    cen = [];
  end
  if (isempty(cen))
    cen = [0 0];
  end
  if (nargin < 4)
    angle1 = 0;
  end
  if (nargin < 5)
    angle2 = pi/2;
  end

  n = max(n, 3);

  % (-pi,pi]
  angle1 = angle(exp(1i*angle1));
  angle2 = angle(exp(1i*angle2));

  if (angle1 <= angle2);
    %disp('lte');
    a1 = angle1;
    a2 = angle2;
  else
    %disp('gt');
    a1 = angle1;
    a2 = (2*pi+angle2);
  end

  dth = ((a2-a1)/(n-1));
  th = (a1:dth:a2)';
  th = angle(exp(1i*th));

  x = cos(th);
  y = sin(th);

  x = R*x + cen(1);
  y = R*y + cen(2);

