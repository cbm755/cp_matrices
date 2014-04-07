function [cpx,cpy,cpz,dist,bdy,s] = cpLineSegment3d(x,y,z, p, q)
%CPLINESEGMENT3D  Closest Point function for a line segment in 3D
%   [cpx,cpy,cpz, dist, bdy, s] = cpLineSegment3d(x,y,z, pt1, pt2)
%
%   Return closest points for a line segment between 'pt1' and 'pt2'
%
%   (outputs 'dist', 'bdy' and 's' are optional).  's' is the
%   parameter of the line.

  if (nargin == 3)
    p = [0 0 0];
    q = [1 0 0];
  end

  [cpx0,cpy0,cpz0, dist0, bdy0, s0] = cpLine(x,y,z, q-p, p);

  cpx = ((s0 >= 0) & (s0 <= 1)) .* cpx0 + ...
        (s0 < 0) .* (p(1)) + (s0 > 1) .* (q(1));
  cpy = ((s0 >= 0) & (s0 <= 1)) .* cpy0 + ...
        (s0 < 0) .* (p(2)) + (s0 > 1) .* (q(2));
  cpz = ((s0 >= 0) & (s0 <= 1)) .* cpz0 + ...
        (s0 < 0) .* (p(3)) + (s0 > 1) .* (q(3));

  if (nargout >= 4)
    dist = sqrt( (x-cpx).^2 + (y-cpy).^2 + (z-cpz).^2 );
  end
  if (nargout >= 5)
    bdy = (s0 < 0) .* 1  +  (s0 > 1) .* 2;
  end
  if (nargout >= 6)
    s = ((s0 >= 0) & (s0 <= 1)) .* s0  +  ...
        (s0 < 0) .* zeros(size(s0))  +  ...
        (s0 > 1) .* ones(size(s0));
  end
