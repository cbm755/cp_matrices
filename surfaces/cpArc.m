function [cpx, cpy, dist, bdy] = cpArc(x, y, R, cen, angle1, angle2)
%CPARC  Closest Point function for a circular arc.
%   [cpx, cpy, dist, bdy] = cpArc(x, y, R, cen, angle1, angle2)
%      A circular arc of radius R centered at 'cen'.  The angle is
%      in [angle1, angle2].
%   'R', 'cen', 'angle1' and 'angle2' all take default values if
%   omitted.
%
%   Code is vectorized: any size/shape for x should work.

  % defaults
  if (nargin < 3)
    R = 1;
  end
  if (nargin < 4)
    cen = [];
  end
  if (isempty(cen))
    cen = [0 0];
  end
  if (nargin < 5)
    angle1 = 0;
  end
  if (nargin < 6)
    angle2 = pi/2;
  end
  % (-pi,pi]
  angle1 = angle(exp(1i*angle1));
  angle2 = angle(exp(1i*angle2));
  %  if (angle2 < angle1)
  %  temp = angle1;
  %  angle1 = angle2;
  %  angle2 = angle1;
  %end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);

  [th, tilde] = cart2pol(x, y);

  %if ( (th >= angle1) && (th <= angle2) )
  [cpx0,cpy0] = pol2cart(th, R);
  %bdy0 = zeros(size(th));
  %dist0 = sqrt( (x-cpx0).^2 + (y-cpy0).^2 );

  [cpx1,cpy1] = pol2cart(angle1, R);
  [cpx2,cpy2] = pol2cart(angle2, R);
  dist1 = sqrt( (x-cpx1).^2 + (y-cpy1).^2 );
  dist2 = sqrt( (x-cpx2).^2 + (y-cpy2).^2 );

  %where0 = ( (th >= angle1) & (th <= angle2) );
  where0 = ( (angle2 >= angle1) & ( (th >= angle1) & (th <= angle2) ) )  |  ...
           ( (angle2 <  angle1) & ( (th >= angle1) | (th <= angle2) ) );
  where1 = dist1 < dist2;

  cpx = where0 .* cpx0  +  (~where0) .* ( ...
      (where1) * cpx1 + (~where1) * cpx2  );
  cpy = where0 .* cpy0  +  (~where0) .* ( ...
      (where1) * cpy1 + (~where1) * cpy2  );
  bdy = where0 * 0  +  (~where0) .* ( ...
      (where1) * 1 + (~where1) * 2  );

  dist = sqrt( (x-cpx).^2 + (y-cpy).^2 );

  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
