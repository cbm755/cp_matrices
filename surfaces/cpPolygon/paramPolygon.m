function [x, y] = paramPolygon(n, poly)
%PARAMPOLYGON  Parameterization of a polygon
%   [x, y] = paramPolygon(n, poly)
%   'n' is the approximate number of points.

  % defaults
  if (nargin < 2)
    poly = makePolyDefault();
  end

  N = size(poly,1);

  % how many points per segment?  first compute the arclenths
  arclens = zeros(1,N);
  for j=1:(N-1)
    arclens(j) = norm(poly(j,:)-poly(j+1,:),2);
  end
  arclens(N) = norm(poly(N,:)-poly(1,:),2);

  m = ceil( n * (arclens/sum(arclens)));

  [x,y] = priv_paramLineSegment(m(1), poly(1,:), poly(2,:));
  for j=2:(N-1)
    [x2,y2] = priv_paramLineSegment(m(j), poly(j,:), poly(j+1,:));
    x = [x; x2(2:end)];
    y = [y; y2(2:end)];
  end
  [x2,y2] = priv_paramLineSegment(m(end), poly(N,:), poly(1,:));
  x = [x; x2(2:end)];
  y = [y; y2(2:end)];


  % private function for the linesegment
  function [x,y] = priv_paramLineSegment(n, p, q)
    x = linspace(p(1), q(1), max(n,2))';
    y = linspace(p(2), q(2), max(n,2))';
  end

end