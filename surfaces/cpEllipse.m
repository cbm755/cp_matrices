function [cpx, cpy, dist] = cpEllipse(x, y, a, b, cen)
%CPELLIPSE  Closest Point function for an ellipse.
%   [cpx, cpy, dist] = cpEllipse(x, y, a, b)
%      An ellipse centered at the origin with major axis a and
%      minor axis b.
%   [cpx, cpy, dist] = cpEllipse(x, y, a, b, CEN)
%      An ellipse centered at CEN with major axis a and minor axis b.
%
% Internally, uses cpParamCurve with Newton solves to find cp

%%
  % defaults
  if (nargin < 4)
    if (nargin == 3)
      error('must provide both or either of a,b');
    end
    a = 1.5;
    b = 0.75;
  end
  if (nargin < 5)
    cen = [0,0];
  end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);

% parameterised curve:
xs = @(t) a*cos(t);
ys = @(t) b*sin(t);

% derivative of parametrisation:
xp = @(t) -a*sin(t);
yp = @(t) b*cos(t);

% second derivative:
xpp = @(t) -a*cos(t);
ypp = @(t) -b*sin(t);

% start and endpts of parameter variable
endpt1 = 0;
endpt2 = 2*pi;

%%
x1d = x(:); y1d = y(:);

nx = length(x1d);     % number of points

% find closest points
cpx = zeros(nx,1);
cpy = zeros(nx,1);
dist = zeros(nx,1);
fail = zeros(nx,1);

for pt = 1:nx
  [cpx(pt), cpy(pt), dist(pt), fail(pt)] = ...
      cpParamCurveClosed(x1d(pt),y1d(pt),xs,ys,xp,yp,xpp,ypp, [endpt1 endpt2]);
  %[t1,t2,t3,t4] = cpParamCurve_2D(x1d(pt),y1d(pt),xs,ys,xp,yp,xpp,ypp,endpt1,endpt2,1);
  %cpx(pt) = t1;
  %cpy(pt) = t2;
  %dist(pt) = t3;
  %fail(pt) = t4;
end

cpx = reshape(cpx, size(x));
cpy = reshape(cpy, size(x));
dist = reshape(dist, size(x));
fail = reshape(fail, size(x));

figure(1);
plot(cpx(~fail), cpy(~fail), 'bx');
hold on;
axis equal;
plot(x(~~fail), y(~~fail), 'rx');

  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
