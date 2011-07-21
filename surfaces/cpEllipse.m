function [cpx, cpy, dist] = cpEllipse(x, y, a, b)
%CPELLIPSE  Closest Point function for an ellipse.
%   [cpx, cpy, dist] = cpEllipse(x, y)
%      An ellipse centered at the origin with major axis a and minor axis b.

% Uses cpParamCurve with Newton solves to find cp

%%

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
      cpParamCurve_2D(x1d(pt),y1d(pt),xs,ys,xp,yp,xpp,ypp,endpt1,endpt2,1);
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
