function [cpx, cpy, dist, exterior, param, nx, ny] = cpBeanInterior(x,y,scale)

  if nargin<3
      scale = 1.0;
  end
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

  sp = cscvn(pts');

  DEBUG = 0;
  [cpx, cpy, dist, param] = cpSpline2D(x, y, sp, 1, DEBUG);

  
  sp1 = fnder(sp);
  S1.type='()';
  S1.subs = {1,':'};
  S2.type='()';
  S2.subs = {2,':'};
  % derivative of parametrisation:
  xp = @(t) subsref(ppval(sp1,t), S1);
  yp = @(t) subsref(ppval(sp1,t), S2);
  
  % compute unit tangent vectors
  tangent_lengths = sqrt(xp(param).^2 + yp(param).^2);
  tx = xp(param)./tangent_lengths;
  ty = yp(param)./tangent_lengths;
  
  % compute unit outward normal vectors
  nx = ty';
  ny = -tx';
  
  % compute normal vectors with lengths to be dist
  Nx = x - cpx;
  Ny = y - cpy;
  
  % If (nx,ny) and (Nx,Ny) have the same direction, then (x,y) is outside
  % the beancurve.
  exterior = abs(Nx - nx.*dist) + abs(Ny - ny.*dist) < 1e-10*dist & dist > 1e-12;
  interior = ~exterior;
  cpx(interior) = x(interior);
  cpy(interior) = y(interior);
  
  dist(interior) = 0;
  
end