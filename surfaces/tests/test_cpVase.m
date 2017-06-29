function [pass, str] = test_cpVase()
  str = 'tests for cpVase';

  pass = [];
  c = 0;
  tol = 3*eps;

  makeplots = 0;

  x1d = linspace(-2,2,15);

  [xx,yy,zz] = meshgrid(x1d,x1d,x1d);

  %lim = [0.15*pi  0.4*pi];
  lim = [];
  [cpx,cpy,cpz,dist,bdy] = cpVase(xx,yy,zz, lim);

  % dist should match computed dist
  dist2 = sqrt((cpx-xx).^2 + (cpy-yy).^2 + (cpz-zz).^2);
  c = c + 1;
  pass(c) = assertAlmostEqual(dist,dist2,tol);

  % min/max matches param function
  [x,y,z] = paramVase(32, lim);
  c = c + 1;
  pass(c) = assertAlmostEqual(max(max(max(cpz))), max(max(z)));
  c = c + 1;
  pass(c) = assertAlmostEqual(min(min(min(cpz))), min(min(z)));

  if makeplots
    figure(1); clf;
    %param = @(N) paramVase(N, lim);
    %porcupine_plot3d(xx,yy,zz, cpx,cpy,cpz, bdy, param);
    surf(x,y,z);
    hold on;
    I = bdy == 0;
    plot3(cpx(I), cpy(I), cpz(I), 'k.');
    hold on;
    axis equal
    I = bdy == 1;
    plot3(cpx(I), cpy(I), cpz(I), 'r.');
    I = bdy == 2;
    plot3(cpx(I), cpy(I), cpz(I), 'b.');
  end
