function [pass, str] = test_cpTorus()
  str = 'cpTorus tests, test a bunch of points';

  pass = [];
  c = 0;

  [cpx,cpy,cpz,sdist] = cpTorus(1.4,0,0);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx,cpy,cpz,sdist], [1.4,0,0,0]);

  [cpx,cpy,cpz,sdist] = cpTorus(-1,0,0);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx,cpy,cpz,sdist], [-1.4,0,0,-0.4]);

  [cpx,cpy,cpz,sdist] = cpTorus(1,0,0);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx,cpy,cpz,sdist], [1.4,0,0,-0.4]);

  [cpx,cpy,cpz,sdist] = cpTorus(0,1,0);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx,cpy,cpz,sdist], [0,1.4,0,-0.4]);

  [cpx,cpy,cpz,sdist] = cpTorus(1,0,0.4);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx,cpy,cpz,sdist], [1,0,0.4,0]);

  % from origin, check cpz is zero and distance is correct
  [cpx,cpy,cpz,sdist] = cpTorus(0,0,0);
  c = c + 1;
  pass(c) = assertAlmostEqual([sqrt(cpx^2+cpy^2),cpz,sdist], [0.6,0,0.6]);

  % from inner medial axis, check distance
  th = linspace(0,2*pi,50);
  [x,y] = pol2cart(th, 1.5);
  z = zeros(size(x));
  [cpx,cpy,cpz,sdist] = cpTorus(x,y,z, 1.5, 0.5);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpz,sdist], ...
                              [zeros(size(cpx)),-0.5*ones(size(sdist))]);

  % check the parameters
  [cpx,cpy,cpz,sdist] = cpTorus(2,0,0,10,2);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx,cpy,cpz,sdist], [8,0,0,6]);

  % check the center shift
  [cpx,cpy,cpz,sdist] = cpTorus(2+1,0+2,0+3,10,2,[1 2 3]);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx,cpy,cpz,sdist], [8+1,0+2,0+3,6]);

  % compute with a meshgrid
  x1d = linspace(-2,2,8);
  [x,y,z] = meshgrid(x1d,x1d,x1d);
  [cpx,cpy,cpz,sdist] = cpTorus(x,y,z, 1.1, 0.3, [1 2 3]);
  % check the distance calc
  mydist = sqrt( (x-cpx).^2 + (y-cpy).^2 + (z-cpz).^2 );
  c = c + 1;
  pass(c) = assertAlmostEqual(mydist, abs(sdist));

  % point in the parameterization all have distance zero
  [x,y,z] = paramTorus(16);
  [cpx,cpy,cpz,sdist] = cpTorus(x,y,z);
  c = c + 1;
  pass(c) = assertAlmostEqual(sdist, 0);

