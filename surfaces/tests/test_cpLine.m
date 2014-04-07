function [pass, str] = test_cpLine()
  str = ['cpLine: various simple 2d and 3d tests'];

  c = 0; pass = [];

  c=c+1;  pass(c) = ~isnan(cpLine(0,0,[1 0],[0 0]));

  [x,y,d,b,s] = cpLine(0,0,[1 0],[0 0]);
  c=c+1;  pass(c) = assertAlmostEqual([x,y,d,b,s],[0,0,0,0,0]);

  [x,y,d,b,s] = cpLine(-1,0,[1 0],[0 0]);
  c=c+1;  pass(c) = assertAlmostEqual([x,y,d,b,s],[-1,0,0,0,-1]);

  [x,y,d,b,s] = cpLine(10,0,[5 0],[0 0]);
  c=c+1;  pass(c) = assertAlmostEqual([x,y,d,b,s],[10,0,0,0,2]);

  %% 3D
  cpf = @(x,y,z) cpLine(x,y,z, [1 0 0], [0 0 0]);
  [x,y,z,d,b,s] = cpf(0,0,0);
  c=c+1;  pass(c) = assertAlmostEqual([x,y,z,d,b,s],[0,0,0,0,0,0]);
  [x,y,z,d,b,s] = cpf(-1,0,0);
  c=c+1;  pass(c) = assertAlmostEqual([x,y,z,d,b,s],[-1,0,0,0,0,-1]);
  [x,y,z,d,b,s] = cpf(10,0,0);
  c=c+1;  pass(c) = assertAlmostEqual([x,y,z,d,b,s],[10,0,0,0,0,10]);
  [x,y,z,d,b,s] = cpf(-1,2,0);
  c=c+1;  pass(c) = assertAlmostEqual([x,y,z,d,b,s],[-1,0,0,2,0,-1]);
  [x,y,z,d,b,s] = cpf(-1,0,3);
  c=c+1;  pass(c) = assertAlmostEqual([x,y,z,d,b,s],[-1,0,0,3,0,-1]);

