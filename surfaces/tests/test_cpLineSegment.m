function [pass, str] = test_cpLineSegment()
  str = ['cpLineSegment: 2d and 3d tests'];

  c = 0;  pass = [];

  %% 2D
  cpf = @(x,y) cpLineSegment2d(x,y, [0 0], [1 0]);

  [x,y,d,b,s] = cpf(0.5,0);
  c=c+1;  pass(c) = assertAlmostEqual([x,y,d,b,s],[0.5,0,0,0,0.5]);
  [x,y,d,b,s] = cpf(0,0);
  c=c+1;  pass(c) = assertAlmostEqual([x,y,d,b,s],[0,0,0,0,0]);
  [x,y,d,b,s] = cpf(1,0);
  c=c+1;  pass(c) = assertAlmostEqual([x,y,d,b,s],[1,0,0,0,1]);
  [x,y,d,b,s] = cpf(-1,0);
  c=c+1;  pass(c) = assertAlmostEqual([x,y,d,b,s],[0,0,1,1,0]);
  [x,y,d,b,s] = cpf(-2,0);
  c=c+1;  pass(c) = assertAlmostEqual([x,y,d,b,s],[0,0,2,1,0]);
  [x,y,d,b,s] = cpf(-3,4);
  c=c+1;  pass(c) = assertAlmostEqual([x,y,d,b,s],[0,0,5,1,0]);
  [x,y,d,b,s] = cpf(2,0);
  c=c+1;  pass(c) = assertAlmostEqual([x,y,d,b,s],[1,0,1,2,1]);


  %% 3D
  cpf = @(x,y,z) cpLineSegment3d(x,y,z, [0 0 0], [1 0 0]);
  [x,y,z,d,b,s] = cpf(0,0,0);
  c=c+1;  pass(c) = assertAlmostEqual([x,y,z,d,b,s],[0,0,0,0,0,0]);
  [x,y,z,d,b,s] = cpf(-1,0,0);
  c=c+1;  pass(c) = assertAlmostEqual([x,y,z,d,b,s],[0,0,0,1,1,0]);
  [x,y,z,d,b,s] = cpf(-3,4,0);
  c=c+1;  pass(c) = assertAlmostEqual([x,y,z,d,b,s],[0,0,0,5,1,0]);
  [x,y,z,d,b,s] = cpf(3,0,0);
  c=c+1;  pass(c) = assertAlmostEqual([x,y,z,d,b,s],[1,0,0,2,2,1]);

