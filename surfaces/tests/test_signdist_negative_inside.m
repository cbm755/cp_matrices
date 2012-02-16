function [pass, str] = test_signdist_negative_inside()
  str = ['some cpFunctions return signed distance: check its negative inside (at origin anyway)'];

  pass = [];
  c = 0;

  %% start with 2D
  list = {@cpCircle @cpEllipse @cpSquare @cpRoundedSquare @cpPolygon};

  for i = 1:length(list)
    f = list{i};
    [cpx, cpy, sd] = f(0, 0);
    c = c + 1;
    if sd < 0
      pass(c) = 1;
    else
      disp('  this one failed:');
      sd
      f
      pass(c) = 0;
    end
  end

  %% now do 3D

  list = {@cpSphere};

  for i = 1:length(list)
    f = list{i};
    [cpx, cpy, cpz, sd] = f(0, 0, 0);
    c = c + 1;
    if sd < 0
      pass(c) = 1;
    else
      disp('  Fail:');
      f
      sd
      pass(c) = 0;
    end
  end
