function [nn,tt1,tt2] = normals_from_diffcp3d(g)
%NORMALS_FROM_DIFFCP3D
%   Helper function for 3D

  % lots of testing/debugging code here, these enable it
  make_tests = 0; make_plots = 0;

  [Dx,Dy,Dz] = firstderiv_matrices(g, 'centered');
  E = interp_matrix(g, 3);

  Jxx = Dx * g.cpx;
  Jxy = Dy * g.cpx;
  Jxz = Dz * g.cpx;

  Jyx = Dx * g.cpy;
  Jyy = Dy * g.cpy;
  Jyz = Dz * g.cpy;

  Jzx = Dx * g.cpz;
  Jzy = Dy * g.cpz;
  Jzz = Dz * g.cpz;

  Jxx = E*Jxx;
  Jxy = E*Jxy;
  Jxz = E*Jxz;

  Jyx = E*Jyx;
  Jyy = E*Jyy;
  Jyz = E*Jyz;

  Jzx = E*Jzx;
  Jzy = E*Jzy;
  Jzz = E*Jzz;

  if (make_tests)
    dx = g.x1d(2)-g.x1d(1);
  end
  if (make_tests && make_plots)
    figure(1); clf;
    sphere(32);
    %plot2d_compdomain(g.x,g.x,g.y,g.dx,g.dx, 1)
    %thp = linspace(0,2*pi,1000);
    %plot(cos(thp), sin(thp), 'k-');
    axis equal
    hold on
  end

  tic

  % codim-1 assumed here
  nn = zeros(length(g.band), 3);
  tt1 = zeros(length(g.band), 3);
  tt2 = zeros(length(g.band), 3);

  data = [];
  for i=1:length(g.band);
    if (1==1)
    %if abs(g.dist(i)) < 1*dx
      A = [Jxx(i) Jxy(i) Jxz(i); ...
           Jyx(i) Jyy(i) Jyz(i); ...
           Jzx(i) Jzy(i) Jzz(i)];
      A = (A + A')/2;
      [V,D] = eig(A);
      ews = diag(D);
      [ews,I] = sort(ews);
      V = V(:,I);

      % normal is first eigenvector, next two span the tangent
      % plane (in codim-1)
      n = V(:,1);
      t1 = V(:,2);
      t2 = V(:,3);
      nn(i,:) = [n(1) n(2) n(3)];
      tt1(i,:) = [t1(1) t1(2) t1(3)];
      tt2(i,:) = [t2(1) t2(2) t2(3)];

      n2o = [g.cpx(i) - g.x(i); g.cpy(i) - g.y(i); g.cpz(i) - g.z(i)];
      n2 = n2o ./ norm(n2o,2);

      if (make_tests)
        nex = [g.x(i); g.y(i); g.z(i)];
        nex = nex ./ norm(nex,2);
        A = [norm(n2o,2) g.dist(i)  ews' ...
             abs(nex'*n2)-1  abs(nex'*n)-1  ...
             abs(nex'*t1)  abs(nex'*t2)];
        disp(A);
        data = [data; A];
        if (make_plots)
          %n2 = mynormr(n2);
          plot3([g.x(i) g.x(i)+n(1)/10], ...
                [g.y(i) g.y(i)+n(2)/10], ...
                [g.z(i) g.z(i)+n(3)/10], 'b-');
          plot3([g.x(i) g.x(i)+n2(1)/10], ...
                [g.y(i) g.y(i)+n2(2)/10], ...
                [g.z(i) g.z(i)+n2(3)/10], 'r-');
          %drawnow()
        end
      end
    end
  end
  toc

  if (make_tests)
    disp(sprintf([...
        'h=%g \t ew1 in [%.3g,%.3g] \t ew in [%.3g,%.3g]' ...
        ' \t n''*nex: %.3g \t t''*nex: %.3g'],  ...
                 g.dx, min(data(:,3)), max(data(:,3)), ...
                 min(min(data(:,4:5)))-1, max(max(data(:,4:5)))-1, max(abs(data(:,7))), ...
                 max(max(abs(data(:,8:9)))) ...
                 ));
  end