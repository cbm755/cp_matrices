function [nn,tt] = normals_from_diffcp2d(g)
%NORMALS_FROM_DIFFCP2D
%   Helper function for 2D

  make_tests = 0; make_plots = 0;

  [Dx,Dy] = firstderiv_matrices(g, 'centered');
  E = interp_matrix(g, 3);

  Jxx = E * (Dx * g.cpx);
  Jxy = E * (Dy * g.cpx);
  Jyx = E * (Dx * g.cpy);
  Jyy = E * (Dy * g.cpy);

  if (make_plots)
    [th,r] = cart2pol(g.x,g.y);
    figure(1); clf;
    %plot2d_compdomain(g.x,g.x,g.y,g.dx,g.dx, 1)
    thp = linspace(0,2*pi,1000);
    plot(cos(thp), sin(thp), 'k-');
    axis equal
    hold on
  end

  nn = zeros(length(g.band), 2);
  tt = zeros(length(g.band), 2);

  for i=1:length(g.band);
    A = [Jxx(i)  Jxy(i);  Jyx(i)  Jyy(i)];
    [V,D] = eig(A);
    ews = diag(D);
    [ews,I] = sort(ews);
    V = V(:,I);

    % normal is first eigenvector, next two span the tangent
    % plane (in codim-1)
    n = V(:,1);
    t = V(:,2);
    nn(i,:) = [n(1) n(2)];
    tt(i,:) = [t(1) t(2)];

    if (make_tests)
      tex = [-sin(th(i)); cos(th(i))];
      nex = [cos(th(i)); sin(th(i))];
      n2 = [g.cpx(i) - g.x(i); g.cpy(i) - g.y(i)];
      n3 = n2 ./ norm(n2,2);
      disp([norm(n2,2) g.dist(i)  t'*n2  t'*n3  t'*n])
      %n2 = mynormr(n2);
      plot([g.x(i) g.x(i)+n3(1)/10], [g.y(i) g.y(i)+n3(2)/10], 'r-');
      plot([g.x(i) g.x(i)+n(1)/10],  [g.y(i) g.y(i)+n(2)/10], 'b-');
      pause(0)
    end
  end


