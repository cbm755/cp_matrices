%% Non-Euclidean Closest Points
%
% The Closest Point Method works with more general nototions of
% distance.  Here we look at an ellipse in both Euclidean and
% non-Euclidean closest points.
%
%% An ellipse

% first, an ellipse
a = 1.6;
b = 0.8;
figure(1); clf;
[xp, yp, thp] = paramEllipse(1000,a,b);
lw = 'linewidth';
plot(xp,yp,'k-', lw, 1);
title('Ellipse');  xlabel('x');  ylabel('y');
xlim([-1.7 1.7])
axis equal
hold on

% closest point functions
ecp = @(x,y) cpEllipse(x,y,a,b);
cp = @(x,y) cpEllipse_nonEuclidean(x,y,a,b);


%% Distance & Closest Points

% make a closest point representation on a meshgrid
g = make_coarse_grid(2, cp, 0.1, [-2 -1 2 1]);
ge = make_coarse_grid(2, ecp, 0.1, [-2 -1 2 1]);

figure(2); clf
plot2d_compdomain(g.dist, g, 2)
plot(xp,yp,'k-',lw,2)
title('non-Euclidean distance')

figure(3); clf
plot2d_compdomain(ge.dist, ge, 3)
plot(xp,yp,'k-',lw,2)
title('Euclidean distance')

%%
% Note the difference in interior of the ellipse: particularly the
% medial axes where trajectories of curves with the same closest
% point collide:

figure(2); clf;
porcupine_plot2d(g.x, g.y, g.cpx, g.cpy)
axis tight
title('non-Euclidean closest points')

figure(3); clf;
porcupine_plot2d(ge.x, ge.y, ge.cpx, ge.cpy)
axis tight
title('Euclidean closest points')






%% Points on a trajectory
% For Euclidean closest points, all points in the normal direction
% to the surface share the same closest point.  For the particular
% non-Euclidean distance used here, a curve sharing a common
% closest point which passes through the points (X,Y) and (x,y) is
% given by
%
% $$ X^{a^2} y^{b^2} = Y^{b^2} x^{a^2}. $$

% choose a point
X = 1.3; Y = 1.2;

% plot the steepest descent curve from that point.  All points on this
% curve share the same closest point.
x = linspace(0, 2, 256);
y = ((Y^(b^2) * x.^(a^2)) / X^(a^2)).^(1/b^2);

figure(1);
plot(x, y, 'm-')
plot(X, Y, 'm*')

% our computed closest point should agree
[cpx, cpy, sd, t] = cp(X,Y);
plot(cpx, cpy, 'mo');
title('Trajectory curve sharing same closest point');
ylim([-1.5 1.5])



%%
% We know points on the surface are given by $x(t) = a\cos(t)$ and
% $y(t) = a\sin(t)$ this gives an equation to solve for $t$ to find
% the closest point for a given $(X,Y)$.  Our cp function uses
% bisection search to find such a $t$.

% add a few more curves
for i=0:10
  X = -2 + i/5; Y = 0 + i/7;
  if i == 10
    y = linspace(0,2,256);
    x = zeros(size(y));
  else
    x = linspace(-2,0,256);
    y = ((Y^(b^2) * (-x).^(a^2)) / (-X)^(a^2)).^(1/b^2);
  end
  col = [i/10 0 1-i/10];
  plot(x, y, '-', 'color', col)
  plot(X, Y, '*', 'color', col)
  [cpx, cpy, sd, t] = cp(X,Y);
  plot(cpx, cpy, 'o', 'color', col);
end
title('Trajectory curves sharing closest points');


%% Closest Point Method calculations: the diffusion equation
% We solve the diffusion equation on the ellipse
%
% $$u_t = \triangle_s u$$

% As usual, we can use a narrow band
g = refine_cpgrid(g, 5);
ge = refine_cpgrid(ge, 5);
% optionally refine a few more times (to do a convergence study)
%g = refine_cpgrid(g, 5);
%ge = refine_cpgrid(ge, 5);
%g = refine_cpgrid(g, 5);
%ge = refine_cpgrid(ge, 5);

% we also need a highly resolved grid which I'll use in lieu of
% exact solutions
gr = refine_cpgrid(ge);
gr = refine_cpgrid(gr, 5);


%%

figure(2); clf;
plot2d_compdomain(g.dist, g, 2);
plot(xp, yp, 'k-', lw, 2);

figure(3); clf;
plot2d_compdomain(ge.dist, ge, 3);
plot(xp, yp, 'k-', lw, 2);

%%
% Notice the two grids look the same and indeed:

length(g.x)
length(ge.x)
nnz(g.x - ge.x)

%%
% This is because the grid is the minimal grid to contain a certain
% stencil which does not depend on the distance norm used.  However,
% the closest points are different:

max(abs(g.cpx - ge.cpx))



%% Differential operators
% Build the discerete differential operators (matrices) needed to
% build a Closest Point Method

[Dxb, Dxf, Dyb, Dyf] = firstderiv_upw1_2d_matrices(g.x1d, g.y1d, g.band);
[Dxc, Dyc] = firstderiv_cen2_2d_matrices(g.x1d, g.y1d, g.band);

L  = laplacian_2d_matrix(g.x1d,  g.y1d,  2, g.band);
Le = laplacian_2d_matrix(ge.x1d, ge.y1d, 2, ge.band);
Lr = laplacian_2d_matrix(gr.x1d, gr.y1d, 2, gr.band);

E  = interp2_matrix(g.x1d,  g.y1d,  g.cpx,  g.cpy,  3, g.band);
Ee = interp2_matrix(ge.x1d, ge.y1d, ge.cpx, ge.cpy, 3, ge.band);
Er = interp2_matrix(gr.x1d, gr.y1d, gr.cpx, gr.cpy, 3, gr.band);

size(E)
size(Er)



%% Laplace-Beltrami consistency
% For a Euclidean closest point function, we have a theorem
% for the Laplacian which states
%
% $$\triangle_S u = \triangle (Eu),$$
%
% for points on the surface [MM 2012].  That is, one can compute the
% standard (cartesian) Laplacian of an extension and it will agree (on
% the surface) with the Laplace-Beltrami operator.
%
% For non-Euclidean closest points, the above is not true.  Instead we
% use the gradient and divergence principles (which are true for
% general closest point functions).  This gives
%
% $$\triangle_S u = \textrm{div} (E \textrm{grad} (Eu)),$$
%
% for points on the surface [RM 2008], [MM 2012].



%% Time-stepping
% Here we use the method-of-lines penalty approach of [vGMM 2013]
% but could also use the Ruuth--Merriman iteration.

Tf = 0.1;
dt = 0.25*g.dx^2;
numsteps = ceil(Tf/dt)
lambda = 4/(g.dx^2);

u = sin(g.cpx);
ue = sin(ge.cpx);

for k=1:numsteps
  % this would be inconsistent---try it!
  %u = u + dt*(E*(L*u)) - dt*lambda*(u - E*u);

  % 1st-order differences
  %ux = E*(Dxf*u);
  %uy = E*(Dyf*u);
  %u = u + dt*(E*(Dxb*ux + Dyb*uy)) - dt*lambda*(u - E*u);

  % 2nd-order
  ux = E*(Dxc*u);
  uy = E*(Dyc*u);
  u = u + dt*(E*(Dxc*ux + Dyc*uy)) - dt*lambda*(u - E*u);

  % Using Euclidean CP
  ue = ue + dt*(Ee*(Le*ue)) - dt*lambda*(ue - Ee*ue);

  if mod(k, 50) == 0
    fprintf('step %d of %d, max=(%g,%g)\n', ...
            k, numsteps, max(abs(u)), max(abs(ue)));
  end
end



%% Reference soln
% Could not compute with other above because it needs a different
% time grid (due to smaller dx).
dt = 0.25*gr.dx^2;
numsteps = ceil(Tf/dt)
lambda = 4/(gr.dx^2);

ur = sin(gr.cpx);
for k=1:numsteps
  ur = ur + dt*(Er*(Lr*ur)) - dt*lambda*(ur - Er*ur);
  if mod(k, 500) == 0
    fprintf('step %d of %d, max=%g\n', k, numsteps, max(abs(ur)));
  end
end



%% Compare on parameterized grid
% We can interpolate back to a set of parameterized points on the
% surface, and compute errors.
% Note xp,yp from linspace on $\theta$: but this is not equispaced
% in arclength.

% Interpolation matrices at the (xp,yp) surface points
Eplot = interp2_matrix(g.x1d, g.y1d, xp, yp, 3, g.band);
Eplote = interp2_matrix(ge.x1d, ge.y1d, xp, yp, 3, ge.band);
Eplotr = interp2_matrix(gr.x1d, gr.y1d, xp, yp, 3, gr.band);

up = Eplot*u;
uep = Eplote*ue;
urp = Eplotr*ur;

figure(2); clf;
plot(thp, up, 'r-');
hold on;
plot(thp, uep, 'b--');
plot(thp, urp, 'k:');
xlim([0 2*pi])
title('solutions');  xlabel('\theta');  ylabel('u');

figure(3); clf;
plot(thp, up-urp, 'r-');
hold on;
plot(thp, uep-urp, 'b--');
xlim([0 2*pi])
err1 = max(abs(up-urp))
err2 = max(abs(uep-urp))
legend(sprintf('euclidean, max=%g',err2), ...
       sprintf('non-eucl, max=%g',err1));
xlabel('\theta');  ylabel('max error');



%% Convergence study
% By running this code for increasing refinements, we can produce
% the following convergence study:
%
%              Error 1         Error 2         Error 3         Error 4
%      dx      eucl cp  conv   non-eucl conv   n-e 1st  conv   n-e E*L
%      ---------------------------------------------------------------
%      0.05    6.68e-5         1.97e-4         3.06e-4         0.0163
%      0.025   1.62e-5  2.04   4.78e-5  2.04   1.40e-4  1.13   0.0163
%      0.0125  3.95e-6  2.04   1.14e-5  2.07   6.72e-5  1.06   0.0163
%
% These results fit the theory of [MM 2012]:
%
% * |Error 1| is the Euclidean closest point calculation, error
%   appears second order accurate as expected.
% * |Error 2| uses non-Euclidean closest points with centered
%   differences, error is 2nd-order.
% * |Error 3| uses non-Euclidean closest points with one-sided
%   differences, 1st-order.
% * |Error 4| uses non-Euclidean closest points with E*L
%   discretization, inconsistent.
%
% Additionally, although both are second-order, the Euclidean closest
% points give a smaller error constant.  But if we redo the Euclidean
% closest point test using E*div(E*grad) instead of E*L, the
% difference mostly goes away:
%
%             Error 5         Error 2
%     dx      ecp EDED conv   non-eucl conv
%     -------------------------------------
%     0.05    1.51e-4         1.97e-4
%     0.025   3.66e-5  2.05   4.78e-5  2.04
%     0.0125  8.81e-6  2.06   1.14e-5  2.07



%% References
% <http://people.maths.ox.ac.uk/macdonald/publications.html>
%
% # [MM 2012]  Thomas März and Colin B. Macdonald, Calculus on surfaces
%   with general closest point functions, SIAM J. Numer. Anal. 2012.
% # [RM 2008]  Steven J. Ruuth, Barry Merriman A simple embedding
%   method for solving partial differential equations on surfaces,
%   J. Comput. Phys. 2008.
% # [vGMM 2013]  Ingrid von Glehn, Thomas März, and Colin B. Macdonald.
%   An embedded method-of-lines approach to solving partial
%   differential equations on surfaces. 2013. Submitted.




