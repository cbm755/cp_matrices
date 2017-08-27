%% Parameters
a = 1.6; % length along x
b = 0.5; % length along y
cpf = @(x, y) cpEllipse(x, y, a, b);
paramf = @(N) paramEllipse(N, a, b);

max_vec=[];
sum_vec=[];
dx=0.1; % grid size


x1d = (-2.0:dx:2.0)';
y1d = x1d;

%% Banding: do calculation in a narrow band around the vase
dim = 2;    % dimension
p = 3;      % interpolation degree
order = 2;  % Laplacian order
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

% meshgrid is only needed for finding the closest points, not afterwards
[xx, yy] = meshgrid(x1d, y1d);

[cpx, cpy, dist] = cpf(xx, yy);

% TODO: currently need to band before cutting b/c SA calc depends on bw
band = find(abs(dist) <= bw*dx);

% keep only closest points in the band
cpx = cpx(band); cpy = cpy(band);
x = xx(band); y = yy(band);

theta=[pi/2];
%theta = linspace(-pi/2,pi/2,1);

for j= 1:length(theta)
  xh=a*cos(theta(j));
  yh=b*sin(theta(j));

  len_wanted=2;
  hole_cen = [xh yh];

  [cpx,cpy,dist,bdy] = cpCutHole2d(x, y, cpf, hole_cen, len_wanted, dx, bw);

  %% Extension of closest point
  E = interp2_matrix(x1d, y1d, cpx, cpy, p, band);

  %% bdy cond;
  E(bdy,:) = -E(bdy,:);

  L = laplacian_2d_matrix(x1d,y1d, order, band, band);


  %% For constant D
  d = 1;
  D = -1/d*ones(length(E),1);

  %% solution of u
  % boundary condition

  gamma = 2*dim/(dx^2);
  I = speye(size(E));
  M = E*L - gamma*(I-E);
  u = M \ D;
  %u = gmres(M,D);

  %% plot
  N = 512;
  [xp, yp, thp] = paramf(N);
  Eplot = interp2_matrix(x1d, y1d, xp(:), yp(:), p, band);

  up = Eplot*u;

  figure(1); clf;
  plot(thp, up);
  xlabel('theta');

  figure(2); clf;
  plot2d_compdomain(u,x,y,dx,dx, 2)
end
