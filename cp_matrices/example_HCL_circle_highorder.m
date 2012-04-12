function example_HCL_circle_highorder()
%% HCL on a circle
% This example solves the Burgers' equation on a 2D circle, using a
% specified velocity field tangent to the circle.
% It needs to be a function so we can use a subroutine

%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../cp_matrices');

% add functions for finding the closest points
addpath('../surfaces');



%%
% 2D example on a circle
% Construct a grid in the embedding space

dx = 0.05;                   % grid size

% make vectors of x, y, positions of the grid
x1d = (-2.0:dx:2.0)';
y1d = x1d;

nx = length(x1d);
ny = length(y1d);



%% Find closest points on the surface
% For each point (x,y), we store the closest point on the circle
% (cpx,cpy)

% meshgrid is only needed for finding the closest points, not afterwards
[xx yy] = meshgrid(x1d, y1d);
% function cpCircle for finding the closest points on a circle
[cpx, cpy, dist] = cpCircle(xx,yy);
% make into vectors
cpxg = cpx(:); cpyg = cpy(:);


%% Banding: do calculation in a narrow band around the circle
dim = 2;
p = 5;          % interp degree (5 for WENO6 interp, 3 for WENO4 interp)
opStenRad = 3;  % the stencil radius of the spatial disc scheme

% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
% For weno, see the figure in [Macdonald & Ruuth 2008]
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((opStenRad+(p+1)/2)^2));
band = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpxg = cpxg(band); cpyg = cpyg(band);
xg = xx(band); yg = yy(band);


%% Function u in the embedding space
% u is a function defined on the grid
[th, r] = cart2pol(xg,yg);
u = 0.45*cos(th) + 0.55;
%u = cos(th);
initialu = u;       % store initial value

w1 = -sin(th);
w2 = cos(th);




%% Construct an interpolation matrix for closest point
% This creates a matrix which interpolates data from the grid x1d y1d,
% onto the points cpx cpy.

disp('Constructing interpolation and laplacian matrices');

E = interp2_matrix(x1d, y1d, cpxg, cpyg, p, band);

% e.g., closest point extension:
%u = E*u;

%% Create differentiation matrices

order = 2;  % Laplacian order: bw will need to increase if changed
L = laplacian_2d_matrix(x1d,y1d, order, band);
[Dxb,Dxf, Dyb,Dyf] = firstderiv_upw1_2d_matrices(x1d,y1d, band);

%% Construct an interpolation matrix for plotting on circle

% plotting grid on circle, using theta as a parameterization
thetas = linspace(0,2*pi,512)';
r = ones(size(thetas));
% plotting grid in Cartesian coords
[xp,yp] = pol2cart(thetas,r);
xp = xp(:); yp = yp(:);
Eplot = interp2_matrix(x1d, y1d, xp, yp, p, band);


figure(1); clf;
figure(2); clf;

cp.dim = 2;
cp.x1d = x1d;
cp.y1d = y1d;
cp.band = band;
cp.cpx = cpxg;
cp.cpy = cpyg;
% notation: Neighbour East = "NE", etc
[NE NW NN NS] = neighbourMatrices(cp, cp.band, cp.band);
% weno interpolation can cache data for better performance:
Weno6Cache = weno6_interp(cp, u, [cp.cpx cp.cpy], true);

%% Time-stepping for the heat equation

Tf = 10;
dt = 0.1*dx;  % TODO: what is CFL for this?
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps

for kt = 1:numtimesteps

  % SSPRK(3,3) stage 1
  unew = u + dt*F(u);
  %u1 = E*unew;
  u1 = weno6_interp(Weno6Cache, unew);

  % stage 2
  unew = 0.75*u + 0.25*u1 + (0.25*dt)*F(u1);
  u2 = weno6_interp(Weno6Cache, unew);

  % stage 3
  unew = (1/3)*u + (2/3)*u2 + ((2/3)*dt)*F(u2);
  u = weno6_interp(Weno6Cache, unew);

  t = kt*dt;

  if ( (kt < 10) || (mod(kt,10) == 0) || (kt == numtimesteps) )
    % plot over computation band
    plot2d_compdomain(u, xg, yg, dx, dx, 1)
    title( ['embedded domain: soln at time ' num2str(t) ...
            ', timestep #' num2str(kt)] );
    xlabel('x'); ylabel('y');
    %hold on
    plot(xp,yp,'k-', 'linewidth', 2);
    %axis equal;  axis tight;

    % plot value on circle
    set(0, 'CurrentFigure', 2);
    clf;
    circplot = Eplot*u;
    plot(thetas, circplot, 'b-+');
    title( ['soln at time ' num2str(t) ', on circle'] );
    xlabel('theta'); ylabel('u');
    hold on;
    grid on
    drawnow();
  end
end



function rhs = F(u)
  scalar_flux = u.^2/2;  % burgers
  f = scalar_flux.*w1;  % find the components
  g = scalar_flux.*w2;

  % Lax-Friedrichs flux splitting
  % TODO: here f' = u so say 1 for now but should search, GLF, LLF, SLLF
  alphax = 1;
  alphay = 1;
  % split intp f^{+} and f^{-}
  fp = 1/2*(f + alphax*u);
  fm = 1/2*(f - alphax*u);
  gp = 1/2*(g + alphay*u);
  gm = 1/2*(g - alphay*u);


  %% Reconstruction of cell-edge values
  % notation: f^{p}_{i+1/2} = fp_iph

  % HCL WENO5
  fp_iph = fd_weno5_1d( [NW*(NW*fp)  NW*fp  fp  NE*fp  NE*(NE*fp)] );
  fm_iph = fd_weno5_1d( [NE*(NE*(NE*fm))  NE*(NE*fm)  NE*fm  fm  NW*fm] );
  gp_jph = fd_weno5_1d( [NS*(NS*gp)  NS*gp  gp  NN*gp  NN*(NN*gp)] );
  gm_jph = fd_weno5_1d( [NN*(NN*(NN*gm))  NN*(NN*gm)  NN*gm  gm  NS*gm] );
  % first-order upwinding
  %fp_iph = fp;
  %fm_iph = NE*fm;
  %gp_jph = gp;
  %gm_jph = NN*gm;

  % numerical flux differencing
  rhs = - ( ...
      Dxb*fp_iph + Dxb*fm_iph + ...
      Dyb*gp_jph + Dyb*gm_jph ...
      );

end
end