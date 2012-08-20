
dx = 0.25;

% TODO: dx = 0.14 is minimum on my machine, need to use refine grid

% make vectors of x, y, positions of the grid
x1d = (-2.0:dx:2.0)';
%y1d = x1d;
%z1d = x1d;

X1d = {x1d x1d x1d x1d x1d};

nx = length(x1d);

tic
[xx{1:5}] = ndgrid(X1d{:});
toc

% closest point to (2D) sphere in 3D
cen = [0  0  0  rand*0.1  rand*0.1];
cpf = @(x) cpSphereInHighDim(x, 1, cen);

tic
[cpX,dist] = cpf(xx);
toc
%tic
%[cpx,cpy,cpz] = cpSphere(xx,yy,zz, 0.75);
%toc
%W1C = dx/3;
%W2C = dx/4;
%cpw1 = W1C*ones(size(cpx));
%cpw2 = W2C*ones(size(cpx));
%dist = sqrt( (xx-cpx).^2 + (yy-cpy).^2 + (zz-cpz).^2 + ...
%             (ww1-cpw1).^2 + (ww2-cpw2).^2);


% banding
dim = 5;
p = 3;       % max interpolation order
stenrad = 1; % max stencil radius for finite differences
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((stenrad+(p+1)/2)^2));
band = find(abs(dist) <= bw*dx);

% TODO: just put in a cpgrid and call a banding routine
%cpgrid1.x1d = X1d;
%cpgrid1.band = find(dist > -100);
%cpgrid1.x = xx;
%cpgrid1.cpx = cpX;
%cpgrid1.dist = dist;
%cpgrid1.cpfun = cpf;
%cpgrid1.dx = dx;

% store closest points in the band;
for d = 1:dim
  cpX{d} = cpX{d}(band);
  xg{d} = xx{d}(band);
end
dist = dist(band);

cpgrid1.x1d = X1d;
cpgrid1.band = band;
cpgrid1.x = xg;
cpgrid1.cpx = cpX;
cpgrid1.dist = dist;
cpgrid1.cpfun = cpf;
cpgrid1.dx = dx;


disp('refining once');
cpgrid2 = refine_gridnd(cpgrid1, bw);

%disp('refining again');
cpgrid3 = refine_gridnd(cpgrid2, bw);

cpgrid = cpgrid2;

stop

% TODO: make the matrix routinues support cell array
cpXtemp = [cpgrid.cpx{1:5}];

tic
L = laplacian_nd_matrix(cpgrid.x1d, 2, cpgrid.band);
toc

tic
%E1 = interpn_matrix(X1d, cpX, 1, band);
[Ei,Ej,Es] = interpn_matrix(cpgrid.x1d, cpXtemp, 1, cpgrid.band);
toc

tic
E1 = sparse(Ei, Ej, Es, length(cpgrid.band), length(cpgrid.band));
toc

tic
[Ei,Ej,Es] = interpn_matrix(cpgrid.x1d, cpXtemp, 2, cpgrid.band);
toc

tic
E2 = sparse(Ei, Ej, Es, length(cpgrid.band), length(cpgrid.band));
toc


%E3 = interpn_matrix(X1d, cpX, p, band);


stop

%% Construct an interpolation matrix for plotting on circle

% plotting grid on circle, using theta as a parameterization
thetas = linspace(0,2*pi,100)';
r = ones(size(thetas));
% plotting grid in Cartesian coords
[xp,yp] = pol2cart(thetas,r);
xp = xp(:);  yp = yp(:);
zp = ZC*ones(size(xp));
wp = WC*ones(size(xp));

Eplot = interpn_matrix({x1d, y1d, z1d, w1d}, [xp yp zp wp], p, band);

figure(1); clf;
figure(2); clf;
figure(3); clf;


%% Time-stepping for the heat equation
Tf = 10;
%dt = 0.2*dx^2;
dt = 0.5*dx;
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps

u = u0;


lambda = -2*dim/dx^2;
I = speye(size(L));
M1 = E1*L - lambda*(I - E3);
A = I - dt*M1;


for kt = 1:numtimesteps
  % explicit Euler timestepping
  %unew = u + dt*(L*u);
  %u = E3*unew;

  % MOL, penalty method
  unew = A \ u;
  u = unew;

  % implicit Euler with bilinear interp p=1
  %unew = u + dt*((E1*L)*unew);
  %unew = A \ u;

  % closest point extension with p=3
  %u = E1*unew;
  %u = unew;

  t = kt*dt;

  if ( (kt < 5) || (mod(kt,numtimesteps/20) == 0) || (kt == numtimesteps) )
    % plot value on circle
    set(0, 'CurrentFigure', 2);
    clf;
    circplot = Eplot*u;
    exactplot = exp(-t)*cos(thetas);
    plot(thetas, circplot);
    title( ['soln at time ' num2str(t) ', on circle'] );
    xlabel('theta'); ylabel('u');
    hold on;
    % plot analytic result
    plot(thetas, exactplot, 'r--');
    legend('numerical soln', 'exact soln', 'Location', 'SouthEast');
    error_circ_inf = max(abs( exactplot - circplot ))

    set(0, 'CurrentFigure', 3);
    clf;
    plot(thetas, circplot - exactplot);
    title( ['error at time ' num2str(t) ', on circle'] );
    xlabel('theta'); ylabel('error');

    %pause(0);
    drawnow();
  end
end
