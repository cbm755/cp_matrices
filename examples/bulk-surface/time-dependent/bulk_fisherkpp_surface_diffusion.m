%% Parameters
d = 0.001;
D = 40*d;
Nu = 0.1;
Mu = 5;
dx = 0.025;    

%% Construct a grid in the embedding space
% make vectors of x, y, positions of the grid
x1d = (-2.0:dx:2.0)';
y1d = x1d;

nx = length(x1d);
ny = length(y1d);


%% Find closest points on the surface
% For each point (x,y), we store the closest point on the circle
% (cpx,cpy)

% meshgrid is only needed for finding the closest points, not afterwards
[xx,yy] = meshgrid(x1d, y1d);
% Finding closest points for the unit disk.
[cpxV,cpyV,distV,bdyV] = cpDisk(xx,yy);
[cpxbarV,cpybarV] = cpbar_2d(xx,yy,@cpDisk);
% cpxbarV = cpxV; cpybarV = cpyV;
% cpxbarV(bdyV) = 2*cpxV(bdyV) - xx(bdyV);
% cpybarV(bdyV) = 2*cpyV(bdyV) - yy(bdyV);
% Finding closest points for the unit circle.
[cpxU, cpyU, distU] = cpCircle(xx,yy);
% make into vectors
cpxgU = cpxU(:); cpygU = cpyU(:);
cpxgbarV = cpxbarV(:); cpygbarV = cpybarV(:);
cpxgV = cpxV(:); cpygV = cpyV(:);


%% Banding: do calculation in a narrow band around the circle
dim = 2;    % dimension
p = 3;      % interpolation degree
order = 2;  % Laplacian order
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
bandV = find(distV <= bw*dx); % create v band
bandU = find(abs(distU) <= bw*dx); % create u band


% store closest points in the band;
cpxgbarV = cpxgbarV(bandV); cpygbarV = cpygbarV(bandV); distV = distV(bandV); bdyV = bdyV(bandV);
cpxgV = cpxgV(bandV); cpygV = cpygV(bandV);
xgV = xx(bandV); ygV = yy(bandV);

cpxgU = cpxgU(bandU); cpygU = cpygU(bandU); distU = distU(bandU);
xgU = xx(bandU); ygU = yy(bandU);


%% Set up initial conditions
[thV, rV] = cart2pol(xgV,ygV);
[thU,rU] = cart2pol(xgU,ygU);
% v = (1 - (rV/2).^2).^2;
% v(rV>1/2) = 0;
radius = 1/16; centerx = 7/8;
ind = (xgV-centerx).^2 + ygV.^2 <= radius^2;
v = zeros(size(xgV));
v(ind) = (1-(xgV(ind)-centerx).^2 - ygV(ind).^2).^2;
u = zeros(length(rU),1);
initialv = v;       % store initial value
initialu = u;

%% parameters for plotting a circle
theta = linspace(0,2*pi,200);
xp = cos(theta); yp = sin(theta);

%% Construct interpolation matrices
disp('Constructing interpolation and laplacian matrices');

EbarV = interp2_matrix(x1d, y1d, cpxgbarV, cpygbarV, p, bandV);
Ev = interp2_matrix(x1d, y1d, cpxgV, cpygV, p, bandV);
Eu = interp2_matrix(x1d, y1d, cpxgU, cpygU, p, bandU);
EcpUbandV = interp2_matrix(x1d, y1d, cpxgU, cpygU, p, bandV);               % interpolating value of v from band of v onto positions of cp's of u
EcpVbandU = interp2_matrix(x1d, y1d, cpxgV, cpygV, p,bandU);   %using value of u and positions of cp v

%% Create Laplacian matrices

Lu = laplacian_2d_matrix(x1d,y1d, order, bandU);
Lv = laplacian_2d_matrix(x1d,y1d, order, bandV);
%% 

figure(1); clf;
figure(2); clf;
figure(3); clf;

%% Time-stepping for the heat equation
Tf = 22;
%dt = 0.1*dx^2*min(1/(4*d),1/(2*D));
dt = 0.25/max(D,d)*dx^2;
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps

outer = bdyV;
inner = ~outer;

Iouter = speye(nnz(outer));

distOuter = spdiags(distV(outer), 0, nnz(outer), nnz(outer) );
A = Iouter - EbarV(outer,outer) + 2*Nu/d*(distOuter*Ev(outer,outer));
M = EbarV(outer,inner) - 2*Nu/d*(distOuter*Ev(outer,inner));

for kt = 1:numtimesteps

    t = kt*dt;
        
    % evolve u
    vOnCircle = EcpUbandV*v; % this is where problem occurs error of 6*10^(-4) in vCircle
    unew = u + D*dt*(Lu*u) + dt*(Nu*vOnCircle - Mu*u); % this is where the problem is carried through error = 0.8
    u = Eu*unew; %u should be correct on circle and extended to band2 now.
    
    % evolve v
    v = v + d*dt*(Lv*v) + dt*(v.*(1-v));
    
    % impose v boundary condition
    uOnCircle = EcpVbandU(outer,:)*u;
    b = M*v(inner) + 2*Mu/d*(distOuter*uOnCircle);   
    v(outer) = A\b;
  
  
if mod(kt,100) ==0 || kt < 10 || kt == numtimesteps
    % plot over computation band
    plot2d_compdomain(v(inner), xgV(inner), ygV(inner), dx, dx, 1)
    
    title( ['v: soln at time ' num2str(t) ...
            ', timestep #' num2str(kt)] );
    xlabel('x'); ylabel('y');
    %hold on
    plot(xp,yp,'k-', 'linewidth', 2);
    %axis equal;  axis tight
    %caxis([0 1])
    
    plot2d_compdomain(u, xgU, ygU, dx, dx, 2)
    title( ['u: soln at time ' num2str(t) ...
            ', timestep #' num2str(kt)] );
    
    plot2d_compdomain(vOnCircle, xgU, ygU, dx, dx, 3)
    title( ['vOnCircle: soln at time ' num2str(t) ...
            ', timestep #' num2str(kt)] );
        
    drawnow();  
end
end
