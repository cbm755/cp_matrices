function [error_inf_u, error_inf_v] = ball_sphere_diffusion(dx)

tic
if nargin < 1
    dx = 0.2;
end
dx
%dx = 0.2

%%
% v_t = \Delta v + f in \Omega
% \alpha v - \beta u + \frac{\partial v}{\partial n} = 0 on S
% u_t = \Delta_S u + \frac{\partial v}{\partial n} + g on S
R = 1;

alpha = 1;
beta = 1;

k = 1.527338738;
A = sin(k) + cos(k)/k - sin(k)/k^2;
exactfnu = @(t,theta,phi) A*exp(-k^2*t)*cos(phi+pi/2);
% When r goes to 0, exactfnv goes to 0.
exactfnv = @(t,theta,phi,r) exp(-k^2*t)*(sin(k*r) - k*r.*cos(k*r))./(max(k*r,1e-15)).^2.*cos(phi+pi/2);

x1d = (-3:dx:3);
y1d = x1d;
z1d = x1d;

%% Find closest points on the surface or in the volume
% meshgrid is only needed for finding the closest points, not afterwards
[xx, yy, zz] = meshgrid(x1d, y1d, z1d);

% Finding closest points for the unit disk.
[cpxV,cpyV,cpzV, distV,bdyV] = cpSphereInterior(xx,yy,zz,R);
% Finding closest points for the unit circle.
[cpxU, cpyU, cpzU, distU] = cpSphere(xx,yy,zz,R);
% make into vectors
cpxgU = cpxU(:); cpygU = cpyU(:); cpzgU = cpzU(:);
cpxgV = cpxV(:); cpygV = cpyV(:); cpzgV = cpzV(:);

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
cpxgV = cpxgV(bandV); cpygV = cpygV(bandV); cpzgV = cpzgV(bandV);
xgV = xx(bandV); ygV = yy(bandV); zgV = zz(bandV);
distV = distV(bandV); bdyV = bdyV(bandV); 

cpxgU = cpxgU(bandU); cpygU = cpygU(bandU); cpzgU = cpzgU(bandU);
xgU = xx(bandU); ygU = yy(bandU); zgU = zz(bandU);

%% Construct interpolation matrices
disp('Constructing interpolation and laplacian matrices');

EcpVbandU = interp3_matrix(x1d, y1d, z1d, cpxgV(bdyV), cpygV(bdyV), cpzgV(bdyV), p, bandU);   %using value of u and positions of cp v

Lv = laplacian_3d_matrix(x1d, y1d, z1d, order, bandV);
cpxv_bar_bdy = 2*cpxgV(bdyV) - xgV(bdyV);
cpyv_bar_bdy = 2*cpygV(bdyV) - ygV(bdyV);
cpzv_bar_bdy = 2*cpzgV(bdyV) - zgV(bdyV);
EbarV_bdy = interp3_matrix(x1d,y1d,z1d,cpxv_bar_bdy,cpyv_bar_bdy,cpzv_bar_bdy,p,bandV);
Iv = speye(size(Lv));
ng = nnz(bdyV);
D = spdiags(distV(bdyV),0,ng,ng);
Av_bdy =  alpha*D*(Iv(bdyV,:) + EbarV_bdy)/2 + (-EbarV_bdy + Iv(bdyV,:))/2;
Eoo_v = Av_bdy(:,bdyV);
Eoi_v = Av_bdy(:,~bdyV);

Lu = laplacian_3d_matrix(x1d,y1d,z1d,order,bandU);
E1u = interp3_matrix(x1d,y1d,z1d, cpxgU,cpygU,cpzgU, 1, bandU);
E3u = interp3_matrix(x1d,y1d,z1d, cpxgU,cpygU,cpzgU, 3, bandU);
Iu = speye(size(Lu));
Au = (E1u*Lu-6/dx^2*(Iu-E3u)) - beta*Iu;
EcpUbandV = interp3_matrix(x1d, y1d, z1d, cpxgU, cpygU, cpzgU, p, bandV);               % interpolating value of v from band of v onto positions of cp's of u

Buv = beta*D*EcpVbandU;

%% Time stepping
Tf = .1;
dt = 0.125*dx^2;
numtimesteps = ceil(Tf/dt)
dt = Tf/numtimesteps

[xp,yp,zp] = paramSphere(256,R);
[th_u_plot, phi_u_plot] = cart2sph(xp(:),yp(:),zp(:));
Eplot = interp3_matrix( x1d, y1d, z1d, xp(:), yp(:), zp(:), p, bandU );

outer = bdyV;
inner = ~outer;

% set up intial u, v
[th_u,phi_u] = cart2sph(cpxgU,cpygU,cpzgU);
u = exactfnu(0,th_u,phi_u);
[th_v, phi_v, r_v] = cart2sph(xgV,ygV,zgV);
v = exactfnv(0,th_v,phi_v,r_v);

ifplot = 0;
if ifplot
    figure(1);
    figure(2);
    figure(3);
end
for kt = 1:numtimesteps

    t = kt*dt;

    % evolve v
    v = v + dt*(Lv*v);
   
    % impose v boundary condition
    v(outer) = Eoo_v \ (Buv*u - Eoi_v*v(inner));
    
    % evolve u
    vOnCircle = EcpUbandV*v;
    u = u + dt*( Au*u + alpha*vOnCircle );
%     unew = u + dt*( Lu*u + alpha*vOnCircle - beta*u );
%     u = E3u*unew; 

    if mod(kt,100) ==0 || kt < 10 || kt == numtimesteps
        vexact = exactfnv(t,th_v,phi_v,r_v);
        errorV = v - vexact;
        error_inf_v = norm(errorV(inner),inf) / norm(vexact,inf);

        uexact = exactfnu(t,th_u_plot,phi_u_plot);
        uplot = Eplot*u;
        errorU = uplot - uexact;
        error_inf_u = norm(errorU,inf) / norm(uexact,inf);

        [t, error_inf_u, error_inf_v]
         
        if ifplot            
            figure(1)
            surf(xp,yp,zp,uplot)
            title( ['u: soln at time ' num2str(t) ...
                    ', timestep #' num2str(kt)] );
    
            drawnow();
        end
    end
end

[error_inf_u, error_inf_v]
toc

end
