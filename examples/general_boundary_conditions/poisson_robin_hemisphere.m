%% Poisson equations on hemisphere \Delta_S u = f
% Robin BC: \partial_n u + \kappa u = g

% Azimuth angle: phi \in [0,2\pi)
% Polar angle: theta \in [0,\pi]

% Solution: u = cos(2*phi)*sin(theta)^2 + sin(3*phi)*sin(theta)^3
% RHS of Poisson eqn: f = -6*cos(2*phi)*sin(theta)^2 - 12*sin(3*phi)*sin(theta)^3
% Robin coeff: kappa = 1
% RHS of Robin BC: g = cos(2*phi) + sin(3*phi)


clear ; close all ; % clc ;


%% Construct a grid in the embedding space

dx = 0.05 ;  % grid spacing
R = 1 ; %Radius of the circle

pad = 5 ;
x1d = ((-R-pad*dx):dx:(R+pad*dx))' ;
y1d = x1d ;
z1d = y1d ; 
nx = length(x1d) ;
ny = length(y1d) ;
nz = length(z1d) ;
[xx, yy, zz] = meshgrid(x1d,y1d,z1d) ;


%% Find closest points on the surface

cpf = @cpHemisphere ;
[cpbarx, cpbary, cpbarz, dist, bdy] = cpbar_3d(xx, yy, zz, cpf, R) ;
[cpx, cpy, cpz, ~, ~] = cpHemisphere(xx, yy, zz) ;


%% Banding

dim = 3;
p = 1 ; % Low interp degree
q = 3 ; % High interp degree
order = 2;  % Laplacian order
bw = 1.0001*sqrt((dim-1)*((q+1)/2)^2 + ((order/2+(q+1)/2)^2)) ; % Bandwidth = bw * dx
if (pad < bw)
  warning('maybe insufficient padding');
end
band = find(abs(dist) <= bw*dx);
outband = find(abs(dist) > bw*dx);

% Last place we need the 3D array
xx_band = xx(band) ;
yy_band = yy(band) ;
zz_band = zz(band) ;
cpbarx_band = cpbarx(band) ;
cpbary_band = cpbary(band) ;
cpbarz_band = cpbarz(band) ; 
dist_band = dist(band) ;
bdy_band = bdy(band) ;

cpx_band = cpx(band) ;
cpy_band = cpy(band) ;
cpz_band = cpz(band) ; 


%% Construct an interpolation matrix for closest point 

% Interpolation matrice for modified closest point 
Ebarqmat = interp3_matrix(x1d,y1d,z1d, cpbarx_band, cpbary_band, cpbarz_band, q, band);

% Interpolation matrice for closest point 
Eqmat = interp3_matrix(x1d,y1d,z1d, cpx_band, cpy_band, cpz_band, q, band);


%% Create Laplacian matrix 

Lmat = laplacian_3d_matrix(x1d,y1d,z1d, order, band, band);


%% plotting grid

[xp,yp,zp] = paramHemisphere(100, R) ;
xp1 = xp(:);  yp1 = yp(:);  zp1 = zp(:);
Eplot = interp3_matrix(x1d,y1d,z1d, xp1,yp1,zp1, p, band);


%% Convert Cartesian coordinates of band points to spherical coordinate


[az_band, el_band] = cart2sph(cpbarx_band, cpbary_band, cpbarz_band) ; 

% Convert to our convention 
az_band = az_band + pi ; % So that [-pi,pi] -> [0,2*pi]
el_band = pi/2 - el_band ; % Colatitude -> Polar



%% Setting Robin problem \Delta u = f with \partial_n u + kappa u = g at bdy


u_exact = sin(el_band).^2.*cos(2*az_band) + sin(el_band).^3.*sin(3*az_band) ; % Exact solution
fvec = -6*sin(el_band).^2.*cos(2*az_band) - 12*sin(el_band).^3.*sin(3*az_band) ; % RHS of Poisson eqn
kappa = 1 ; % Robin coefficient
gvec = cos(2*az_band) + sin(3*az_band) ; % RHS of inhomogeneous Robin BC


%% Setting up the linear system

% n =  [0;0;-1] ; % Unit outward conormal

dvec = zeros(size(band)) ; 
for k = 1:length(band)
    if bdy_band(k) == 1      
        xtmp = xx_band(k) ; 
        ytmp = yy_band(k) ; 
        ztmp = zz_band(k) ;
        cpxtmp = cpx_band(k) ; 
        cpytmp = cpy_band(k) ; 
        cpztmp = cpz_band(k) ;         
        cpbarxtmp = cpbarx_band(k) ; 
        cpbarytmp = cpbary_band(k) ; 
        cpbarztmp = cpbarz_band(k) ;
        w = [xtmp; ytmp; ztmp] - [cpxtmp; cpytmp; cpztmp] ;
        n = [cpxtmp; cpytmp; cpztmp] - [cpbarxtmp; cpbarytmp; cpbarztmp] ; 
        n = n/norm(n) ;
        if norm([cpx_band(k)-cpbarx_band(k); cpy_band(k)-cpbary_band(k); cpz_band(k)-cpbarz_band(k)],2) >= 1e-4
            dvec(k) = 2*(w'*n) ;
        end
    end
end

numpts = length(band) ;
Dmat = spdiags(dvec, 0, numpts, numpts) ;

gamma = 2*dim/dx^2 ; % penalty parameter
Imat = speye(size(Ebarqmat)) ;
Amat = Eqmat*Lmat - gamma*(Imat - Ebarqmat - Dmat*(-kappa*Eqmat)) ;

u = Amat \ (fvec - gamma*(dvec.*gvec)) ;


%% Plot numerical and exact solution and their diffference

figure(1) ; clf ; 
subplot(1,2,1)
sphplot = Eplot*u ;
sphplot = reshape(sphplot, size(xp)) ;
surf(xp, yp, zp, sphplot) ; hold on ;
axis equal ; axis tight ; 
xlabel('x') ; ylabel('y') ; zlabel('z')
title('Numerical solution')
shading interp
colorbar('North')


subplot(1,2,2)
sphplot = Eplot*u_exact ;
sphplot = reshape(sphplot, size(xp)) ;
surf(xp, yp, zp, sphplot) ; hold on ;
axis equal ; axis tight ; 
xlabel('x') ; ylabel('y') ; zlabel('z')
title('Exact solution')
shading interp
colorbar('North')



%% Error analysis

diff = abs(u - u_exact) ;
upper = find(zz_band > 0) ;
diff = diff(upper) ; 
[error, idx] = max(diff) ; 
error = error / max(abs(u_exact)) ;
disp([dx, error])

figure(2) ; clf ; 
sphplot = Eplot*abs(u - u_exact) ;
sphplot = reshape(sphplot, size(xp)) ;
surf(xp, yp, zp, sphplot) ; hold on ;
axis equal ; axis tight ; 
xlabel('x') ; ylabel('y') ; zlabel('z')
title('error')
shading interp
colorbar('North')