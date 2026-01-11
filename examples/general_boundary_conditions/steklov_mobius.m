%% Steklov problem on mobius strip

clear ; close all ; % clc ; 


%% Construct a grid in the embedding space

% % grid size
% dx = 0.05 ;
% 
% % make vectors of x, y, positions of the grid
% % x1d = (-1.8:dx:1.8)' ;
% % y1d = x1d ;
% % z1d = x1d ;
% 
% 
% [xx, yy, zz] = meshgrid(x1d, y1d, z1d) ;
% 
% % Find closest points on the surface
% R = 1 ; % radius of center circle
% T = 0.35 ; % thickness
% cpf = @cpMobiusStrip ; 
% [cpbarx, cpbary, cpbarz, dist, bdy] = cpbar_3d(xx, yy, zz, cpf, R, T) ;
% [cpx, cpy, cpz, ~, ~] = cpf(xx, yy, zz) ;

dx = 0.025 ;

load(['mobius_dx=', num2str(dx), '.mat']) ;


%% Banding: do calculation in a narrow band around the surface

dim = 3 ;    % dimension
p = 1 ;      % interpolation degree of differential operator
q = 3 ;      % interpolation degree of penalty term
order = 2 ;  % Laplacian order
bw = 1.0001*sqrt((dim-1)*((max(p,q)+1)/2)^2 + ((order/2+(max(p,q)+1)/2)^2));
band = find(abs(dist) < bw*dx);


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

Lmat = laplacian_3d_matrix(x1d, y1d, z1d, order, band) ;


%% Construct an interpolation matrix for plotting on sphere

% % plotting grid on sphere, based on parametrization
[xp,yp,zp] = paramMobiusStrip(200) ;
xp1 = xp(:) ; yp1 = yp(:) ; zp1 = zp(:) ;

% % Eplot is a matrix which interpolates data onto the plotting grid
Eplot = interp3_matrix(x1d, y1d, z1d, xp1, yp1, zp1, 3, band) ;


%% Setting up the linear system

dvec = zeros(size(band)) ; 
for k = 1:length(band)
    if bdy_band(k) == 1  


        % % Compute unit outward conormal
        % u = atan2(cpy_band(k), cpx_band(k) + R*T/2) ;
        % n = [R*T*cos(u/2)*cos(u) ; 
        %      R*T*cos(u/2)*sin(u) ; 
        %      R*T*sin(u/2)] ;
        % n = n/norm(n) ;
        
        % Approximate unit outward conormal
        n = [cpx_band(k); cpy_band(k); cpz_band(k)] - [cpbarx_band(k); cpbary_band(k); cpbarz_band(k)]  ;
        n = n/norm(n,2) ;


        % Compute vector pointing from surface to boundary grid
        w = [xx_band(k); yy_band(k); zz_band(k)] - [cpx_band(k); cpy_band(k); cpz_band(k)] ;
        % Make boundary normal outward
        if dot(n, w) < 0 
            n = -n ; 
        end

        if norm([cpx_band(k)-cpbarx_band(k); cpy_band(k)-cpbary_band(k); cpz_band(k)-cpbarz_band(k)],2) >= 1e-4
            dvec(k) = 2*(w'*n) ;
        end

    end
end

% axis equal ; axis tight ; axis off

numpts = length(band) ;
Dmat = spdiags(dvec, 0, numpts, numpts) ;

gamma = 2*dim/dx^2 ; % penalty parameter
Imat = speye(size(Ebarqmat)) ;
Amat = Ebarqmat*Lmat - gamma*(Imat - Ebarqmat) ;
Bmat = -(gamma)*(Dmat*Eqmat) ;


%% Compute Steklov spectrum 

numrow = 5 ; 
numcol = 6 ; 

[V,D] = eigs(Amat, Bmat, numrow*numcol + 1, 0.1);

D = diag(D);
[Lambda,I] = sort(abs(D));
disp(Lambda)


figure(2) ; clf ;
colormap('cool')

tiledlayout(numrow, numcol)
idxset = 2:numrow*numcol+1 ; 

for i = 1:length(idxset)
    idx = idxset(i) ;
    lambda = Lambda(idx);
    eigenvec = V(:,I(idx));
    
    uplot = Eplot*Eqmat*real(eigenvec);
    uplot = reshape(uplot, size(xp));
    nexttile
    surf(xp, yp, zp, uplot);    
    xlabel('x') ; ylabel('y') ; zlabel('z');
    title(['ev = ' num2str(lambda)]);
    axis equal ; axis off ;
    shading interp
end

set(gcf, 'Position', [1           1        1512         773])

