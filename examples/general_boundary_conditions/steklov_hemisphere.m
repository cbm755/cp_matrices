%% Steklov problem on hemisphere 

clear ;
close all ; 


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
p = 1 ;      % Low interp degree (for plotting)
q = 3 ;      % High interp degree (for computation)
order = 2 ;  % Laplacian order
bw = 1.0001*sqrt((dim-1)*((q+1)/2)^2 + ((order/2+(q+1)/2)^2));
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
xp1 = xp(:) ;  yp1 = yp(:) ;  zp1 = zp(:) ;

% % Eplot is a matrix which interpolates data onto the plotting grid
Eplot = interp3_matrix(x1d,y1d,z1d, xp1,yp1,zp1, 3, band) ;


%% Setting up linear system

n = [0;0;-1] ;
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
        w =[xtmp; ytmp; ztmp] - [cpxtmp; cpytmp; cpztmp] ;
        n = [cpxtmp; cpytmp; cpztmp] - [cpbarxtmp; cpbarytmp; cpbarztmp] ;
        n = n/norm(n) ;
        if norm([cpx_band(k)-cpbarx_band(k); cpy_band(k)-cpbary_band(k); cpz_band(k)-cpbarz_band(k)],2) >= 1e-4
            dvec(k) = 2*(w'*n) ;
        end
    end
end

N = length(band) ;
Dmat = spdiags(dvec, 0, N, N) ;

gamma = 2*dim/dx^2 ; % penalty parameter
Imat = speye(size(Ebarqmat)) ;
Amat = Eqmat*Lmat - gamma*(Imat - Ebarqmat) ;
Bmat = -gamma*(Dmat*Eqmat) ;


%% Compute Steklov spectrum 

numrow = 5 ; 
numcol = 6 ; 

[V,D] = eigs(Amat, Bmat, numrow*numcol + 1, 0.1);

D = diag(D);
[Lambda,I] = sort(abs(D));
disp(Lambda)


figure(1) ; clf ;
colormap('cool')

tiledlayout(numrow, numcol)
idxset = 2:numrow*numcol+1 ; 

for i = 1:length(idxset)
    idx = idxset(i) ;
    lambda = Lambda(idx);
    eigenvec = V(:,I(idx));
    
    uplot = Eplot*(Eqmat*real(eigenvec));
    uplot = reshape(uplot, size(xp));
    nexttile
    surf(xp, yp, zp, uplot);    
    xlabel('x') ; ylabel('y') ; zlabel('z');
    title(['ev = ' num2str(lambda)]);
    axis equal ; axis off ;
    shading interp
end

set(gcf, 'Position', [1           1        1512         773])



