%% Gray--Scott on mobius strip
% Forward Euler time-stepping

clear ; close all ; % clc ;

rng(1989)

dx = 0.05;
% run make_mobius_grid to make this file
load(['mobius_dx=', num2str(dx), '.mat']) ;


%% Banding: do calculation in a narrow band around the surface

dim = 3 ;    % dimension
p = 1 ;      % Low interp degree
q = 3 ;      % High interp degree
order = 2 ;  % Laplacian order
bw = 1.0001*sqrt((dim-1)*((q+1)/2)^2 + ((order/2+(q+1)/2)^2));
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

% Interpolation matrice for modified closest point (Low and high order)
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


%% Parameters and nonlinearities for Gray--Scott

% F = 0.045 ;  k = 0.06 ;  Du = 8e-5 ;  Dv = 4e-5 ; % spot (NBC)
% F = 0.052 ;  k = 0.063 ;  Du = 8e-5 ;  Dv = 4e-5 ; % stripe (NBC)
% F = 0.037 ;  k = 0.06 ; Du = 8e-5 ; Dv = 4e-5 ; % Labyrinthine (NBC)
% F = 0.020 ;  k = 0.052 ;  Du = 8e-5 ;  Dv = 4e-5 ; % chaos
F = 0.010 ;  k = 0.042 ;  Du = 8e-5 ;  Dv = Du/2.5 ; % spiral wave with kappa = 70

% Robin coefficient
kappa = 10 ; 

% Reaction terms
f = @(u,v) (-u.*v.*v + F*(1-u)) ;
g = @(u,v) ( u.*v.*v - (F+k)*v) ;


%% Initial conditions - small perturbation from steady state

um = real(1 - sqrt(1 - 4*(F+k)^2/F))/2 ;
vm = (F + k)/um ;
uvec = um + 0.1*(rand(size(xx_band)) - 0.5) ;
vvec = vm + 0.1*(rand(size(xx_band)) - 0.5) ;


%% Matrices for time-stepping

dvec = zeros(size(band)) ; 
for k = 1:length(band)
    if bdy_band(k) == 1  
        % Compute unit outward conormal
        % u = atan2(cpy_band(k), cpx_band(k) + R*T/2) ;
        % n = [R*T*cos(u/2)*cos(u) ; 
        %      R*T*cos(u/2)*sin(u) ; 
        %      R*T*sin(u/2)] ;
        % n = n/norm(n) ;

        % % Compute unit outward conormal
        n = [cpx_band(k); cpy_band(k); cpz_band(k)] - [cpbarx_band(k); cpbary_band(k); cpbarz_band(k)]  ;
        n = n/norm(n) ;

        % Compute vector pointing from surface to boundary grid
        w = [xx_band(k); yy_band(k); zz_band(k)] - [cpx_band(k); cpy_band(k); cpz_band(k)] ;
        % Make boundary normal outward
        if dot(n, w) < 0 
            n = -n ; 
        end 

	% See Remark 2 of [WongMacdonaldLee2026]
        if norm([cpx_band(k)-cpbarx_band(k); cpy_band(k)-cpbary_band(k); cpz_band(k)-cpbarz_band(k)],2) >= 1e-4
            dvec(k) = 2*(w'*n) ;
        end
    end
end


numpts = length(band) ;
Dmat = spdiags(dvec, 0, numpts, numpts) ;

gamma = 2*dim*max(Du,Dv)/(dx^2);
Imat = speye(size(Ebarqmat)) ;
Au = Du*(Eqmat*Lmat) - gamma*(Imat - Ebarqmat - Dmat*(-kappa*Eqmat)) ;
Av = Dv*(Eqmat*Lmat) - gamma*(Imat - Ebarqmat - Dmat*(-kappa*Eqmat)) ;


%% Record video

videoname = ['GS_mobius_kappa=', num2str(kappa), '_dx=', num2str(dx)] ;
writerObj = VideoWriter(videoname, 'MPEG-4') ;
writerObj.FrameRate = 30 ; 

% Open the video writer
open(writerObj);


%% Time-stepping

Tf = 10000 ; % Simulation time
dt = .1 * (1/max(Du,Dv))*dx^2 ; % Time step
numtimesteps = ceil(Tf/dt) ; % Number of time steps
dt = Tf/numtimesteps ; % Adjust for integer number of steps


figure(1) ; clf ;
sphplot = Eplot*uvec;
sphplot = reshape(sphplot, size(xp));
Hplot = surf(xp, yp, zp, sphplot);
axis equal ; axis off ;
view(-38.9, 33.6)
shading interp
colorbar
colormap('jet')
clim([0,1])
writeVideo(writerObj, getframe(gcf)) ; 


t = 0 ;
kt  = 0 ;

while t < Tf

  uvec0 = uvec ;
  vvec0 = vvec ;

  uvec = uvec0 + dt*( Eqmat*f(uvec0,vvec0) + Au*uvec0 ) ;
  vvec = vvec0 + dt*( Eqmat*g(uvec0,vvec0) + Av*vvec0 ) ;
  

  if ( (mod(kt,25)==0) || (kt<=10) || (kt==numtimesteps) )
    disp([kt t]);
    sphplot = Eplot*uvec;
    sphplot = reshape(sphplot, size(xp));
    set(0, 'CurrentFigure', 1);
    set(Hplot, 'CData', sphplot);
    title( ['u at time ' num2str(t) ', kt= ' num2str(kt)] );
    drawnow() ;
    writeVideo(writerObj, getframe(gcf)) ;     
  end

  t = t + dt;
  kt = kt + 1 ; 

end


close(writerObj)
