%function [error_inf_u, error_inf_v] = bulk_ellipse_diffusion(dx)

% [err_u, err_v, ratio_u, ratio_v] = conv_study(@bulk_ellipse_diffusion,5,0.4);
dx = 0.05
tic
%if nargin < 1
%    dx = 0.2;
%end

%%
% -\Delta v + v = f in \Omega
% \alpha v - \beta u + \frac{\partial v}{\partial n} = 0 on S
% -\Delta_S u + u + \frac{\partial v}{\partial n} = g on S

alpha = 1;
beta = 1;

a = sqrt(3); b = sqrt(2);
m = 1; n = 2;
exactfnv = @(t,x,y) beta*(exp(-(m^2+n^2)*t).*sin(m*x).*sin(n*y));
dxfnv = @(t,x,y) beta*m*(exp(-(m^2+n^2)*t).*cos(m*x).*sin(n*y));
dyfnv = @(t,x,y) beta*n*(exp(-(m^2+n^2)*t).*sin(m*x).*cos(n*y));
nxfn = @(x,y) x./sqrt(x.^2/a^4 + y.^2/b^4)/a^2;
nyfn = @(x,y) y./sqrt(x.^2/a^4 + y.^2/b^4)/b^2;
                             
exactfnu = @(t,x,y) 1/beta * ( alpha*exactfnv(t,x,y) + dxfnv(t,x,y).*nxfn(x,y) + dyfnv(t,x,y).*nyfn(x,y) );
phi = @(x,y) x.^2/a^2 + y.^2/b^2 - 1;                             
has_t = true;
laplace_beltrami_fun = laplace_beltrami_ls2d(phi,exactfnu,has_t);
% sourcefnu = u_t - \Delta_S u + \beta u - \alpha v
sourcefnu = @(t,x,y) -(m^2+n^2)*exactfnu(t,x,y) - laplace_beltrami_fun(t,x,y) + beta*exactfnu(t,x,y) - alpha*exactfnv(t,x,y);

x1d = (-4:dx:4);
y1d = x1d;

%% Find closest points on the surface or in the volume
% meshgrid is only needed for finding the closest points, not afterwards
[xx, yy] = meshgrid(x1d, y1d);

% Finding closest points for the unit disk.
[cpxV,cpyV,distV,bdyV] = cpEllipseInterior(xx,yy,[a b]);
% Finding closest points for the unit circle.
[cpxU, cpyU, distU] = cpEllipse(xx,yy,a,b);
% make into vectors
cpxgU = cpxU(:); cpygU = cpyU(:);
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
cpxgV = cpxgV(bandV); cpygV = cpygV(bandV);
xgV = xx(bandV); ygV = yy(bandV);
distV = distV(bandV); bdyV = bdyV(bandV);

cpxgU = cpxgU(bandU); cpygU = cpygU(bandU);
xgU = xx(bandU); ygU = yy(bandU);

%% Construct interpolation matrices
disp('Constructing interpolation and laplacian matrices');

EcpVbandU = interp2_matrix(x1d, y1d, cpxgV(bdyV), cpygV(bdyV), p,bandU);   %using value of u and positions of cp v

Lv = laplacian_2d_matrix(x1d,y1d, order, bandV);
cpxv_bar_bdy = 2*cpxgV(bdyV) - xgV(bdyV);
cpyv_bar_bdy = 2*cpygV(bdyV) - ygV(bdyV);
EbarV_bdy = interp2_matrix(x1d,y1d,cpxv_bar_bdy,cpyv_bar_bdy,p,bandV);
Iv = speye(size(Lv));
ng = nnz(bdyV);
D = spdiags(distV(bdyV),0,ng,ng);
Av_bdy =  alpha*D*(Iv(bdyV,:) + EbarV_bdy)/2 + (-EbarV_bdy + Iv(bdyV,:))/2;
Eoo_v = Av_bdy(:,bdyV);
Eoi_v = Av_bdy(:,~bdyV);

Lu = laplacian_2d_matrix(x1d,y1d, order, bandU);
E1u = interp2_matrix(x1d,y1d,cpxgU, cpygU, 1, bandU);
E3u = interp2_matrix(x1d,y1d,cpxgU, cpygU, 3, bandU);
Iu = speye(size(Lu));
Au = (E1u*Lu-4/dx^2*(Iu-E3u)) - beta*Iu;
EcpUbandV = interp2_matrix(x1d, y1d, cpxgU, cpygU, p, bandV);               % interpolating value of v from band of v onto positions of cp's of u

Buv = beta*D*EcpVbandU;

%% Time stepping
Tf = 0.5;
dt = 0.25*dx^2;
numtimesteps = ceil(Tf/dt)
dt = Tf/numtimesteps
[xp,yp] = paramEllipse(256,a,b);
xp = xp(:); yp = yp(:);
Eplot = interp2_matrix( x1d, y1d, xp, yp, p, bandU );

outer = bdyV;
inner = ~outer;

% set up intial u, v
u = exactfnu(0,cpxgU,cpygU);
v = exactfnv(0,xgV,ygV);

t = 0;
ifplot = 0;

if 1==0
h1 = figure(1);
h2 = figure(2);

plot2d_compdomain(v(inner), xgV(inner), ygV(inner), dx, dx, 1)
hold on
plot(xp',yp','-k','LineWidth',2);
axis equal;  axis tight
xlabel('\fontsize{15}x'); ylabel('\fontsize{15}y');
set(gca,'Fontsize',12)
print(h1,'-dpng','-r450','bulk-ellipse-v-0.png');  


plot2d_compdomain(u, xgU, ygU, dx, dx, 2)
hold on
plot(xp',yp','-k','LineWidth',2);
axis equal;  axis tight
xlabel('\fontsize{15}x'); ylabel('\fontsize{15}y');
set(gca,'Fontsize',12)
print(h2,'-dpng','-r450','bulk-ellipse-u-0.png');  
end

for kt = 1:numtimesteps

    t = kt*dt;
    % evolve v
    v = v + dt*(Lv*v);
   
    % impose v boundary condition
    v(outer) = Eoo_v \ (Buv*u - Eoi_v*v(inner));
    
    % evolve u
    vOnCircle = EcpUbandV*v; 
    u = u + dt*( Au*u + alpha*vOnCircle + sourcefnu(t,cpxgU,cpygU) );
%     unew = u + dt*( Lu*u + alpha*vOnCircle - beta*u + sourcefnu(t,cpxgU,cpygU) ); 
%     u = E3u*unew; 
    
    if mod(kt,100) ==0 || kt < 10 || kt == numtimesteps
        vexact = exactfnv(t,xgV(inner),ygV(inner));
        errorV = v(inner) - vexact;
        error_inf_v = norm(errorV,inf) / norm(vexact,inf);

        %uexact = exactfnu(thetas);
        uexact = exactfnu(t,xp,yp);
        errorU = Eplot*u - uexact;
        error_inf_u = norm(errorU,inf) / norm(uexact,inf);

        %[t, error_inf_u, error_inf_v]
         
        if ifplot
            % plot over computation band
            plot2d_compdomain(v(inner), xgV(inner), ygV(inner), dx, dx, 1)
            %plot2d_compdomain(v, xgV, ygV, dx, dx, 1)

            title( ['v: soln at time ' num2str(t) ...
                    ', timestep #' num2str(kt)] );
            xlabel('x'); ylabel('y');
            hold on
            plot(xp',yp','-k','LineWidth',2);
            axis equal;  axis tight
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
end

if 1==0
h3 = figure(3);
h4 = figure(4);

plot2d_compdomain(v(inner), xgV(inner), ygV(inner), dx, dx, 3)
hold on
plot(xp',yp','-k','LineWidth',2);
axis equal;  axis tight
xlabel('\fontsize{15}x'); ylabel('\fontsize{15}y');
set(gca,'Fontsize',12)
print(h3,'-dpng','-r450','bulk-ellipse-v-Tf.png');  


plot2d_compdomain(u, xgU, ygU, dx, dx, 4)
hold on
plot(xp',yp','-k','LineWidth',2);
axis equal;  axis tight
xlabel('\fontsize{15}x'); ylabel('\fontsize{15}y');
set(gca,'Fontsize',12)
print(h4,'-dpng','-r450','bulk-ellipse-u-Tf.png');  
end


[error_inf_u, error_inf_v]
toc

%end