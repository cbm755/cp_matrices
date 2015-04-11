function [error_inf_u, error_inf_v] = bulk_surface_poisson(dx)

if nargin < 1
    dx = 0.1;
end
%%
% -\Delta v + v = f in \Omega
% \alpha v - \beta u + \frac{\partial v}{\partial n} = 0 on S
% -\Delta_S u + u + \frac{\partial v}{\partial n} = g on S
R = sqrt(2);

alpha = 1;
beta = 1;

exactfnv = @(x,y) beta*exp(-x.*(x-1).*y.*(y-1));
dxfnv = @(x,y) exactfnv(x,y).*(1-2*x).*y.*(y-1);
dyfnv = @(x,y) exactfnv(x,y).*(1-2*y).*x.*(x-1);
nxfn = @(x,y) x./sqrt(x.^2+y.^2);
nyfn = @(x,y) y./sqrt(x.^2+y.^2);
rhsfnv = @(x,y) - exactfnv(x,y).*( (1-2*x).^2.*y.^2.*(y-1).^2 - 2*y.*(y-1) + ...
                                 (1-2*y).^2.*x.^2.*(x-1).^2 - 2*x.*(x-1) ) + exactfnv(x,y);
                             
exactfnu = @(x,y) 1/beta * ( alpha*exactfnv(x,y) + dxfnv(x,y).*nxfn(x,y) + dyfnv(x,y).*nyfn(x,y) );
phi = @(x,y) x.^2 + y.^2 - R^2;                             
laplace_beltrami_fun = laplace_beltrami_ls2d(phi,exactfnu);
rhsfnu = @(x,y) -laplace_beltrami_fun(x,y) + exactfnu(x,y) - alpha*exactfnv(x,y) + beta*exactfnu(x,y);


% k = 2; q = 5;
% exactfnv = @(theta,r) r.^q.*sin(k*theta);
% rhsfnv = @(theta,r) - (q^2-k^2)*r.^(q-2).*sin(k*theta) + exactfnv(theta,r);
% exactfnu = @(theta) 1/beta * ( alpha*R^q + q*R^(q-1) )*sin(k*theta) ;
% rhsfnu = @(theta) k^2/beta*(alpha*R^q+q*R^(q-1))/R^2*sin(k*theta) + (1+beta)*exactfnu(theta) - alpha*exactfnv(theta,R);


x1d = (-4:dx:4);
y1d = x1d;

%% Find closest points on the surface or in the volume
% meshgrid is only needed for finding the closest points, not afterwards
[xx, yy] = meshgrid(x1d, y1d);

% Finding closest points for the unit disk.
[cpxV,cpyV,distV,bdyV] = cpCircleInterior(xx,yy,R);
% Finding closest points for the unit circle.
[cpxU, cpyU, distU] = cpCircle(xx,yy,R);
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
Av = -Lv + Iv;
Av(bdyV,:) =  alpha*D*(Iv(bdyV,:) + EbarV_bdy)/2 + (-EbarV_bdy + Iv(bdyV,:))/2;


Lu = laplacian_2d_matrix(x1d,y1d, order, bandU);
E1u = interp2_matrix(x1d,y1d,cpxgU, cpygU, 1, bandU);
E3u = interp2_matrix(x1d,y1d,cpxgU, cpygU, 3, bandU);
Iu = speye(size(Lu));
Au = -(E1u*Lu-4/dx^2*(Iu-E3u)) + (1+beta)*Iu;
EcpUbandV = interp2_matrix(x1d, y1d, cpxgU, cpygU, p, bandV);               % interpolating value of v from band of v onto positions of cp's of u

Buv = sparse(length(bandV),length(bandU)); 
Buv(bdyV,:) = -beta*D*EcpVbandU;
A = [Av, Buv; -alpha*EcpUbandV, Au];

%[thV,rV] = cart2pol(xgV,ygV);
%f = rhsfnv(thV,rV);
f = rhsfnv(xgV,ygV);
f(bdyV) = 0;

%thU = cart2pol(cpxgU,cpygU);
%g = rhsfnu(thU);
g = rhsfnu(cpxgU,cpygU);

rhs = [f; g];

soln = A \ rhs;

v = soln(1:length(bandV));
u = soln(length(bandV)+1:end);

%[th,r] = cart2pol(xgV(~bdyV),ygV(~bdyV));
%vexact = exactfnv(th,r);
vexact = exactfnv(xgV(~bdyV),ygV(~bdyV));
errorV = v(~bdyV) - vexact;
error_inf_v = norm(errorV,inf) / norm(vexact,inf);

thetas = linspace(0, 2*pi, 100)';
r = R*ones( size(thetas) );
[xp, yp] = pol2cart(thetas, r);
xp = xp(:); yp = yp(:);
Eplot = interp2_matrix( x1d, y1d, xp, yp, p, bandU );

%uexact = exactfnu(thetas);
uexact = exactfnu(xp,yp);
errorU = Eplot*u - uexact;
error_inf_u = norm(errorU,inf) / norm(uexact,inf);

[error_inf_u, error_inf_v]

end











