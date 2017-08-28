%% Project Vase with constant SA linechart

clear, close all;

%% Parameters
tol = 10^-2;

% cpf = @cpVase;
% paramf = @paramVase;


a=1.5; % length along z
b=1; % length along xy plane
ab=[a,b];

cpf = @(x,y,z) cpEllipsoid(x, y, z, [a b], [], 'z');
paramf = @(N) paramEllipsoid(N, [a b], [], 'z');


max_vec=[];
sum_vec=[];
dx=0.1; % grid size
%SA_scale=0.05;

%[theta ,beta]=meshgrid(linspace(0,2*pi,5),linspace(0,pi,5));

%theta=[0];
theta= linspace(-pi/2,pi/2,1);
beta=0;

% xh=linspace(0,1,3);
% yh=linspace(0,1,3);
% zh=linspace(0,1,3);

for j= 1:length(theta)

% theta=-pi/2;%south pole
% beta= 0;%changing polar 

xh=b*cos(theta(j)).*cos(beta);
yh=b*cos(theta(j)).*sin(beta);
zh=a*sin(theta(j));

% SA_wanted=[0;0; 0.1; 0];
% holes = [0 0 a;
%         0 0 -a;
%          xh yh zh;
%          0 b 0];
SA_wanted=0.1;
holes = [xh yh zh];


%% Construct a grid in the embedding space

%dx = 0.025;                   % grid size

% make vectors of x, y, positions of the grid
x1d = (-2.0:dx:2.0)';
y1d = x1d;
z1d = x1d;

nx = length(x1d);
ny = length(y1d);
nz = length(z1d);


%% Banding: do calculation in a narrow band around the vase
dim = 3;    % dimension
p = 3;      % interpolation degree
order = 2;  % Laplacian order
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

%% Find closest points on the surface
% For each point (x,y,z), we store the closest point on the sphere
% (cpx,cpy,cpz)
% meshgrid is only needed for finding the closest points, not afterwards
[xx,yy,zz] = meshgrid(x1d, y1d, z1d);
xx = xx(:);  yy = yy(:);  zz = zz(:);

[cpx,cpy,cpz,dist] = cpf(xx, yy, zz);

band = find(abs(dist) <= bw*dx);
cpx = cpx(band); cpy = cpy(band); cpz = cpz(band);
x = xx(band); y = yy(band); z = zz(band);


tic
[cpx,cpy,cpz,dist,bdy] = cpCutHole(x, y, z, cpf, holes, SA_wanted, dx, bw);
toc
%[cpx,cpy,cpz,dist] = cpFromTriSlow(hole_cen(1), hole_cen(2), hole_cen(3), Faces, Vertices);
% make into vectors


%% Extension of closest point
tic
E = interp3_matrix(x1d, y1d, z1d, cpx, cpy, cpz, p, band);

% bdy cond;
E(bdy,:) = -E(bdy,:);
 
%% solve PDE using xg, yg, zg
L = laplacian_3d_matrix(x1d,y1d,z1d, order, band, band);
toc

%% For constant D
d = 1;
D = -1/d*ones(length(E),1);

 
%% solution of u
% boundary condition

gamma = 2*3/((dx)^2);
I=speye(size(E));
M=E*L-gamma*(I-E);
tic
u = M \ D;
toc
%u = gmres(M,D);


%% plot

N=100;
[xp, yp, zp] = paramf(N);
Eplot = interp3_matrix(x1d, y1d, z1d, xp(:), yp(:), zp(:), p, band);

up = Eplot*u;
up=reshape(up,size(xp));

%% Plot ellipsoid with trap
surf(xp,yp,zp,max(up,0))
colormap('jet')
shading interp
axis equal;
set(gcf,'color','w');

%% Max of u
%  max_up=max(max(up));
%  max_vec=[max_vec max_up]
 
%% Expected MFPT
% sum_u=sum(u)/sum(bdy(:));
%  sum_vec=[sum_vec sum_u]
end


% h=plot(theta(:),max_vec,'o-')
%  set(h,'linewidth',4)
%  set(gca,'fontsize',20)

%% Max of u
%  plot(theta(:),max_vec,'o-')
% xlabel('\theta')
% ylabel('Maximum MFPT')
% axis equal;
% set(gcf,'color','w');

%% Expected MFPT
%  plot(theta(:),sum_vec,'o-')
%  xlabel('\theta')
%  ylabel('Expected MFPT')
%  %axis equal;
%  set(gcf,'color','w');











