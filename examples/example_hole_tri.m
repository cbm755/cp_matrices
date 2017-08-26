%% Project bunny 3d

%clear, close all;

%% Parameters
tol = 10^-5;

%cen=[0,0,0]; % center of vase
%a=1.5; % length along z
%b=1; % length along xy plane
%ab=[a,b];
%lim=[-pi/2, pi/2]; 
n=8;
%[xh,yh,zh] = paramVase(n, lim, ab, cen);
max_vec=[];
bdy_vec=[];




% bunny
%hole_cen = [-0.2518716,0.2134296,0.07213944];
% pig
%hole_cen = [0.5, -3, 0.5];
hole_cen = [0.1, 3, 0.5];%bunny ear
%hole_cen = [0.5, 0.5, 0.5];%pig ear
%hole_cen = [-0.5, -0.7, 0.3];%bunny back
%hole_cen = [0.3, 0.3,-0.5];%pig hip


%% Construct a grid in the embedding space
if (1==1)
dx=0.2;
% ply file contains the triangles
disp('reading plyread');
%PlyFile = 'bunny.ply';
%PlyFile = 'pig_loop2.ply';
PlyFile = 'annies_pig.ply';
[Faces, Vertices] = plyread(PlyFile, 'tri');

disp('running tri2cp');
[IJK,DIST,CP,XYZ] = tri2cp(Faces, Vertices, dx, -2);
i = IJK(:,1);
j = IJK(:,2);
k = IJK(:,3);
dist = DIST;
cpx = CP(:,1);
cpy = CP(:,2);
cpz = CP(:,3);
x = XYZ(:,1);
y = XYZ(:,2);
z = XYZ(:,3);

x1d=-2.0:dx:2.0;
y1d=x1d;
z1d=x1d;
nx=length(x1d);
ny=length(y1d);
nz=length(z1d);
end



cpf = @(x,y,z) cpFromTriSlow(x, y, z, Faces, Vertices);


%% Debugging plot
V = Vertices;
figure(1); clf; trisurf(Faces, V(:,1), V(:,2), V(:,3));
shading flat
[xp,yp,zp] = sphere(20);
xp = rh*xp + cx;
yp = rh*yp + cy;
zp = rh*zp + cz;
hold on;
surf(xp, yp, zp, 10+zp);


figure(2)
%porcupine_plot3d(x,y,z, cpx,cpy,cpz, bdy, [], 2)




%% Banding: do calculation in a narrow band around the vase
dim = 3;    % dimension
p = 3;      % interpolation degree
order = 2;  % Laplacian order
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
% band = find(abs(dist) <= bw*dx);

% TODO: bw is needed by cpCutHole, but here it is determined implicitly by the tri2cp code...


[cpx, cpy, cpz, dist, bdy] = cpCutHole(x,y,z,cpf, hole_cen, 1, 1e-10, dx, bw);

M = bdy;

figure(4); clf;
plot3(cpx(M), cpy(M), cpz(M), 'k.', 'markersize', 1);
hold on;
trisurf(Faces, V(:,1), V(:,2), V(:,3));
shading flat


band = sub2ind([ny,nx,nz], j,i,k);
% Eplot = interp3_matrix(x1d, y1d,z1d, xp(:), yp(:),zp(:), p, band);
 
%% Extension of closest point
E = interp3_matrix(x1d, y1d, z1d, cpx, cpy, cpz, p, band);

% bdy cond;
% using the new "cp bar" definition; 
E(bdy,:) = -E(bdy,:);
 
%% solve PDE using xg, yg, zg
L = laplacian_3d_matrix(x1d,y1d,z1d, order, band, band);

%% For constant D
d = 1;
D = -1/d*ones(length(E),1);

%% solution of u
% boundary condition
gamma = 2*3/((dx)^2);
I=speye(size(E));
M=E*L-gamma*(I-E);
u = M \ D;


 
% surf(xg,yg,zg,u);
%% plot


xp = Vertices(:,1);
yp = Vertices(:,2);
zp = Vertices(:,3);
% N=100;
% [xp,yp,zp] = paramVase(N, lim, ab, cen);
Eplot = interp3_matrix(x1d, y1d,z1d , xp(:), yp(:),zp(:), p, band);

up = Eplot*u;
up=reshape(up,size(xp));

figure(5); clf;
trisurf(Faces,xp,yp,zp,max(up,0))
material dull
camlight left
axis equal;
set(gcf,'color','w');
%surf(xp,yp,zp,max(up,0))
colormap('jet')
material dull
shading interp
axis equal;
hold on;
plot3(cpx(bdy), cpy(bdy), cpz(bdy), 'k.', 'markersize', 1)

 max_up=max(max(up));
 max_vec=[max_vec max_up]
 
 %% Max of u
 plot(theta(:),max_vec,'o-')
xlabel('\theta')
ylabel('Maximum MFPT')
axis equal;
set(gcf,'color','w');

