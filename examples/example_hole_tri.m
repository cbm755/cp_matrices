%% Example of cutting a hole in a triangulated shape


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
dx=0.1/2;
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
figure(1); clf;
title('initial hole location')
trisurf(Faces, V(:,1), V(:,2), V(:,3));
shading flat
[xp,yp,zp] = sphere(20);
xp = 0.1*xp + hole_cen(1);
yp = 0.1*yp + hole_cen(2);
zp = 0.1*zp + hole_cen(3);
hold on;
surf(xp, yp, zp, 3+zp);
axis equal




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


tic
[cpx, cpy, cpz, dist, bdy] = cpCutHole(x, y, z, cpf, hole_cen, 0.2, dx, bw, {cpx, cpy,cpz});
toc

figure(4); clf;
title('outline of hole')
plot3(cpx(bdy), cpy(bdy), cpz(bdy), 'k.', 'linewidth', 3);
hold on;
trisurf(Faces, V(:,1), V(:,2), V(:,3));
alpha(0.8)
shading flat
view([-144 12])
axis equal
camlight left


band = sub2ind([ny,nx,nz], j,i,k);

toc
%% Extension of closest point
E = interp3_matrix(x1d, y1d, z1d, cpx, cpy, cpz, p, band);

% bdy cond;
% using the new "cp bar" definition; 
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

 
% surf(xg,yg,zg,u);
%% plot


xp = Vertices(:,1);
yp = Vertices(:,2);
zp = Vertices(:,3);
Eplot = interp3_matrix(x1d, y1d,z1d , xp(:), yp(:),zp(:), p, band);

up = Eplot*u;
up=reshape(up,size(xp));

figure(5); clf;
trisurf(Faces,xp,yp,zp,max(up,0))
material dull
camlight left
axis equal;
set(gcf,'color','w');
colormap('jet')
material dull
shading interp
axis equal;
hold on;
plot3(cpx(bdy), cpy(bdy), cpz(bdy), 'r.', 'linewidth', 3)


figure(5); clf
title('value of bdy')
trisurf(Faces,xp,yp,zp,Eplot*double(bdy))
hold on
plot3(cpx(bdy), cpy(bdy), cpz(bdy), 'r.', 'linewidth', 3)
shading flat
view([-144 12])
axis equal
