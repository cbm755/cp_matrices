% TODO

dx = 0.25;

% column vector of the 1d grid
x1d = (-2.0:dx:2.0)';

% grid in each dimension, we're going to compute on the tensor product
X1d = {x1d x1d x1d x1d x1d};

tic
[xx{1:5}] = ndgrid(X1d{:});
toc

% closest point to (2D) sphere in 3D
cen = [0  0  0  rand*0.1  rand*0.1];
cpf = @(x) cpSphereInHighDim(x, 1, cen);

tic
[cpX,dist] = cpf(xx);
toc
%tic
%[cpx,cpy,cpz] = cpSphere(xx,yy,zz, 0.75);
%toc
%W1C = dx/3;
%W2C = dx/4;
%cpw1 = W1C*ones(size(cpx));
%cpw2 = W2C*ones(size(cpx));
%dist = sqrt( (xx-cpx).^2 + (yy-cpy).^2 + (zz-cpz).^2 + ...
%             (ww1-cpw1).^2 + (ww2-cpw2).^2);


% banding
dim = 5;
p = 3;       % max interpolation order
stenrad = 1; % max stencil radius for finite differences
bw = rm_bandwidth(dim, p, stenrad);
band = find(abs(dist) <= bw*dx);

% keep only stuff in the band
for d = 1:dim
  cpX{d} = cpX{d}(band);
  xg{d} = xx{d}(band);
end
dist = dist(band);

% TODO: just put in a cpgrid and call a banding routine?

g1.dim = dim;
g1.dx = dx;
g1.x1d = X1d;
g1.cpfun = cpf;
g1.band = band;
g1.x = xg;
g1.cpx = cpX;
g1.dist = dist;



disp('refining once');
g2 = refine_cpgrid_bw(g1, bw);

%disp('refining again');
%g3 = refine_cpgrid_bw(g2, bw);

g = g2;


tic
L = laplacian_nd_matrix(g.x1d, 2, g.band);
toc

tic
[Ei,Ej,Es] = interpn_matrix(g.x1d, g.cpx, 1, g.band);
toc

% TODO: can call directly now
%E1 = interpn_matrix(g.x1d, g.cpx, 1, g.band);

tic
E1 = sparse(Ei, Ej, Es, length(g.band), length(g.band));
toc

tic
[Ei,Ej,Es] = interpn_matrix(g.x1d, g.cpx, 2, g.band);
toc

tic
E2 = sparse(Ei, Ej, Es, length(g.band), length(g.band));
toc

disp('timing for matrix mult');
n = length(g.band)
u = rand(n,1);
tic; v=L*u; toc
tic; v=E1*u; toc
tic; v=E2*u; toc


