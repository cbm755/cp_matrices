
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
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((stenrad+(p+1)/2)^2));
band = find(abs(dist) <= bw*dx);

% keep only stuff in the band
for d = 1:dim
  cpX{d} = cpX{d}(band);
  xg{d} = xx{d}(band);
end
dist = dist(band);

% TODO: just put in a cpgrid and call a banding routine?

cpgrid1.dim = 5;
cpgrid1.x1d = X1d;
cpgrid1.band = band;
cpgrid1.x = xg;
cpgrid1.cpx = cpX;
cpgrid1.dist = dist;
cpgrid1.cpfun = cpf;
cpgrid1.dx = dx;


disp('refining once');
cpgrid2 = refine_gridnd(cpgrid1, bw);

%disp('refining again');
cpgrid3 = refine_gridnd(cpgrid2, bw);

cpgrid = cpgrid2;


% TODO: make the matrix routinues support cell array
cpXtemp = [cpgrid.cpx{1:5}];

tic
L = laplacian_nd_matrix(cpgrid.x1d, 2, cpgrid.band);
toc

tic
%E1 = interpn_matrix(X1d, cpX, 1, band);
[Ei,Ej,Es] = interpn_matrix(cpgrid.x1d, cpXtemp, 1, cpgrid.band);
toc

tic
E1 = sparse(Ei, Ej, Es, length(cpgrid.band), length(cpgrid.band));
toc

tic
[Ei,Ej,Es] = interpn_matrix(cpgrid.x1d, cpXtemp, 2, cpgrid.band);
toc

tic
E2 = sparse(Ei, Ej, Es, length(cpgrid.band), length(cpgrid.band));
toc

disp('timing for matrix mult');
n = length(cpgrid.band)
u = rand(n,1);
tic; v=L*u; toc
tic; v=E1*u; toc
tic; v=E2*u; toc


