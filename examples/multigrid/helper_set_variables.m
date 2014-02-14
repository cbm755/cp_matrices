function [X1D, Y1D, XCP, YCP, BAND, Mc, Lc, Ec, V, F, A, BDYG] = ...
    helper_set_variables(x0, x1, y0, y1, dx, dx_coarsest, dim, p, order, rhsfn, cpf, has_boundary)

%% set 

% X,Y:     x,y coordinates defining the square domain

% XCP,YCP: closest points' coordinates of each level of V-cycle

% BAND:   store the innerbandfull (with respect to the whole domain square)
%         for the use of getting tranformation matrices (they are used to do 
%         prolongation and restriction between each neighbouring level of V-cycle,
%         see function 'helper_set_TM')

% number of vcycle levels:
n_level = round(log(dx_coarsest/dx)/log(2)) + 1;
% for judging which points need boundary condition modification
tol = 1e-10;

% allocate space for all the cells
DX = cell(n_level,1);
X1D = cell(n_level,1);
Y1D = cell(n_level,1);
BAND = cell(n_level,1);
XCP = cell(n_level,1);
YCP = cell(n_level,1);
Mc = cell(n_level,1);
Lc = cell(n_level,1);
Ec = cell(n_level,1);
V = cell(n_level,1);
F = cell(n_level,1);
A = cell(n_level,1);
BDYG = cell(n_level,1);


x1d = (x0:dx_coarsest:x1)';
y1d = (y0:dx_coarsest:y1)';

[xx yy] = meshgrid(x1d,y1d);
if has_boundary
    [cpx, cpy, dist, bdy] = cpf(xx,yy);
else
    [cpx, cpy, dist] = cpf(xx,yy);
    bdy = [];
end

bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
band = find(dist <= bw*dx_coarsest);

cpxg = cpx(band); cpyg = cpy(band);
distg = dist(band); 
if isempty(bdy)
    bdyg = [];
else
    bdyg = bdy(band);
end

E = interp2_matrix_band(x1d, y1d, cpxg, cpyg, p, band);
L = laplacian_2d_matrix(x1d, y1d, order, band, band);
M = lapsharp(L, E);

i = n_level;

X1D{i} = x1d;
Y1D{i} = y1d;
XCP{i} = cpxg;
YCP{i} = cpyg;
BAND{i} = band;
Mc{i} = M;
Lc{i} = L;
Ec{i} = E;
DX{i} = dx_coarsest;
DISTG{i} = distg;
BDYG{i} = bdyg;

[thg, rg] = cart2pol(cpxg,cpyg);
% in fact we do not need to set the rhs of each level of vcycle,
% but for purpose of debug, we do this here.
F{i} = rhsfn(thg);
V{i} = zeros(size(F{i}));

%    A{i} = helper_set_A(x1d,y1d,band);

for i = n_level-1:-1:1
    [BAND{i}, xg, yg, XCP{i}, YCP{i}, DISTG{i}, BDYG{i}, DX{i}, X1D{i}, Y1D{i}] = ...
        refine_grid(1, cpf, DX{i+1}, X1D{i+1}, Y1D{i+1}, bw, BAND{i+1}, DISTG{i+1}, BDYG{i+1});
    Ec{i} = interp2_matrix_band(X1D{i}, Y1D{i}, XCP{i}, YCP{i}, p, BAND{i});
    Lc{i} = laplacian_2d_matrix(X1D{i}, Y1D{i}, order, BAND{i}, BAND{i});
    Mc{i} = lapsharp(Lc{i}, Ec{i});
      
    [thg, rg] = cart2pol(XCP{i},YCP{i});
    % in fact we do not need to set the rhs of each level of vcycle,
    % but for purpose of debug, we do this here.
    F{i} = rhsfn(thg);
    V{i} = zeros(size(F{i}));
end

end

