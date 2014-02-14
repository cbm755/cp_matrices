%% test geometric multigrid method to solve poisson equation on a sphere

%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../../cp_matrices');

% add functions for finding the closest points
addpath('../../surfaces');

% add notay amg
addpath('/scratch/cheny1/opt/AGMG_3.1.1/Matlab')

x0 = -3;
x1 = 3;

%%
% 2D example on a circle
% Construct a grid in the embedding space

dx = 0.4; % grid size
%dx = 0.05
%dx = 0.00625;  % lots of memory...
dx_coarsest = 0.4;   % coarsest grid size
x1d_coarsest = (x0:dx_coarsest:x1)';
y1d_coarsest = x1d_coarsest;
z1d_coarsest = x1d_coarsest;

dy = dx;
dz = dx;

dim = 3;  % dimension
p = 3;    % interpolation order
order = 2;  % Laplacian order: bw will need to increase if changed

bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

n1 = 3;
n2 = 3;

p_f2c = 1;
p_c2f = 1;

w = 1;

cpf = @cpSphere;

has_boundary = false;

%[a_x1d, a_y1d, a_xcp, a_ycp, a_band, Mc, Lc, Ec, V, F, A, a_bdyg] = ...
%    helper_set_variables(x0, x1, y0, y1, dx, dx_coarsest, dim, p, order, rhsfn, cpf, has_boundary);

disp('building cp grids ... ')
[a_band, a_xcp, a_ycp, a_zcp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_z1d, a_xg, a_yg, a_zg] = ...
    build_mg_cpgrid3d(x1d_coarsest, y1d_coarsest, z1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary);

n_level = length(a_band);

disp('building cp matrices ... ')
%[Mc, Lc, Ec] = build_mg_cpmatrix3d(a_band, a_xcp, a_ycp, a_zcp, a_x1d, a_y1d, a_z1d, p, order);
Mc = cell(n_level,1);
Lc = cell(n_level,1);
Ec = cell(n_level,1);

t0 = cputime;
MAX = 10000;
MAX1 = 100;
for cnt = 0.8*MAX:MAX

    alpha = cnt/MAX;

    for cnt1 = MAX1:MAX1
        beta = cnt1/MAX1 * (1-alpha); 
        gamma = 1 - alpha - beta;
        for i = 1:1:n_level
            ddx = a_x1d{i}(2) - a_x1d{i}(1);
            Ec{i} = interp3_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, a_xcp{i}, a_ycp{i}, a_zcp{i}, p);
            Ec{i} = Ec{i}(:, a_band{i});
            %Lc{i} = laplacian_3d_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, order, a_band{i}, a_band{i});
            Lc{i} = laplacian_wider_stencil_3d_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, order, alpha, beta, gamma, a_band{i}, a_band{i});
    
            E = interp3_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, a_xcp{i}, a_ycp{i}, a_zcp{i}, 1);
            E = E(:,a_band{i});
            gamma1 = 6 * ( alpha + beta/3 + gamma/2 );
            GAMMA = findGammaMatrix(alpha,gamma1,a_x1d{i},a_y1d{i},a_z1d{i},a_xcp{i},a_ycp{i},a_zcp{i},p);
            Mc{i} = E*Lc{i} - GAMMA*(speye(size(E))-Ec{i});
        end

        % test whether M-matrix property holds
        for i = 1:1:length(Mc)
            diagM = diag(diag(Mc{i}));
            % flag1 = diagM > 0;
            % nnz(flag1)
            tmp = Mc{i} - diagM;
            flag2 = tmp < 0;
            a = nnz(flag2)
        end

        if (a == 0)
            [alpha, beta, gamma] 
        end

     end
end
time = cputime - t0