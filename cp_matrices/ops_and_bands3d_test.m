function [L, E, R, innerband, outerband, innerbandfull, outerbandfull, cpx,cpy,cpz, xg,yg,zg, innerInOuter] = ...
      ops_and_bands3d_test(x1d,y1d,z1d, cpxinit,cpyinit,cpzinit, xinit, yinit, zinit, band_init, p, order)
% innerband/outerband: indices of the inner/outer band w.r.t. the initial
%              grid.
% innerbandfull/outerbandfull: indices of the inner/outer band w.r.t. the
%              full meshgrid square.

%% Construct an interpolation matrix for closest point
% This creates a matrix which interpolates data from the grid x1d y1d,
% onto the points cpx cpy.
disp('Constructing interpolation matrix')
tic;
Etemp = interp3_matrix_test(x1d, y1d, z1d, cpxinit, cpyinit, cpzinit, p);
toc;
disp('done')

disp('finding innerband')
Etemp1 = Etemp(:,band_init);
tic; [i1,j1,S1] = find(Etemp1); toc;
tic; innerband = unique(j1); toc;
innerbandfull = band_init(innerband);
disp('done')


%% Create Laplacian matrix for heat equation
% in general, want the biggest finite difference stencil here so the
% others fit too
disp('Constructing laplacian matrix')
tic;
Ltemp = laplacian_3d_matrix_test(x1d,y1d,z1d, order, innerbandfull, band_init);
toc;
disp('done')

%% The outerband
% We find a narrow outerband by using the column-space of L.
disp('finding outerbanb')
tic; [i,j,S] = find(Ltemp); toc;
tic; outerband = unique(j); toc;
% indices are into band_init (the original columns of L),
% look them up to get the outerband in terms of the
% meshgrid(x1d,y1d) indices.
outerbandfull = band_init(outerband);
disp('done')

disp('returning final L & E')
tic;
L = Ltemp(:, outerband);
E = Etemp(outerband, innerbandfull);
clear Etemp
toc;
disp('done')
%% Restriction Operator
% constructing the innerband w.r.t. the outerband: innerInOuter, and then
% construct R used to extract inner values from an outer band vector.
disp('finding innerInOuter')
tic; Etemp2 = Etemp1(:,outerband); toc;
tic; [i2,j2,S2] = find(Etemp2); toc;
tic; innerInOuter = unique(j2);toc;
tic; clear Etemp1 Etemp2; toc;
disp('done')


disp('building up the Restriction operators')
tic;
li = length(innerband);
R = sparse((1:li)',innerInOuter,ones(li,1),li,length(outerband),li);
toc;
disp('done')

xg = xinit(innerband);
yg = yinit(innerband);
zg = zinit(innerband);
cpx = cpxinit(innerband);
cpy = cpyinit(innerband);
cpz = cpzinit(innerband);

end
