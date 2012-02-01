function [L, E, R, innerband, outerband, innerbandfull, outerbandfull, cpx, cpy, xg, yg, innerInOuter] = ...
      ops_and_bands2d_test(x1d,y1d, cpxinit,cpyinit, xginit,yginit, band_init, p, order)
% innerband/outerband: indices of the inner/outer band w.r.t. the initial
%              grid.
% innerbandfull/outerbandfull: indices of the inner/outer band w.r.t. the
%              full meshgrid square.

%% Construct an interpolation matrix for closest point
% This creates a matrix which interpolates data from the grid x1d y1d,
% onto the points cpx cpy.
disp('Constructing interpolation matrix')
tic;
[Etemp Ei Ej] = interp2_matrix_test(x1d, y1d, cpxinit, cpyinit, p);
toc;
disp('done')

%tic; [Y I] = sort(Ej(:,1)); toc;
%tic; [tmp I1] = unique(Ej(I,1)); toc;
%tic; M = Ej(I,:); toc;
%tic; iband_test = unique(M(I1,:)); toc;

%tic; iband_test2 = unique(Ej(:)); toc;

disp('finding innerband')
tic; Etemp1 = Etemp(:,band_init); toc;
tic; [i1,j1,S1] = find(Etemp1); toc;
tic; innerband = unique(j1); toc;
innerbandfull = band_init(innerband);
disp('done')

cpx = cpxinit(innerband);
cpy = cpyinit(innerband);
xg = xginit(innerband);
yg = yginit(innerband);

%% Create Laplacian matrix for heat equation
% in general, want the biggest finite difference stencil here so the
% others fit too
disp('Constructing laplacian matrix')
tic;
Ltemp = laplacian_2d_matrix_test(x1d,y1d, order, innerbandfull, band_init);
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

%% Restriction Operator
% constructing the innerband w.r.t. the outerband: innerInOuter, and then
% construct R used to extract inner values from an outer band vector.
disp('finding innerInOuter')
tic; Etemp2 = Etemp1(:,outerband); toc;
tic; [i2,j2,S2] = find(Etemp2); toc;
tic; innerInOuter = unique(j2);toc;
disp('done')

disp('building up the final returning operators')
tic;
li = length(innerband);
R = sparse((1:li)',innerInOuter,ones(li,1),li,length(outerband),li);

L = Ltemp(:, outerband);
E = Etemp(outerband, innerbandfull);
toc;
disp('done')

%clear Ltemp Etemp  % optional, erase the originals


%% Could instead regenerate everything, now that we have the two bands
%cpxin3 = cpx(innerbandfull); cpyin3 = cpy(innerbandfull);
%cpxout3 = cpx(outerbandfull); cpyout3 = cpy(outerbandfull);
%E3 = interp2_matrix_band(x1d,y1d, cpxout3, cpyout3, p, innerbandfull);
%L3 = laplacian_2d_matrix(x1d,y1d, order, innerbandfull, outerbandfull);


%% Outer band grid points and closest point
%cpxout = cpxinit(outerband);
%cpyout = cpyinit(outerband);
%xout = xinit(outerband);
%yout = yinit(outerband);

%cpxin = cpxinit(innerband);
%cpyin = cpyinit(innerband);
%xin = xinit(innerband);
%yin = yinit(innerband);

% alternatively, can use R to extract inner band from outer
%cpxin = R*cpxout;  cpyin = R*cpyout;
%xin = R*xout;  yin = R*yout;
