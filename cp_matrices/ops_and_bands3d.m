function [L, E, R, iband, oband, ibandfull, obandfull] = ...
      ops_and_bands3d(x1d,y1d,z1d, xinit,yinit,zinit, ...
                      cpxinit,cpyinit,cpzinit, band_init, p, order)
% iband/oband: indices of the inner/outer band w.r.t. the initial
%              grid.
% ibandfull/obandfull: indices of the inner/outer band w.r.t. the
%              full meshgrid cube.


%% Construct an interpolation matrix for closest point
% This creates a matrix which interpolates data from the grid x1d y1d,
% onto the points cpx cpy.
disp('Constructing interpolation matrix');
% various alternatives, here we use the "full" E matrix and find
% the support of the column space.
Etemp = interp3_matrix(x1d,y1d,z1d, cpxinit,cpyinit,cpzinit, p);
tic; [i,j,S] = find(Etemp); toc;
tic; ibandfull = unique(j); toc;
% find index of the inner band in the original band
iband = zeros(size(ibandfull));
for i=1:length(ibandfull)
  I = find(ibandfull(i) == band_init);
  iband(i) = I;
end


%% Create Laplacian matrix for heat equation
% in general, want the biggest finite difference stencil here so the
% others fit too
Ltemp = laplacian_3d_matrix(x1d,y1d,z1d, order, ibandfull, band_init);


%% The outerband
% We find a narrow outerband by using the column-space of L.
tic; [i,j,S] = find(Ltemp); toc;
tic; oband = unique(j); toc;
% indices are into band_init (the original columns of L),
% look them up to get the outerband in terms of the
% meshgrid(x1d,y1d) indices.
obandfull = band_init(oband);


%% we can extract the L and E matrices from the larger ones
L = Ltemp(:, oband);
E = Etemp(oband, ibandfull);
%clear Ltemp Etemp  % optional, erase the originals


%% Could instead regenerate everything, now that we have the two bands
%cpxin3 = cpx(ibandfull); cpyin3 = cpy(ibandfull); cpzin3 = cpz(ibandfull);
%cpxout3 = cpx(obandfull); cpyout3 = cpy(obandfull); cpzout3 = cpz(obandfull);
%E3 = interp3_matrix(x1d,y1d,z1d, cpxout3,cpyout3,cpzout3, p, ibandfull);
%L3 = laplacian_3d_matrix(x1d,y1d,z1d, order, ibandfull, obandfull);


%% Restriction operator
% used to extract inner values from an outer band vector.  There
% is probably a slick loop-free way to do this.
innerInOuter = zeros(size(ibandfull));
R = sparse([],[],[],length(ibandfull),length(obandfull),length(ibandfull));
for i=1:length(ibandfull)
  I = find(obandfull == ibandfull(i));
  innerInOuter(i) = I;
  R(i,I) = 1;
end


%% Outer band grid points and closest point
%cpxout = cpxinit(oband);
%cpyout = cpyinit(oband);
%cpzout = cpzinit(oband);
%xout = xinit(oband);
%yout = yinit(oband);
%zout = zinit(oband);

%cpxin = cpxinit(iband);
%cpyin = cpyinit(iband);
%cpzin = cpzinit(iband);
%xin = xinit(iband);
%yin = yinit(iband);
%zin = zinit(iband);

% alternatively, can use R to extract inner band from outer
%cpxin = R*cpxout;  cpyin = R*cpyout;  cpzin = R*cpzout;
%xin = R*xout;  yin = R*yout;  zin = R*zout;



