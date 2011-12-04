function [L, E, R, innerband, outerband, innerbandfull, outerbandfull] = ...
      ops_and_bands2d(x1d,y1d, xinit,yinit, cpxinit,cpyinit, band_init, p, order)
% innerband/outerband: indices of the inner/outer band w.r.t. the initial
%              grid.
% innerbandfull/outerbandfull: indices of the inner/outer band w.r.t. the
%              full meshgrid square.

%% Construct an interpolation matrix for closest point
% This creates a matrix which interpolates data from the grid x1d y1d,
% onto the points cpx cpy.
disp('Constructing interpolation matrix');
% various alternatives, here we use the "full" E matrix...
% TODO: some subtly here: some entries of Etemp may be zero not by
% structure but just due to the particular interpolate weights.
% Maybe its safer to extract this information during construction
% of E (see python code).  OTOH, that might effect the numerical
% rank of the M...
Etemp = interp2_matrix(x1d, y1d, cpxinit, cpyinit, p);
[i,j,S] = find(Etemp);
innerbandfull = unique(j);
innerband = zeros(size(innerbandfull));
for i=1:length(innerbandfull)
  I = find(innerbandfull(i) == band_init);
  innerband(i) = I;
end


%% Create Laplacian matrix for heat equation
% in general, want the biggest finite difference stencil here so the
% others fit too
Ltemp = laplacian_2d_matrix(x1d,y1d, order, innerbandfull, band_init);


%% The outerband
% We find a narrow outerband by using the column-space of L.
tic; [i,j,S] = find(Ltemp); toc;
tic; outerband = unique(j); toc;
% indices are into band_init (the original columns of L),
% look them up to get the outerband in terms of the
% meshgrid(x1d,y1d) indices.
outerbandfull = band_init(outerband);


L = Ltemp(:, outerband);
E = Etemp(outerband, innerbandfull);
%clear Ltemp Etemp  % optional, erase the originals


%% Could instead regenerate everything, now that we have the two bands
%cpxin3 = cpx(innerbandfull); cpyin3 = cpy(innerbandfull);
%cpxout3 = cpx(outerbandfull); cpyout3 = cpy(outerbandfull);
%E3 = interp2_matrix_band(x1d,y1d, cpxout3, cpyout3, p, innerbandfull);
%L3 = laplacian_2d_matrix(x1d,y1d, order, innerbandfull, outerbandfull);


%% Restriction operator
% used to extract inner values from an outer band vector.  There
% is probably a slick loop-free way to do this.
innerInOuter = zeros(size(innerbandfull));
R = sparse([],[],[],length(innerbandfull),length(outerbandfull),length(innerbandfull));
for i=1:length(innerbandfull)
  I = find(outerbandfull == innerbandfull(i));
  innerInOuter(i) = I;
  R(i,I) = 1;
end

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
