function [Mip1 Mim1 Mjp1 Mjm1 Mkp1 Mkm1] = neighbourMatrices(cp, band1, band2)
%NEIGHBOURMATRICES  matrices to find neighbouring grid values
%   [E W N S] = neighbourMatrices(cp)
%      2D call, E=East=i+1, etc
%   [Mip1 Mim1 Mjp1 Mjm1 Mkp1 Mkm1] = neighbourMatrices(cp)
%      3D call
%
%   TODO: call the outputs E W N S U D?
%
%   TODO: WIP

  if (nargin < 2)
    % todo: iband / oband if they exist?
    %band1 = iband;
    %band2 = oband;
    band1 = cp.band;
    band2 = cp.band;
  end

  dim = cp.dim;
  %band = cp.band;

  if (dim >= 2)
    x1d = cp.x1d;
    y1d = cp.y1d;
    Nx = length(x1d);
    Ny = length(y1d);
  end
  if (dim >= 3)
    z1d = cp.z1d;
    Nz = length(z1d);
  end
  if (dim < 2) | (dim > 3)
    warning('dimension not implemented');
  end

  %dx = x1d(2) - x1d(1);
  %relpt = x1d(1);

  if (dim == 2)
    [j, i] = ind2sub([Ny Nx], band1);
    ip1 = i+1;  im1 = i-1;
    jp1 = j+1;  jm1 = j-1;
    ip2 = i+2;  jp2 = j+2;
    % im1 becomes negative, or ip1 too large, these could go outside the
    % meshgrid.  sub2ind() has checks for this and will give an error.
    % TODO: this could still be troublesome because it means relpt
    % and Nx might need to be larger than necessary for the band.
    % Could apply some limits here.
    bandip1 = sub2ind([Ny Nx], j,   ip1);
    bandim1 = sub2ind([Ny Nx], j,   im1);
    bandjp1 = sub2ind([Ny Nx], jp1, i);
    bandjm1 = sub2ind([Ny Nx], jm1, i);
  elseif ( dim==3 )
    [j, i, k] = ind2sub([Ny Nx Nz], band1);
    ip1 = i+1;  im1 = i-1;
    jp1 = j+1;  jm1 = j-1;
    kp1 = k+1;  km1 = k-1;
    bandip1 = sub2ind([Ny Nx Nz], j,   ip1, k);
    bandim1 = sub2ind([Ny Nx Nz], j,   im1, k);
    bandjp1 = sub2ind([Ny Nx Nz], jp1, i,   k);
    bandjm1 = sub2ind([Ny Nx Nz], jm1, i,   k);
    bandkp1 = sub2ind([Ny Nx Nz], j,   i,   kp1);
    bandkm1 = sub2ind([Ny Nx Nz], j,   i,   km1);
  end

  tic;
  nzmax = length(band1);
  Mip1 = sparse([], [], [], length(band1), length(band2), nzmax);
  Mim1 = sparse([], [], [], length(band1), length(band2), nzmax);
  Mjp1 = sparse([], [], [], length(band1), length(band2), nzmax);
  Mjm1 = sparse([], [], [], length(band1), length(band2), nzmax);
  if (dim == 3)
    Mkp1 = sparse([], [], [], length(band1), length(band2), nzmax);
    Mkm1 = sparse([], [], [], length(band1), length(band2), nzmax);
  end

  % TODO: might be faster to build the index permutation vectors
  % and then construct the matrices all-at-once.
  for c = 1:length(band1)
    %I = find(band == bandip1(c));
    %if ~isempty(I)
    %  Mip1(c, I) = 1;
    %end
    % if I is empty, it does the right thing
    I = find(band2 == bandip1(c));   Mip1(c, I) = 1;
    I = find(band2 == bandim1(c));   Mim1(c, I) = 1;
    I = find(band2 == bandjp1(c));   Mjp1(c, I) = 1;
    I = find(band2 == bandjm1(c));   Mjm1(c, I) = 1;
    if (dim == 3)
      I = find(band2 == bandkp1(c));   Mkp1(c, I) = 1;
      I = find(band2 == bandkm1(c));   Mkm1(c, I) = 1;
    end
  end
  toc

