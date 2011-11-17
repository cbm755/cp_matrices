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
  Nx = length(cp.x1d);
  Ny = length(cp.y1d);
  if (dim >= 3)
    Nz = length(cp.z1d);
  end
  if (dim < 2) | (dim > 3)
    warning('dimension not implemented');
  end

  if (dim == 2)
    meshgridsz = Nx*Ny;
    [j, i] = ind2sub([Ny Nx], band1);
    ip1 = i+1;  im1 = i-1;
    jp1 = j+1;  jm1 = j-1;
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
    meshgridsz = Nx*Ny*Nz;
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

  Mip1 = findInBand(bandip1, band2, meshgridsz);
  Mim1 = findInBand(bandim1, band2, meshgridsz);
  Mjp1 = findInBand(bandjp1, band2, meshgridsz);
  Mjm1 = findInBand(bandjm1, band2, meshgridsz);
  if (dim == 3)
    Mkp1 = findInBand(bandkp1, band2, meshgridsz);
    Mkm1 = findInBand(bandkm1, band2, meshgridsz);
  end


