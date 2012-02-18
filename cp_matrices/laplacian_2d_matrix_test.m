function L = laplacian_2d_matrix_test(x,y, order, band1, band2)
%LAPLACIAN_2D_MATRIX  Build a 2D discrete Laplacian
%   ORDER: 2 or 4 for 2nd or 4th-order
%   Does no error checking up the equispaced nature of x,y,z
%
%   TODO: explain the roles of band1 and band2
%
%   TODO: currently assumes dx=dy, but should be an easy fix
%
%   TODO: issue with two bands: currently need extra padding in
%   wherever we call this from: should fix this
%
%   To use ndgrid ordering pass "true" as the final argument

  % input checking
  [temp1, temp2] = size(x);
  if ~(  (ndims(x) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('x must be a vector, not e.g., meshgrid output');
  end
  [temp1, temp2] = size(y);
  if ~(  (ndims(y) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('y must be a vector, not e.g., meshgrid output');
  end

  dx = x(2)-x(1);
  dy = y(2)-y(1);
  if ~assertAlmostEqual(dx, dy, 100*eps)
    error('this routine requires dx == dy');
  end
  %ddx = [dx  dy];
  %dim = length(ddx);
  Nx = round( (x(end)-x(1)) / dx ) + 1;
  Ny = round( (y(end)-y(1)) / dy ) + 1;
  %ptL = [x(1) y(1)];
  %ptH = [x(end) y(end)];

  %warning('only works for meshgrid, nd_grid not implemented')
  if (order == 2)
    weights = [-4 1 1 1 1] / dx^2;
    stencil = [0  Ny -Ny 1 -1];
    %weights = [1 1 -4 1 1] / dx^2;
    %stencil = [-Ny -1 0 1 Ny];
     
  elseif (order == 4)
    weights = [-5.0 ...
               (-1.0/12.0)  (4.0/3.0)  (4.0/3.0)  (-1.0/12.0) ...
               (-1.0/12.0)  (4.0/3.0)  (4.0/3.0)  (-1.0/12.0) ...
              ] / dx^2;
    stencil = [0 -2*Ny -Ny Ny 2*Ny -2 -1 1 2]; 
  else
    error(['order ' num2str(order) ' not implemented']);
  end

  L = helper_diff_matrix2d_test(x, y, band1, band2, weights, stencil);
