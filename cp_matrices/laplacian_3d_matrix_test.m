function L = laplacian_3d_matrix_test(x,y,z, order, band1, band2, time_dependent, dt, lambda)
%LAPLACIAN_3D_MATRIX  Build a 3D discrete Laplacian
%   ORDER: 2 or 4 for 2nd or 4th-order
%   Does no error checking up the equispaced nature of x,y,z
%
%   TODO: explain the roles of band1 and band2
%
%   To use ndgrid ordering pass "true" as the final argument

  if (nargin <= 5)
    band2 = band1;
	time_dependent = false;
  end
  
  % specifying time_dependent, but did not how much implicit is the scheme
  if (nargin == 7)
      if time_dependent == true
	  % corresponding to the total implicit scheme
      error('is time-dependent, but have not specified dt and lambda')
      end
  end

  % input checking
  [temp1, temp2] = size(x);
  if ~(  (ndims(x) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('x must be a vector, not e.g., meshgrid output');
  end
  [temp1, temp2] = size(y);
  if ~(  (ndims(y) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('y must be a vector, not e.g., meshgrid output');
  end
  [temp1, temp2] = size(z);
  if ~(  (ndims(z) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('z must be a vector, not e.g., meshgrid output');
  end

  dx = x(2)-x(1);
  dy = y(2)-y(1);
  dz = z(2)-z(1);
  if ~assertAlmostEqual([dx dx], [dy dz], 100*eps)
    error('this routine requires dx == dy == dz');
  end
  %ddx = [dx  dy  dz];
  %dim = length(ddx);
  Nx = round( (x(end)-x(1)) / dx ) + 1;
  Ny = round( (y(end)-y(1)) / dy ) + 1;
  Nz = round( (z(end)-z(1)) / dz ) + 1;
  %ptL = [x(1) y(1) z(1)];
  %ptH = [x(end) y(end) z(end)];

  if (order == 2)
    weights = [-6 1 1 1 1 1 1] / dx^2;
    stencil = [0 Ny -Ny 1 -1 Nx*Ny -Nx*Ny];
  elseif (order == 4)
    weights = [-15.0/2.0 ...
               (-1.0/12.0)  (4.0/3.0)  (4.0/3.0)  (-1.0/12.0) ...
               (-1.0/12.0)  (4.0/3.0)  (4.0/3.0)  (-1.0/12.0) ...
               (-1.0/12.0)  (4.0/3.0)  (4.0/3.0)  (-1.0/12.0) ...
              ] / dx^2;
    stencil = [0, -2*Ny,    -Ny,    Ny,    2*Ny, ...
                  -2,       -1,     1,     2,    ...
                  -2*Nx*Ny, -Nx*Ny, Nx*Ny, 2*Nx*Ny];
  else
    error(['order ' num2str(order) ' not implemented']);
  end

  if nargin >= 7 && time_dependent
	  weights =  - lambda*dt*weights;
      weights(1) = weights(1) + 1;
  end

  L = helper_diff_matrix3d_test(x, y, z, band1, band2, weights, stencil);
