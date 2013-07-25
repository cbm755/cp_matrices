function w = LagrangeWeights1D(xg, x, dx, N)
%LAGRANGEWEIGHTS1D   Barycentric Lagrange interpolation weights
%   w = LagrangeWeights1D(xg, x, dx, N)
%   1D Barycentric Lagrange interpolation weights.
%   xg is the base grid point, the left end grid point
%     (see also "B" in the diagram in findGridInterpBasePt.m)
%   x is the point of interpolation
%   dx is the grid spacing (uniform)
%   N is the number of points in the stencil (degree + 1)
%
%   Only works with scalar inputs, if you want to call this in a
%   loop, check whether LAGRANGEWEIGHTS1D_VEC would be better.

  % Is exactly on a grid point, then return binary weights
  for j=0:(N-1)
    if (x == (xg+j*dx))
      w = zeros(1,N);
      w(j+1) = 1;
      return
    end
  end

  % Calculate the basic weights: faster to hardcode than use nchoosek
  switch N
    case 1, w = 1;
    case 2, w = [1, -1];
    case 3, w = [1, -2, 1];
    case 4, w = [1, -3, 3, -1];
    case 5, w = [1, -4, 6, -4, 1];
    case 6, w = [1, -5, 10, -10, 5, -1];
    case 7, w = [1, -6, 15, -20,  15,  -6,   1];
    case 8, w = [1, -7, 21, -35,  35, -21,   7,  -1];
    otherwise
      warning('LagrangeWeights: should hardcode more weights');
      w = zeros(1,N);
      for i=1:N
        w(i) = (-1)^(i-1) * nchoosek(N-1,i-1);
      end
  end

  %for j=0:(N-1)
  %  w(j+1) = w(j+1) / ( x - (xg + j*dx) );
  %end

  % next two lines should be same, but in Octave a:b is treated
  % differently somehow than [a:b] and this makes my unit test
  % against LagrangeWeights1D_vec fail.
  %w = w ./ ( x - (xg + (0:(N-1))*dx) );
  w = w ./ ( x - (xg + [0:(N-1)]*dx) );
  w = w / sum(w);
end

