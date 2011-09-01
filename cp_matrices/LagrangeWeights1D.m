function w = LagrangeWeights1D(xg, x, dx, N)
% 1D Barycentric Lagrange interpolation weights.
% xg is the base grid point, the left end grid point, B in the diagram
% in findGridInterpBasePt.m

  % Is exactly on a grid point, then return binary weights
  for j=0:(N-1)
    if (x == (xg+j*dx))
      w = zeros(1,N);
      w(j+1) = 1;
      return
    end
  end

  % I think its faster to hardcode than use nchoosek

  %# for i in range(0,N):
  %  #    w[i] = (-1)**i * comb(N-1,i)
  if     (N == 4), w = [1, -3, 3, -1];
  elseif (N == 5), w = [1, -4, 6, -4, 1];
  elseif (N == 1), w = 1;
  elseif (N == 2), w = [1, -1];
  elseif (N == 3), w = [1, -2, 1];
  elseif (N == 6), w = [1, -5, 10, -10, 5, -1];
  elseif (N == 7), w = [1, -6, 15, -20,  15,  -6,   1];
  elseif (N == 8), w = [1, -7, 21, -35,  35, -21,   7,  -1];
  else error('need to hardcode more weights');
  end

  for j=0:(N-1)
    w(j+1) = w(j+1) / ( x - (xg + j*dx) );
  end
  w = w / sum(w);
  return

end

