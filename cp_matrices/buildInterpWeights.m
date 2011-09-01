function [weights, varargout] = buildInterpWeights(X, relpt, dx, p)
%Helper function for INTERP*_MATRIX()
%
%   Build the 2D/3D weights to interpolate at a point X for a grid
%   measured relative to its lower-left corner "RELPT" with grid
%   spacing DX (which can be a vector, differing in each dimension)

  % Xgrid is the B-pt (the base grid point below X, see figure in
  % findGridInterpBasePt)
  [I, Xgrid] = findGridInterpBasePt(X, p, relpt, dx);

  dim = length(X);
  % 1D Stencil Width
  N = p+1;

  if (dim == 2)
    xweights = LagrangeWeights1D(Xgrid(1), X(1), dx(1), N);
    yweights = LagrangeWeights1D(Xgrid(2), X(2), dx(2), N);
    ii = I(1) + (0:(N-1));
    jj = I(2) + (0:(N-1));
    [iig, jjg] = meshgrid(ii, jj);
    [xweights, yweights] = meshgrid(xweights, yweights);
    weights = xweights .* yweights;
    weights = weights(:);
    iig = iig(:);
    jjg = jjg(:);
    varargout = {iig, jjg};

  elseif (dim == 3)
    xw = LagrangeWeights1D(Xgrid(1), X(1), dx(1), N);
    yw = LagrangeWeights1D(Xgrid(2), X(2), dx(2), N);
    zw = LagrangeWeights1D(Xgrid(3), X(3), dx(3), N);
    ii = I(1) + (0:(N-1));
    jj = I(2) + (0:(N-1));
    kk = I(3) + (0:(N-1));
    [iig, jjg, kkg] = meshgrid(ii, jj, kk);
    [xw, yw, zw] = meshgrid(xw, yw, zw);
    weights = xw .* yw .* zw;
    weights = weights(:);
    iig = iig(:);
    jjg = jjg(:);
    kkg = kkg(:);
    varargout = {iig, jjg, kkg};

  else
    error(['Dimension ' num2str(dim) ' not yet implemented']);
  end

  sumw = sum(weights);
  if ( abs(sumw - 1.0) > 32*eps)
    err = sumw - 1.0
    error('interp weight problem, didn''t sum close enough to 1');
  end
