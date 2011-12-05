function [varargout] = buildInterpPattern(X, relpt, dx, p)

% Helper Function for INTERP*_PATTERN
% Similar as buildInterpWeights, but do not compute weights.
  
% Xgrid is the B-pt (the base grid point below X, see figure in
% findGridInterpBasePt)
  [I, Xgrid] = findGridInterpBasePt(X, p, relpt, dx);

  dim = length(X);
  % 1D Stencil Width
  N = p+1;

  if (dim == 2)
    ii = I(1) + (0:(N-1));
    jj = I(2) + (0:(N-1));
    [iig, jjg] = meshgrid(ii, jj);
    iig = iig(:);
    jjg = jjg(:);
    varargout = {iig, jjg};

  elseif (dim == 3)
    ii = I(1) + (0:(N-1));
    jj = I(2) + (0:(N-1));
    kk = I(3) + (0:(N-1));
    [iig, jjg, kkg] = meshgrid(ii, jj, kk);
    iig = iig(:);
    jjg = jjg(:);
    kkg = kkg(:);
    varargout = {iig, jjg, kkg};

  else
    error(['Dimension ' num2str(dim) ' not yet implemented']);
  end
  
end