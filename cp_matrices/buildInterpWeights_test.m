function [varargout1, varargout2] = buildInterpWeights_test(X, relpt, dx, p, stencil)
%Helper function for INTERP*_MATRIX()
%
%   Build the 2D/3D weights to interpolate at a point X for a grid
%   measured relative to its lower-left corner "RELPT" with grid
%   spacing DX (which can be a vector, differing in each dimension)

  % Xgrid is the B-pt (the base grid point below X, see figure in
  % findGridInterpBasePt)
  [I, Xgrid] = findGridInterpBasePt_test(X, p, relpt, dx);

  npt = size(X,1);
  dim = size(X,2);
  % 1D Stencil Width
  N = p+1;

  if (dim == 2)
    xweights = LagrangeWeights1D_vec(Xgrid(:,1), X(:,1), dx(1), N);
    yweights = LagrangeWeights1D_vec(Xgrid(:,2), X(:,2), dx(2), N);
    iig = repmat(I(:,1),1,N^2) + repmat(stencil(1,:),npt,1);
    jjg = repmat(I(:,2),1,N^2) + repmat(stencil(2,:),npt,1);
    varargout1 = {xweights, yweights};
    varargout2 = {iig, jjg};

  elseif (dim == 3)
    xw = LagrangeWeights1D_vec(Xgrid(:,1), X(:,1), dx(1), N);
    yw = LagrangeWeights1D_vec(Xgrid(:,2), X(:,2), dx(2), N);
    zw = LagrangeWeights1D_vec(Xgrid(:,3), X(:,3), dx(3), N);
    iig = repmat(I(:,1),1,N^3) + repmat(stencil(1,:),npt,1);
    jjg = repmat(I(:,2),1,N^3) + repmat(stencil(2,:),npt,1);
    kkg = repmat(I(:,3),1,N^3) + repmat(stencil(3,:),npt,1);
    varargout1 = {xw, yw, zw};
    varargout2 = {iig, jjg, kkg};
  else
    error(['Dimension ' num2str(dim) ' not yet implemented']);
  end
  
end