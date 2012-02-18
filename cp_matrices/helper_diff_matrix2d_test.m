function L = helper_diff_matrix2d_test(x, y, band1, band2, weights, stencil)
%HELPER function: not intended for external use
%
% use_ndgrid: false for meshgrid, true for ndgrid
%
%   TODO: explain the roles of band1 and band2
%
%   TODO: this should work with dx != dy
%
%   TODO: issue with two bands: currently need extra padding in
%   wherever we call this from: should fix this

  % TODO: global variable may not be the best way to do this
  global ICPM2009BANDINGCHECKS

  Nx = length(x);
  Ny = length(y);
  band1 = full(band1(:));
  
  StencilSize = length(weights);

  
  %tic
  % Like a finite element code, find all the entries in 3 lists, then
  % insert them all at once while creating the sparse matrix
  
  Li = repmat((1:length(band1))',1,StencilSize);
  Lj = repmat(band1,1,StencilSize) + repmat(stencil,length(band1),1);
  Ls = repmat(weights,length(band1),1);
   
  
  Li_vec = reshape(Li, 1, length(band1)*StencilSize);
  Lj_vec = reshape(Lj, 1, length(band1)*StencilSize);
  Ls_vec = reshape(Ls, 1, length(band1)*StencilSize);
  
  L = sparse(Li_vec(:), Lj_vec(:), Ls_vec(:), length(band1), Nx*Ny);

  % If we're using careful banding a la iCPM2009 then as a sanity
  % check all of the columns outside of band2 should be zero.
  % (This won't be true when using just one band and its probably
  % not true in the dual-band but using bandwidth estimates)
%   if (~isempty(ICPM2009BANDINGCHECKS)) && (ICPM2009BANDINGCHECKS)
%     Lout = L(:,setdiff(1:(Nx*Ny),band2));
%     if (nnz(Lout) > 0)
%       nnz(Lout)
%       error('Lost some non-zero coefficients (from outside the outerband)');
%     end
%   end

  % remove columns not in band2 (the outerband)
  L = L(:,band2);
  
  %Ltime = toc;
end
