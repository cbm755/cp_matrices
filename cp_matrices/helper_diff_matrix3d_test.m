function L = helper_diff_matrix3d_test(x, y, z, band1, band2, weights, stencil)
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

  % TODO: a global variable is probably not a good way to do this.
  % global ICPM2009BANDINGCHECKS

  Nx = length(x);
  Ny = length(y);
  Nz = length(z);
  band1 = full(band1(:));
  
  StencilSize = length(weights);

  % Like a finite element code, find all the entries in 3 lists, then
  % insert them all at once while creating the sparse matrix
  Li = repmat((1:length(band1))',1,StencilSize);
  Lj = repmat(band1,1,StencilSize) + repmat(stencil,length(band1),1);
  Ls = repmat(weights,length(band1),1);

  L = sparse(Li, Lj, Ls, length(band1), Nx*Ny*Nz);

  % If we're using careful banding a la iCPM2009 then as a sanity
  % check all of the columns outside of band2 should be zero.
%   if (~isempty(ICPM2009BANDINGCHECKS)) && (ICPM2009BANDINGCHECKS)
%     Lout = L(:, setdiff(1:(Nx*Ny*Nz),band2));
%     if (nnz(Lout) > 0)
%       nnz(Lout)
%       error('Lost some non-zero coefficients (from outside the outerband)');
%     end
%   end

  % remove columns not in band2 (the outerband)
  L = L(:, band2);
  clear Li Lj Ls

end
