function A = findInBand(I, band, meshgridsz)
%FINDINDBAND  Find some given grid points in the band
%   A helper function for common closest point grid operators.
%
%   A = findInBand(I, band, Nx*Ny*Nz) Find the linear indices of the
%   meshgrid in the band.  Returns a matrix A that when multiplied by
%   u, picks out the those parts indexed by I.
%
%   In theory this is easy to do, but no so trivial to do quickly in
%   matlab.
%
%   TODO: is this generally useful?  what features should it have?
%   error checking?
%   TODO: should also return the "index permutation" so u(I) = A*u.

  % a fast vectorized approach, build the big matrix and then
  % discard columns
  Ei = (1:length(I))';
  Ej = I;
  A = sparse(Ei, Ej, ones(size(Ei)), length(I), meshgridsz);
  A = A(:,band);

  % a much slower approach:
  %Ei = (1:length(I))';
  %Ej = zeros(size(Ei));
  %Es = ones(size(Ei));
  %for i = 1:length(I)
  %  I2 = find(band == I(i));
  %  if isempty(I2)
  %    error('can''t find');
  %  end
  %  Ej(i) = I2;
  %end
  %B = sparse(Ei, Ej, Es, length(I), length(band));
