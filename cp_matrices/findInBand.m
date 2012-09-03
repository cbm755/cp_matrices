function A = findInBand(I, band, meshgridsz)
%FINDINBAND  Find some given grid points in the band
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
%   TODO: could also return the "index permutation" so u(I) = A*u.


  SaveMemory = true;

  % A fast vectorized approach: build big matrix then discard columns.
  % Potentially uses a lot of memory, sometimes is slightly faster.
  if ~SaveMemory
    Ei = (1:length(I))';
    Ej = I;
    A = sparse(Ei, Ej, ones(size(Ei)), length(I), meshgridsz);
    A = A(:,band);
  end


  % This approach uses less memory but might be slightly slower in
  % some cases.
  if SaveMemory
    % TODO: could cache this, store inside the cpgrid object etc.
    invbandmap = make_invbandmap(meshgridsz, band);
    Ei = (1:length(I))';
    Ej = invbandmap(I);
    nz = Ej ~= 0;
    Ei = Ei(nz);
    Ej = Ej(nz);
    A = sparse(Ei, Ej, 1, length(I), length(band));
  end


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
