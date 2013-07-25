function ww = LagrangeWeights1D_vec(xg, x, dx, N)
%LAGRANGEWEIGHTS1D_VEC   Barycentric Lagrange interpolation weights
%   w = LagrangeWeights1D_vec(xg, x, dx, N)
%   1D Barycentric Lagrange interpolation weights.
%   xg is the base grid point, the left end grid point
%     (see also "B" in the diagram in findGridInterpBasePt.m)
%   x is the point of interpolation
%   dx is the (uniform) grid spacing
%   N is the number of points in the stencil (degree + 1)
%
%   This is the vectorized version: xg, x, dx can be column
%   vectors.  N must be a scalar.  If you want to do lots of scalar
%   calls in a loop, than LAGRANGEWEIGHTS1D is faster.
%
%   Notes: this function will very likely divide by zero but this
%   is intended.  In Octave, you can disable the warning with
%     warning('off', 'Octave:divide-by-zero')

% If this is the first time you've looked at this, its rather amazing
% that Barycentric Lagrange works in practice.  But it does: see the
% literature.

%if exist('octave_config_info', 'builtin')
%    % if octave, then temporarily ensure div-by-zero warning is disabled
%    q = warning('query', 'Octave:divide-by-zero');
%    if strcmp(q.state, 'on')
%      warning('off', 'Octave:divide-by-zero')
%    end
%  end


  %% Calculate the basic weights
  % Faster to hardcode than use nchoosek (although doesn't matter
  % much if you're using this vectorized rather then in a loop)
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
      w = zeros(1,N);
      for i=1:N
        w(i) = (-1)^(i-1) * nchoosek(N-1,i-1);
      end
  end

  M = length(x);

  %% compute the interp weights assuming no div by zero
  denom = zeros(M, N);  % we want this later to detect DbZ
  ww = zeros(M, N);
  for j=1:N
    denom(:,j) =  x - (xg + (j-1)*dx);
    ww(:,j) = w(j) ./ denom(:,j);
  end
  % repmat version, (although the loop for denom is much better)
  %vec = 0:(N-1);
  %denom = x(:,ones(N,1)) - (xg(:,ones(N,1)) + vec(ones(size(x)),:)*dx);
  %ww = repmat(w, M,1);
  %ww = ww ./ denom;


  %% find infs from division by zero
  infmask = (denom == 0);
  rows = any(infmask,2);
  ww(rows,:) = infmask(rows,:);


  %% normalize
  rowsums = sum(ww,2);
  for j=1:N
    ww(:,j) = ww(:,j) ./ rowsums;
  end
  % a repmat approach is roughly twice as slow
  %ww = ww ./ repmat(sum(ww,2), 1, N);

  %  if exist('octave_config_info', 'builtin')
  %  if strcmp(q.state, 'on')
  %    warning('on', 'Octave:divide-by-zero')
  %  end
  %end
end

