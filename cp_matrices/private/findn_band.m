function Ej = findn_band(xs, xi, p)
%FINDN_BAND  TODO
%   Ej = FINDN_BAND(X, XI, P)
%   Ej = FINDN_BAND({x y z ... w}, {xi yi zi ... wi}, P)
%   Build a minimal list of points required to perform degree P
%   interpolatation at the points in R^n, specified by xi, yi, ... wi.
%   P defaults to 3.  The points are indexed on the linear index of
%   the input points (ngrid order).
%
%   This is the same as the Ej (columns) of the interpolation matrix
%     [Ei,Ej,Es] = interpn_matrix(...)
%     Ej = unique(Ej);
%   except for a small savings in not computing Ei and Es (may 10%,
%   20%).  TODO: is it worth having this code then?

  if ~iscell(xs)
    error('expected a cell array of {x1d,y1d,...}');
  end
  if ~iscell(xi)
    error('expected a cell array of {cpx,cpy,...}');
  end

  if (nargin < 3)
    p = [];
  end

  if (isempty(p))
    p = 3;  % default interp degree
  end

  T1 = cputime();
  dim = length(xs);
  Ns = zeros(1, dim);
  ddx = zeros(1, dim);
  ptL = zeros(1, dim);
  for d=1:dim
    Ns(d) = length(xs{d});
    ddx(d) = xs{d}(2)-xs{d}(1);
    ptL(d) = xs{d}(1);
  end
  M = prod(Ns);
  if (M > 1e15)
    error('too big to use doubles as indicies: implement int64 indexing')
  end

  % TODO: only need I, faster to copy-paste a bit of it?  Probably not
  [Ibpt, tilde] = findGridInterpBasePt_vec(xi, p, ptL, ddx);

  %xw = {};
  %for d=1:dim
  %  xw{d} = LagrangeWeights1D_vec(Xgrid(:,d), xi(:,d), ddx(d), Nsten);
  %end

  % in my 4D tests, putting "unique" inside the loop was about 5 times faster.
  Ej = [];
  NN = (p+1)*ones(1,d);
  for s=1:prod(NN);
    [ii{1:dim}] = ind2sub(NN, s);
    for d=1:dim
      gi{d} = (Ibpt{d} + ii{d} - 1);
    end
    si = sub2ind(Ns, gi{:});
    Ej = unique([Ej; si(:)]);
  end
  T1 = cputime() - T1;
  %Ej = unique(Ej(:));

end

