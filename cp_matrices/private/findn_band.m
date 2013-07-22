function Ej = findn_band(xs, xi, p)
%FINDN_BAND  TODO
%   Ej = FINDN_BAND(X, XI, P)
%   Ej = FINDN_BAND({X Y Z ... W}, [XI YI ZI ... WI], P)
%   Ej = FINDN({X Y Z ... W}, {XI YI ZI ... WI}, P)
%   Build a minimal list of points required to perform degree P
%   interpolatation at the points in R^n, specified by XI, YI, ... WI.
%   P defaults to 3.  The points are indexed on the linear index of
%   the input points (ngrid order).
%
%   This is the same as the Ej (columns) of the interpolation matrix
%     [Ei,Ej,Es] = interpn_matrix(...)
%     Ej = unique(Ej);
%   except for a small savings in not computing Ei and Es (may 10%,
%   20%).  TODO: is it worth having this code then?

  if ~iscell(xs)
    error('expected a cell array of {x1d,y1d,etc}');
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

  Nsten = p+1;
  EXTSTENSZ = (Nsten)^dim;

  if iscell(xi)
    % convert xi to a matrix of column vectors.  TODO: is it better to
    % to change findGridInterpBasePt_vec to support cell array?
    xi2 = zeros(length(xi{1}), dim);
    for d=1:dim
      xi2(:,d) = xi{d};
    end
    xi = xi2;
  end

  Ni = length(xi(:,1));

  Ej = repmat((1:Ni)',1,EXTSTENSZ);

  % TODO: only need I, is this any faster?  (probably not much)
  %if (mod(p,2) == 0)  % even
  %  I = round( ( x - repmat(relpt,size(x,1),1) ) ./ repmat(dx,size(x,1),1) ) + 1  -  p/2;
  %else % odd
  %  I = floor( ( x - repmat(relpt,size(x,1),1) ) ./ repmat(dx,size(x,1),1) ) + 1  -  (p-1)/2;
  %end
  [Ibpt, Xgrid] = findGridInterpBasePt_vec(xi, p, ptL, ddx);

  %xw = {};
  %for d=1:dim
  %  xw{d} = LagrangeWeights1D_vec(Xgrid(:,d), xi(:,d), ddx(d), Nsten);
  %end

  NN = Nsten*ones(1,d);
  for s=1:prod(NN);
    [ii{1:dim}] = ind2sub(NN, s);
    for d=1:dim
      gi{d} = (Ibpt(:,d) + ii{d} - 1);
    end
    Ej(:,s) = sub2ind(Ns, gi{:});
  end
  T1 = cputime() - T1;
  Ej = unique(Ej(:));

end

