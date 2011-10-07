function [varargout] = cpLineSegment(varargin)
%CPLINESEGMENT  Closest Point function for a line segment in 2D/3D
%   [cpx,cpy, dist, bdy] = cpLine(x,y,pt1,pt2)
%   [cpx,cpy,cpz, dist, bdy] = cpLine(x,y,z, pt1, pt2)
%   A line segment between 'pt1' and 'pt2'

  %% last input is a vector, use it determine dimension
  vec = varargin{nargin};
  dim = length(vec);

  if (nargin ~= dim + 2)
    error('wrong number of input arguments');
  end

  x = {};
  for j=1:dim
    x{j} = varargin{j};
  end
  p = varargin{dim+1};
  q = varargin{dim+2};

  if (dim == 2)
    [cp0_1,cp0_2,dist0,bdy0,s0] = cpLine(x{:}, q-p, p);
    cp0 = {cp0_1  cp0_2};
  elseif (dim == 3)
    [cp0_1,cp0_2,cp0_3,dist0,s0] = cpLine(x{:}, p-q, p);
    cp0 = {cp0_1  cp0_2  cp0_3};
  else
    error('dimension not implemented');
  end

  for j=1:dim
    cp{j} = ((s0 >= 0) & (s0 <= 1)) .* cp0{j} + ...
            (s0 < 0) .* (p(j)) + (s0 > 1) .* (q(j));
  end
  max(max(bdy0));
  bdy = ((s0 >= 0) & (s0 <= 1)) .* 0 + ...
        (s0 < 0) .* (1) + (s0 > 1) .* (2);

  %% outputs
  for j=1:dim
    varargout{j} = cp{j};
  end
  if (nargout > dim)
    %dist = norm(cp - x, 2)
    dist = zeros(size(x{1}));
    for j=1:dim
      dist = dist + (x{j} - cp{j}) .^ 2;
    end
    dist = sqrt(dist);
    varargout{dim+1} = dist;
  end
  if (nargout > dim + 1)
    varargout{dim+2} = bdy;
  end
