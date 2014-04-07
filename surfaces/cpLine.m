function [varargout] = cpLine(varargin)
%CPLINE  Closest Point function for a line in 2D or 3D
%   [cpx,cpy, dist, bdy] = cpLine(x,y,dir)
%   [cpx,cpy, dist, bdy] = cpLine(x,y,dir,pt)
%   [cpx,cpy,cpz, dist, bdy] = ...
%   [cpx,cpy,cpz, dist, bdy, s] = ...
%   'pt' specifies any point on the line and 'dir' specifies the
%   tangential direction.
%
%   Outputs: 'bdy' will always be 0.  's' gives the
%   parameter value measuring distance of the closest point from pt
%   in direction 'dir' (measured in units of magnitude 'dir').


  %% last input is a vector, use it determine dimension
  vec = varargin{nargin};
  dim = length(vec);

  if ~ ((nargin == dim + 1) | (nargin == dim + 2))
    error('wrong number of input arguments');
  end

  %% other inputs
  x = {};
  for j=1:dim
    x{j} = varargin{j};
  end

  dir = varargin{dim+1};
  if (norm(dir) < 100*eps)
    error('direction has norm almost zero')
  end
  if (nargin == dim+2) % there is a point
    p = varargin{dim+2};
  else
    p = zeros(size(dir));
  end

  %% Project
  %s = dot(x - p,  dir) / norm(dir, 2)^2;
  s = zeros(size(x{1}));
  for j=1:dim
    s = s + ( x{j} - p(j) ) .* dir(j);
  end
  s = s / (norm(dir,2)^2);

  %% closest points
  cp = {};
  for j=1:dim
    cp{j} = p(j) + dir(j)*s;
  end

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
    bdy = zeros(size(x{1}));
    varargout{dim+2} = bdy;
  end
  if (nargout > dim + 2)
    param = s;
    varargout{dim+3} = s;
  end
