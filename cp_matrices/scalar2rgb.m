function RGB = scalar2rgb(u, caxes, colmap, int255)
%SCALAR2RGB  Lookup a list of scalars in a colormap
%  RGB = SCALAR2RGB(U) returns RGB for each u, one per row.
%  RGB = SCALAR2RGB(U, [umin,umax]) set a umin,umax.  u outside
%          this range will be clipped.
%  RGB = SCALAR2RGB(U, [], COLMAP)  use a custom colormap COLMAP.
%  RGB = SCALAR2RGB(U, [], [], true)  to return integer RGB in 0...255.
%
%  originally by Gustavo Chavez, some modifications by Colin Macdonald

%  todo: return array rather than separate r g b vectors?

if (nargin < 4)
  int255 = false;
end
if (  (nargin < 3)  ||  isempty(colmap)  )
  %colmap = colormap('default');  % causes an error
  colmap = colormap(jet);
end
if ( (nargin < 2)  ||  isempty(caxes) )
  u0 = min(u);
  u1 = max(u);
else
  u0 = caxes(1);
  u1 = caxes(2);
end

% scale to [0,1]
y = (u - u0) / (u1 - u0);
% clip
y = max( min(y,1), 0);

% length of colormap
N = size(colmap,1);

n = ceil(y*(N-1)+1);

RGB = colmap(n,:);

if int255
  RGB = round(RGB*255);
end