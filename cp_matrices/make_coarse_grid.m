function g = make_coarse_grid(dim, cpfun, dx, bb)
%MAKE_COARSE_GRID  Make a coarse grid, suitable for refinement
%   A convenience function to construct a cpgrid object
%   consisting of a simple meshgrid.  No banding is performed.
%
%   g = make_coarse_grid(dim, cpfun, dx)
%   g = make_coarse_grid(dim, cpfun, dx, bb)
%   The bounding box 'bb' of the grid defaults to [-2 2]^dim
%
%   Example:
%   g = make_coarse_grid(2, @cpCircle, 0.1, [-2 -2 2 2])
%   plot2d_compdomain
%
%   TODO: works for meshgrids in 2D and 3D, n-D is still TODO
%
%   TODO: a form like following would be convenient
%   g = make_coarse_grid(dim, cpfun, x1d, y1d, ...)?
%
%   TODO: should document which fields g can/should have
%     - invbandmap
%     - dx  (votes to delete)
%     - how to deal with inner/outer bands?

  if nargin < 4
    bb = [-2*ones(1,dim)  2*ones(1,dim)];
  end

  g.dim = dim;
  g.dx = dx;   % TODO: some votes for not keeping this in the
               % struct, instead maybe a helper function that
               % gives a scalar or vector of dx's.

  if length(dx) > 1
    warning('TODO: some work to do on refinement of nonequal dx')
    warning('TODO: this function currently also doesn''t support it')
  end

  if g.dim == 2
    g.x1d = bb(1):dx:bb(3);
    g.y1d = bb(2):dx:bb(4);
    [x y] = meshgrid(g.x1d, g.y1d);
    x = x(:);  y = y(:);
    g.cpfun = cpfun;
    g.band = (1:length(x))';
    g.x = x;
    g.y = y;
    [g.cpx, g.cpy, g.dist] = cpfun(x, y);
  elseif g.dim == 3
    g.x1d = bb(1):dx:bb(4);
    g.y1d = bb(2):dx:bb(5);
    g.z1d = bb(3):dx:bb(6);
    [x y z] = meshgrid(g.x1d, g.y1d, g.z1d);
    x = x(:);  y = y(:);  z = z(:);
    g.cpfun = cpfun;
    g.band = (1:length(x))';
    g.x = x;
    g.y = y;
    g.z = z;
    [g.cpx, g.cpy, g.cpz, g.dist] = cpfun(x, y, z);
  else
    error('todo: not implemented');
  end
