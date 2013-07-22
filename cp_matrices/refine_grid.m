function [varargout] = refine_grid(varargin)
%REFINE_GRID   Make a finer CP representation
%   Given a banded cp grid, perform one or more steps of grid refinement.
%   Works in 2D or 3D.
%
%   DEPRECATED: use refine_cpgrid() and refine_cpgrid_bw() instead.
%
%   [cpGrid2] = refine_grid(M, cpGrid)
%      This call refines the grid specified by 'cpGrid' through 'M'
%      steps of refinement (each step halves the grid spacing).
%
%   Full 2D version:
%      [band, xg,yg, cpxg,cpyg, distg, dx, x1d,y1d] = ...
%         refine_grid(M, cpf, dx0, x1d0, y1d0, bw, band0, dist0)
%
%      [band, xg,yg, cpxg,cpyg, distg, bdyg, dx, x1d,y1d] = ...
%         refine_grid(M, cpf, dx0, x1d0, y1d0, bw, band0, dist0,
%         bdy0)
%
%      [band, xg,yg, cpxg,cpyg, distg, bdyg, x1d,y1d] = ...
%         refine_grid(M, cpf, dx0, x1d0, y1d0, bw, band0, dist0,
%         bdy0, use_ndgrid)
%
%      Inputs:
%         M: steps of refinement.
%
%         cpf: cp function handle.
%
%         dx0: grid spacing of current grid.  Must be same in x and y.
%
%         x1d0,y1d0: the 1d grids, these form a "theoretical"
%                    meshgid() (which is not created).
%
%         bw: bandwidth formula (will be multipled by dx), see
%            [Ruuth Merriman 2008].
%
%         band0: the band of the current grid.
%
%         dist0: distance field of the current grid.
%
%         bdy0: optional, for use with open surfaces.
%
%         use_ndgrid: optional, use ndgrid-style orderering instead of
%                     meshgrid.
%
%   3D version:
%      [band, x,y,z, cpx,cpy,cpz, dist, x1d,y1d,z1d] = ..4.
%         refine_grid(M, cpf, dx0, x1d0,y1d0,z1d0, bw, band0, dist0)
%
%   See "example_refine_grid.m".

  warning('Deprecated: see refine_cpgrid_bw(), refine_cpgrid()')

  %[nargin  nargout]
  %varargout = {};
  %for k=1:nargout
  %  varargout{k} = [];
  %end

  if (nargin == 2)
    cpf = cpGrid.cpf;
    dx = cpGrid.dx;
    dim = cpGrid.dim;
    x1d = cpGrid.x1d;
    y1d = cpGrid.y1d;
    if (dim == 3)
      z1d = cpGrid.z1d;
    end
    bw = cpGrid.bw;
    band = cpGrid.band;  % which one?
    dist = cpGrid.dist;
    bdy = cpGrid.bdy;  % mechanism for this?
    use_ndgrid = cpGrid.use_ndgrid;
    warning('cpGrid class functionality not implemeted yet');
  else
    [M, cpf, dx, x1d, y1d] = varargin{1:5};

    t1 = varargin{6};

    % Detect 2D/3D using the sixth input (this might fail if z1d has unit
    % length)
    if (length(t1) == 1)
      dim = 2;
      bw = t1;
      c = 7;
    else
      dim = 3;
      z1d = t1;
      bw = varargin{7};
      c = 8;
    end

    [band, dist] = varargin{c:(c+1)};
    if (nargin < c+2)
      bdy = [];
      withbdy = 0;
    else
      bdy = varargin{c+2};
      withbdy = 1;
    end
    if (nargin < c+3)
      use_ndgrid = 0;
    else
      use_ndgrid = varargin{c+3};
    end
  end


  %input = varargin(2:end);
  if (dim == 3)
    %[varargout{:}] = refine_grid3d(input{:});
    %[band, xg,yg,zg, cpxg,cpyg,cpzg, distg, bdyg, dx, x1d,y1d,z1d] = ...
    %    refine_grid3d(input{:});
    for k=1:M
      A = cputime();
      [band, xg,yg,zg, cpxg,cpyg,cpzg, dist, bdy, dx, x1d,y1d,z1d] = ...
          refine_grid3d(cpf, dx, x1d,y1d,z1d, bw, band, dist, bdy, use_ndgrid);
      fprintf('level %d (dx=%g) processing time: %g\n', k, dx, cputime()-A);
    end
    varargout = {band, xg,yg,zg, cpxg,cpyg,cpzg, dist, bdy, dx, x1d,y1d,z1d};
  else
    for k=1:M
      A = cputime();
      [band, xg,yg, cpxg,cpyg, dist, bdy, dx, x1d,y1d] = ...
          refine_grid2d(cpf, dx, x1d,y1d, bw, band, dist, bdy, use_ndgrid);
      if (M>1)
        fprintf('level %d (dx=%g) processing time: %g\n', k, dx, cputime()-A);
      else
        fprintf('dx=%g processing time: %g\n', dx, cputime()-A);
      end
    end
    if (withbdy)
      varargout = {band, xg,yg, cpxg,cpyg, dist, bdy, dx, x1d,y1d};
    else
      varargout = {band, xg,yg, cpxg,cpyg, dist, dx, x1d,y1d};
    end
  end
