function [varargout] = refine_grid(varargin)
%REFINE_GRID   Make a finer CP representation
%   Given a banded cp grid, perform one level of grid refinement.
%   Works in 2D or 3D.
%
%   2D version:
%      [band, xg,yg, cpxg,cpyg, distg, x1d,y1d] = ...
%         refine_grid(cpf, dx0, x1d0, y1d0, bw, band0, dist0)
%      cpf: cp function handle.
%      dx0: grid spacing of current grid.
%      x1d0,y1d0: the 1d grids, these which form a "theoretical"
%                 meshgid() which is not created.
%      bw: bandwidth formula (will be multipled by dx), see
%          [Ruuth Merriman 2008].
%      band0: the band of the current grid.
%      dist0: distance field of the current grid.
%
%   3D version:
%      [band, x,y,z, cpx,cpy,cpz, dist, x1d,y1d,z1d] = ...
%         refine_grid(cpf, dx0, x1d0,y1d0,z1d0, bw, band0, dist0)
%
%   See "example_refine_grid.m".
%
%   TODO: add user toggle for ndgrid and test.

  %[nargin  nargout]
  varargout = {};
  for k=1:nargout
    varargout{k} = [];
  end

  if (nargin == 7 )
    [varargout{:}] = refine_grid2d(varargin{:});
  elseif (nargin == 8)
    [varargout{:}] = refine_grid3d(varargin{:});
  else
    error('incorrect calling');
  end
