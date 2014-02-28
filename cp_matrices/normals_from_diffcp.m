function [varargout] = normals_from_diffcp(g)
%NORMALS_FROM_DIFFCP   normals from differentiation of closest points
%   Assuming g is a grid structure for a codim-1 surface in 3D:
%      [n,t1,t2] = normals_from_diffcp(g);
%   returns the normal and two tangent plane vectors
%
%   TODO: would be good if this supported using x - cp for points
%   far from the surface.
%
%   TODO: documentation.

%varargout = {1:nargout};
  if g.dim == 3
    [varargout{1:nargout}] = normals_from_diffcp3d(g);
  elseif g.dim == 2
    [varargout{1:nargout}] = normals_from_diffcp2d(g);
  else
    error('dimension not implemented');
  end