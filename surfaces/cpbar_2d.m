function [cpx,cpy, dist, bdy] = cpbar_2d(x,y, f, varargin)
%CPBAR_2D  Calculate the cpbar fcn used for boundary conditions.
%   "cpbar" refers to the double closest point projection described
%   in [Macdonald, Brandman, Ruuth, 2011] which is used to
%   implement higher-order accurate boundary conditions.
%
%   "f" is a handle to a cp function which takes x,y and whatever is
%   in varargin as inputs.  It must return "cpx,cpy, dist, bdy".  See
%   cpSemicircle.m for example.
%
%   Notes:
%     * dist refers to the distance to the original closest point,
%       not the cpbar point.  This means is safe to use this for
%       banding
%     * bdy is non-zero for points where the original closest point
%       hit a boundary
%
%   Code is vectorized: any size/shape for x should work, provided
%   function handle f is vectorized as well.


[cpx,cpy, dist, bdy] = f(x, y, varargin{:});

lbdy = logical(bdy);

% project into the interior
x1 = (2*cpx(lbdy) - x(lbdy));
y1 = (2*cpy(lbdy) - y(lbdy));

% recompute just those ones (it should be ok to do all of them, but
% this should be faster).
[cpx(lbdy),cpy(lbdy), tilde, tilde] = f(x1, y1, varargin{:});

