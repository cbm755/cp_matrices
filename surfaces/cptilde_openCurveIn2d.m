function [cpx,cpy, dist, bdy] = cptilde_openCurveIn2d(x,y, f, n1, n2, varargin)
%CPTILDE_OPENCURVEIN2D  Calculate the cp_tilde fcn used for boundary conditions.
%
%   "f" is a handle to a cp function which takes x,y and whatever is
%   in varargin as inputs.  It must return "cpx,cpy, dist, bdy".  See
%   cpSemicircle.m for example.
%
%   "n1" is the unit normal vector at one of the end points of the curve;
%   "n2" is the unit normal vector at the other end point of the curve.
%   Both "n1" and "n2" should be 1-by-2 vectors.
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


x = x(:);
y = y(:);
[cpx,cpy, dist, bdy] = f(x, y, varargin{:});
cpxy = [cpx cpy];
xy = [x y];

lbdy = cell(2,1);
n = {n1,n2};
N = cell(2,1);
XY = cell(2,1);

for i = 1:2
lbdy{i} = logical(bdy == i);

% For the bdy points, find their mirror points to the normal vector;
% trying to compute mirror_x = 2*dot(x,n) - x.
N{i} = repmat(n{i},nnz(lbdy{i}),1);
% shift, make the first end points to be the origin.
XY{i} = xy(lbdy{i},:)-cpxy(lbdy{i},:);
dot_xy1_n1 = sum(XY{i}.*N{i},2);
XY{i} = 2*repmat(dot_xy1_n1,1,2).*N{i} - XY{i}; 
% shift back
XY{i} = XY{i} + cpxy(lbdy{i},:);
end

% recompute just those ones (it should be ok to do all of them, but
% this should be faster).
for i = 1:2
[cpx(lbdy{i}),cpy(lbdy{i}), disttemp, bdytemp] = f(XY{i}(:,1), XY{i}(:,2), varargin{:});
end

