function [xg, yg, dist] = get_proj_on_conormal_2d(x,y,cpx,cpy,nx,ny)

% TODO: documentation.
% Suppose we have an open curve. (x,y) has closest point (cpx,cpy) which 
% happens to be one of the endpoints, (nx,ny) is the conormal at that endpoint.
% We want to compute (xg,yg) as the projection of (x,y) on the conormal.

% The code is vectorized.

% shift (x,y) so that (cpx,cpy) is origin.
x = x - cpx;
y = y - cpy;

% calculate (xg,yg) in the shifted coord system.
inner_prod = x.*nx + y.*ny;
xg = inner_prod.*nx;
yg = inner_prod.*ny;

dist = sqrt((xg-x).^2 + (yg-y).^2);

% shift back
xg = xg + cpx;
yg = yg + cpy;

end