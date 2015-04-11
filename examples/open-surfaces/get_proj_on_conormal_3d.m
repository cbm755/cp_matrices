function [xg, yg, zg, dist] = get_proj_on_conormal_3d(x,y,z,cpx,cpy,cpz,nx,ny,nz)

% TODO: documentation.
% Suppose we have an open curve. (x,y) has closest point (cpx,cpy) which 
% happens to be one of the endpoints, (nx,ny) is the conormal at that endpoint.
% We want to compute (xg,yg) as the projection of (x,y) on the conormal.

% The code is vectorized.

% shift (x,y,z) so that (cpx,cpy,cpz) is origin.
x = x - cpx;
y = y - cpy;
z = z - cpz;

% calculate (xg,yg,zg) in the shifted coord system.
inner_prod = x.*nx + y.*ny + z.*nz;
xg = inner_prod.*nx;
yg = inner_prod.*ny;
zg = inner_prod.*nz;

dist = sqrt((xg-x).^2 + (yg-y).^2 + (zg-z).^2);

% shift back
xg = xg + cpx;
yg = yg + cpy;
zg = zg + cpz;

end