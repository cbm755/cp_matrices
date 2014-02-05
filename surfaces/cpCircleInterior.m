function [cpx, cpy, dist, bdy] = cpCircleInterior(x, y, R, cen, inOutToggle)

%% find the closest points to a disk or exterior of the disk

% default
if (nargin < 3)
	R = 1;
end
if (nargin < 4)
	cen = [0,0];
end
if (nargin < 5)
	inOutToggle = 0;
end

% find closest point of a circle
x = x - cen(1);
y = y - cen(1);

[th, r] = cart2pol(x,y);
[cpx,cpy] = pol2cart(th,R);

sdist = r - R;

cpx = cpx + cen(1);
cpy = cpy + cen(2);

% find closest point of a disc
sdist = (-1)^(inOutToggle)*sdist;

% We want the boundary points to be "strictly" outside or inside the disk.
bdy = (sdist > eps);

cpx = bdy .* cpx + ~bdy .* x;
cpy = bdy .* cpy + ~bdy .* y;

dist = bdy .* sdist;

