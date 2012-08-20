%example_nd_circle_in_4d()

X1d = {x1d y1d z1d w1d};
cpX = {cpx cpy cpz cpw};

cp0.x1d = X1d;
cp0.band = band;
cp0.x = {xg yg zg wg};
cp0.cpx = cpX;
cp0.dist = dist;
% dx dependecy here doesn't work (e.g., dx/2 doesn't do what I want).
% A random value probably safer.
cen = [0 0 0.1*rand 0.1*rand]
cp0.cpfun = @(x) cpCircleInHighDim(x, 1, cen);
cp0.dx = dx;

cp1 = refine_gridnd(cp0, bw)

cp2 = refine_gridnd(cp1, bw)

