function a = AreaPoly(poly)
T = diff([poly; poly(1,:)]);
h = poly(:,1) .* T(:,2) - poly(:,2) .* T(:,1);
a = .5 * sum(h);


