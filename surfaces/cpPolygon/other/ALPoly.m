function L = ALPoly(poly)
T = diff([poly; poly(1,:)]);
h = sqrt(T(:,1).^2 + T(:,2).^2);
L = sum(h);