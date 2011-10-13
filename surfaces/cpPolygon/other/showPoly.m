function showPoly(poly)
polyc = [poly; poly(1,:)];

% tangents
T = diff(polyc);
h = sqrt(T(:,1).^2 + T(:,2).^2);
T = T./(h * ones(1,2));

% bisectors
B = .5 * (T + [T(end,:); T(1:end-1,:)]);
B = [-B(:,2) B(:,1)];

hold on
plot(polyc(:,1),polyc(:,2),'r.-')
quiver(poly(:,1),poly(:,2),T(:,1),T(:,2))
quiver(poly(:,1),poly(:,2),B(:,1),B(:,2))
plot(poly(1,1),poly(1,2),'g*')
hold off