function [cpx, cpy] = helper_CPopPoly(x, y, poly)
%helper function, use cpPolygon instead

cpx = 0*x;
cpy = 0*y;
d = Inf*ones(size(x));

m = size(poly,1);
poly = [poly; poly(1,:)];

% tangents
T = diff(poly);
h = sqrt(T(:,1).^2 + T(:,2).^2);
T = T./(h * ones(1,2));

for k=1:m
    hx = x(:)-poly(k,1);
    hy = y(:)-poly(k,2);
    hn = sqrt(hx.^2 + hy.^2);

    alpha = hx .* T(k,1) + hy .* T(k,2);
    dd = abs(- hx .* T(k,2) + hy .* T(k,1));

    candx = poly(k,1) + alpha * T(k,1);
    candy = poly(k,2) + alpha * T(k,2);

    ind = find( (0 <= alpha) & (alpha <= h(k)) );
    i = find( dd(ind) < d(ind) );
    d(ind(i)) = dd(ind(i));
    cpx(ind(i)) = candx(ind(i));
    cpy(ind(i)) = candy(ind(i));


    indc = setdiff([1:length(d(:))]',ind);
    i = find( hn(indc) < d(indc) );
    d(indc(i)) = hn(indc(i));
    cpx(indc(i)) = poly(k,1);
    cpy(indc(i)) = poly(k,2);
end



% %test
% figure; axis([-1.5 1.5 -1.5 1.5])
% hold on
% polyc = [poly; poly(1,:)];
% plot(polyc(:,1),polyc(:,2),'r.-')
% quiver(poly(:,1),poly(:,2),T(:,1),T(:,2))
% quiver(poly(:,1),poly(:,2),B(:,1),B(:,2))
% hold off