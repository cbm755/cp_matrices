function d = helper_sDistPoly(x, y, poly)
%helper function, use cpPolygon instead

d = Inf*ones(size(x));
s = ones(size(x));
%sig = sign(AreaPoly(poly));

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
    dd = -hx .* T(k,2) + hy .* T(k,1);
    dda = abs(dd);

    ind = find( (0 <= alpha) & (alpha <= h(k)) );
    i = find( dda(ind) < d(ind) );
    d(ind(i)) = dda(ind(i));
    s(ind(i)) = sign(dd(ind(i)));

    indc = setdiff([1:length(d(:))]',ind);
    i = find( hn(indc) < d(indc) );
    d(indc(i)) = hn(indc(i));

    km1 = mod(k-2,m)+1;
    W = [-T(k,2) -T(km1,2) ; T(k,1) T(km1,1)];
    Winv = W^(-1);
    c1 = Winv(1,1) * hx + Winv(1,2) *hy;
    s(indc(i)) = sign( c1(indc(i)) );
end
d = s.*d;