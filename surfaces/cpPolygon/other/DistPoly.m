function d = DistPoly(x,y,poly)
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
        
    ind = find( (0 <= alpha) & (alpha <= h(k)) );
    d(ind) = min(d(ind),dd(ind));
    
    indc = setdiff([1:length(d(:))]',ind);
    d(indc) = min(d(indc),hn(indc));
end


