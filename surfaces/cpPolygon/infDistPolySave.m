function d = infDistPoly(x,y,poly)
d = 0*x;

m = size(poly,1);
% polyorig = poly;
% polyh = [poly ; poly(1,:)];
% poly = zeros(2*m,2);
% poly(2:2:2*m,:) = .5*(polyh(1:end-1,:) + polyh(2:end,:));
% poly(1:2:2*m-1,:) = polyorig;
% m = size(poly,1);


% tangents
T = diff([poly; poly(1,:)]);
h = sqrt(T(:,1).^2 + T(:,2).^2) * ones(1,2);
T = T./h;

% bisectors
B = .5 * (T + [T(end,:); T(1:end-1,:)]);
B = [-B(:,2) B(:,1)];


% closest poly corner
hx = [];
hy = [];
for k=1:m
    hx = [hx x(:)-poly(k,1)];
    hy = [hy y(:)-poly(k,2)];
end
h = hx.^2 + hy.^2;
[~,ind]=min(h,[],2);
indm1 = mod(ind-2,m)+1;


cc = [x(:)-poly(ind,1) , y(:)-poly(ind,2)];
alpha = cc(:,1) .* T(indm1,1) + cc(:,2) .* T(indm1,2);
beta = cc(:,1) .* T(ind,1) + cc(:,2) .* T(ind,2);
cpa = [cc(:,1) - alpha.* T(indm1,1) , cc(:,2) - alpha.* T(indm1,2)];
dcpa = sqrt( cpa(:,1).^2 + cpa(:,2).^2 );
cpb = [cc(:,1) - beta.* T(ind,1) , cc(:,2) - beta.* T(ind,2)];
dcpb = sqrt( cpb(:,1).^2 + cpb(:,2).^2 );

% detc = B(ind,1) .* cc(:,2) - B(ind,2) .* cc(:,1);
% indch = find(detc < 0);
% d(indch) = dcpb(indch);
% 
% indch = setdiff([1:length(x(:))],indch);
% d(indch) = dcpa(indch);
% 
% detTr = T(ind,1) .* cc(:,2) - T(ind,2) .* cc(:,1);
% detTl = T(indm1,1) .* cc(:,2) - T(indm1,2) .* cc(:,1);
% inds = find( (detTr < 0) & (detTl < 0) );
% d(inds) = -d(inds);

detc = B(ind,1) .* cc(:,2) - B(ind,2) .* cc(:,1);
indch = find(detc < 0);
d(indch) = dcpb(indch);

dets = T(ind,1) .* cc(:,2) - T(ind,2) .* cc(:,1);
inds = find( dets(indch) < 0 );
d(indch(inds)) = -d(indch(inds));


indch = setdiff([1:length(x(:))],indch);
d(indch) = dcpa(indch);

dets = T(indm1,1) .* cc(:,2) - T(indm1,2) .* cc(:,1);
inds = find( dets(indch) < 0 );
d(indch(inds)) = -d(indch(inds));


%test
% figure; axis([-1.5 1.5 -1.5 1.5])
% hold on
% polyc = [poly; poly(1,:)];
% plot(polyc(:,1),polyc(:,2),'r.-')
% quiver(poly(:,1),poly(:,2),T(:,1),T(:,2))
% quiver(poly(:,1),poly(:,2),B(:,1),B(:,2))
% hold off