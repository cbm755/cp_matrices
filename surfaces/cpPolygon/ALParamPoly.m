function [gx,gy]=ALParamPoly(s,poly)
L = ALPoly(poly);
s = mod(s,L);
gx = s;
gy = s;
T = diff([poly; poly(1,:)]);
h = sqrt(T(:,1).^2 + T(:,2).^2);
T = T./(h*ones(1,2));

m = length(h);
L1 = 0;
for k=1:m
    L2 = L1 + h(k);
    ind = find( (L1 <= s) & (s <= L2));
    gx(ind) = poly(k,1) + (s(ind)-L1) * T(k,1);
    gy(ind) = poly(k,2) + (s(ind)-L1) * T(k,2);
    L1 = L2;
end

% %test
% figure; axis([-1.5 1.5 -1.5 1.5])
% hold on
% polyc = [poly; poly(1,:)];
% plot(polyc(:,1),polyc(:,2),'r.-')
% quiver(poly(:,1),poly(:,2),T(:,1),T(:,2))
% plot(gx,gy)
% hold off