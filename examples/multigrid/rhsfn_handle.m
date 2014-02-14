function [f] = rhsfn_handle(th,r,n)

rhsfn_p = @(th,r) ( n*(n-1)*sin(th).^(n-2).*cos(th).^2 - n*sin(th).^n )./ (r.^2);
f = rhsfn_p(th,r);
flag = th<0;
f(flag) = -f(flag);
%uexactfn = @(th) abs(sin(th)).^n;

end
