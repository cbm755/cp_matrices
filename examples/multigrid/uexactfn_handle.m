function [u] = uexactfn_handle(th)
tol = 1e-10;
flag = th<0;
th(flag) = 2*pi+th(flag);
flag1 = th<pi;
flag2 = th>pi;
flag3 = abs(th-pi)<tol;
u = zeros(size(th));
u(flag1) = th(flag1);
u(flag2) = 2*pi-th(flag2);
u(flag3) = pi;
end
