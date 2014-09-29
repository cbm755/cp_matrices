function [cpx, dist, bdy] = cpInterval(x,interval)
    if nargin < 2
        interval = [-1,1];
    end
    a = interval(1);
    b = interval(2);
    mid = (a+b)/2;
    cpx = zeros(size(x));
    cpx(x<mid) = a;
    cpx(x>=mid) = b;
    bdy = x<a | x>b;
    dist = abs(cpx-x);
end