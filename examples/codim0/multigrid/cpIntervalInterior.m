function [cpx, dist, bdy] = cpIntervalInterior(x,interval)
    if nargin < 2
        interval = [-1,1];
    end
    a = interval(1);
    b = interval(2);
    cpx = x;
    cpx(x<a) = a;
    cpx(x>b) = b;
    bdy = x<=a+10*eps | x>=b-10*eps;
    dist = abs(cpx-x);
end
    
   