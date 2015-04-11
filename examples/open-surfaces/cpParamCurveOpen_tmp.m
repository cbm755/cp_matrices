function [cpx, cpy, dist, bdy] = cpParamCurveOpen_tmp(x,y,xs,ys,xp,yp,xpp,ypp,endpt)
    nx = length(x);     % number of points

    cpx = zeros(nx,1);
    cpy = zeros(nx,1);
    dist = zeros(nx,1);
    bdy = zeros(nx,1);
    param = zeros(nx,1);
    
    for pt = 1:nx
        [cpx(pt), cpy(pt), dist(pt), bdy(pt)] = cpParamCurveOpen_oldloop(x(pt),y(pt),xs,ys,xp,yp,xpp,ypp,endpt);
    end
end