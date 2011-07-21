function [cpx, cpy, dist, varargout] = cpParamCurve_2D(x,y,xs,ys,xp,yp,xpp,ypp,endpt1,endpt2,isclosed)
% matlab version of function in python code ParamCurve
% Find closest points of parameterised curves using Newton's method

% x,y : point to find closest point for
% xs,ys : parameterised curve
% xp,yp : first derivative of curve function
% xpp,ypp : second derivative of curve function
% endpt1,endpt2: endpts for curve with boundaries / endvalues of parameter
% isclosed: true for closed curves

% TODO: think about endpoints / curves without bndries

% TODO: currently in 2D, generalise to 3D

%%

% distance squared from (x,y) to curve
d2 = @(t,x,y) (xs(t) - x).^2 + (ys(t) - y).^2;

% derivative of dist wrt t
g = @(t,x,y) 2*(xs(t) - x).*xp(t) + 2*(ys(t) - y).*yp(t);

% second derivative:
gp = @(t,x,y) 2*xp(t).*xp(t) + 2*(xs(t) - x).*xpp(t) + 2*yp(t).*yp(t) + 2*(ys(t) - y).*ypp(t);

ss = linspace(endpt1, endpt2, 10000);
dd = d2(ss, x, y);
[dontcare, i] = min(dd);
s_guess = ss(i);

%%
% Newton's method

tol = 1e-14;
endpttol = 0.1;
maxn = 255;

n = 1;
outsideCounter = 0;
s = s_guess;

while (true)
   n = n + 1;
   
   % Newton's method
   if (abs(gp(s,x,y)) > tol)
       snew = s - g(s,x,y) / gp(s,x,y);
   else
     % second deriv is zero
     snew = s;
   end
       
   if (~isclosed && (snew - endpt1 < -endpttol) || (snew - endpt2 > endpttol))
       outsideCounter = outsideCounter + 1;
   end
   
   if (outsideCounter >= 5)
       fail = 1;
       fprintf('*** probably diverging: (s_guess,n,s,snew)= %f, %d, %f, %f \n', s_guess,n,s,snew);
       break;
   elseif (n > maxn)
       fail = 1;
       fprintf('too many iterations: (n,snew,x,y) = %d, %f, %f, %f \n', n, snew, x, y);
       g(s,x,y)
       gp(s,x,y)
       break;
   elseif (abs(s-snew) < tol)
       if (isclosed || ((endpt1 <= snew) && (snew <= endpt2)))
            fail = 0;
            %fprintf('converged: (cpx,cpy,s,n)= %f, %f, %f, %d \n',xs(snew),ys(snew),snew,n);
       else
           fail = 1;
           fprintf('converged outside interval (cpx,cpy,s,n)= %f, %f, %f, %d \n',xs(snew),ys(snew),snew,n);
       end;
       break;       
   end
   
   % update
   s = snew;
end


if (fail) 
   % not converged: use endpt value
   
   dpt1 = d2(endpt1, xs(endpt1), ys(endpt1));
   dpt2 = d2(endpt2, xs(endpt2), ys(endpt2));
   [dpt,i] = min([dpt1 dpt2]);
   dist = sqrt(dpt);
   endpt = [endpt1, endpt2];
   cpx = xs(endpt(i)); cpy = ys(endpt(i));
   
   varargout = {1};
       
else   % converged
    
    cpx = xs(s); cpy = ys(s);
    dist = sqrt((cpx - x).^2 + (cpy - y).^2);
    
    varargout = {0};

end

