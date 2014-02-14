function [cpx, cpy, cpz, f] = cp_ls_test_bak(xg,yg,zg)

obj = @(X) dist2(X,xg,yg,zg);
con = @(X) nonlinfcn(X);
opts = optimset('GradObj','on',...%'Hessian','user-supplied',...
    'GradConstr','on','Algorithm','interior-point','display','iter');
X0 = [xg, yg, zg];
X0 = X0';
X0 = X0(:);
[X,f]= fmincon(obj,X0,[],[],[],[],[],[],con,opts);

n = length(X);
cpx = X(1:1:n/3);
cpy = X(n/3+1:1:2*n/3);
cpz = X(2*n/3+1:1:n);

end

function [d2 grad_d2 hess_d2] = dist2(X,x0,y0,z0)

n = length(X);
m = n/3;
vec_d2 = (X(1:1:m)-x0).^2 + (X(m+1:1:2*m)-y0).^2 + (X(2*m+1:1:n)-z0).^2;
d2 = sum(vec_d2);

grad_d2 = 2*[X(1:1:m)-x0; X(m+1:1:2*m)-y0; X(2*m+1:1:n)-z0];

hess_d2 = 2*speye(3*length(x0));

end

function [c,ceq,grad_c,grad_ceq] = nonlinfcn(X)
% We only have equality constraint, non-equality constraints are empty.
c = [];
grad_c = [];

% level set: phi = (x-z^2)^2 + y^2 + z^2 - 1
% grad_phi = [2*(x-z^2)  2*y  4z^3-4*x*z+2*z]
n = length(X);
m = n/3;
ceq = (X(1:m)-X(2*m+1:n).^2).^2 + X(m+1:2*m).^2 + X(2*m+1:n).^2 - 1;
% gradc_{i,j} = [\partial c_j / \partial x_i].

gi = repmat((1:m)',1,3);
gj = reshape((1:n),3,m)'; 
gij = [2*(X(1:m)-X(2*m+1:n).^2) 2*X(m+1:2*m) 4*X(2*m+1:n).^2-4*X(1:m).*X(2*m+1:n)+2*X(2*m+1:n) ];
grad_ceq = sparse(gi,gj,gij,m,n);
grad_ceq = grad_ceq';
end
