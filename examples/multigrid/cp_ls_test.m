function [cpx, cpy, cpz, dist] = cp_ls_test(xg,yg,zg,cpf_tom,n)

%profile on

% Divide the points into several groups, and in each group, minize the sum
% of the squared distances. Hope this would be faster than do a for loop to 
% minimize for each grid point.

% 'n' is the number of points in each group. Can't be too large, other wise
% there are too many constraints, and too many components in the sum lead
% extra erros.
if nargin < 5
    n = 20;
end

cpx = xg;
cpy = yg;
cpz = zg;

len = length(xg);
m = floor(len/n);
res = mod(len,n); 

%[xg0, yg0, zg0] = cpf_tom(xg,yg,zg);
xg0 = xg;
yg0 = yg;
zg0 = zg;

% rearrange the the vectors a bit in order to use parfor
xg_mat = reshape(xg0(1:m*n),n,m);
yg_mat = reshape(yg0(1:m*n),n,m);
zg_mat = reshape(zg0(1:m*n),n,m);

cpx_mat = reshape(cpx(1:m*n),n,m);
cpy_mat = reshape(cpy(1:m*n),n,m);
cpz_mat = reshape(cpz(1:m*n),n,m);

tolc = 1e-16;
tolx = 1e-14;

parfor i = 1:m
    
    obj = @(X) dist2( X, xg_mat(:,i), yg_mat(:,i), zg_mat(:,i) );
    con = @(X) nonlinfcn(X);
    opts = optimset('GradObj','on',...%'Hessian','user-supplied',...
        'GradConstr','on','Algorithm','interior-point',...
        'TolCon',tolc,'TolX',tolx,'display','off');
    X0 = [xg_mat(:,i), yg_mat(:,i), zg_mat(:,i)];
    X0 = X0';
    X0 = X0(:);
    X = fmincon(obj,X0,[],[],[],[],[],[],con,opts);

    N = length(X);
    cpx_mat(:,i) = X(1:3:N-2);
    cpy_mat(:,i) = X(2:3:N-1);
    cpz_mat(:,i) = X(3:3:N);
end

cpx(1:m*n) = reshape(cpx_mat,m*n,1);
cpy(1:m*n) = reshape(cpy_mat,m*n,1);  
cpz(1:m*n) = reshape(cpz_mat,m*n,1);

if res >0 
    obj = @(X) dist2( X, xg(m*n+1 : len), yg(m*n+1 : len), zg(m*n+1 : len) );
    con = @(X) nonlinfcn(X);
    opts = optimset('GradObj','on',...%'Hessian','user-supplied',...
        'GradConstr','on','Algorithm','Interior-point',...
        'TolCon',tolc,'TolX',tolx,'display','off');
    X0 = [xg(m*n+1 : len), yg(m*n+1 : len), zg(m*n+1 : len)];
    X0 = X0';
    X0 = X0(:);
    X = fmincon(obj,X0,[],[],[],[],[],[],con,opts);

    N = length(X);
    cpx(m*n+1 : len) = X(1:3:N-2);
    cpy(m*n+1 : len) = X(2:3:N-1);
    cpz(m*n+1 : len) = X(3:3:N);
end

ind = xg == 0 & yg == 0 & zg == 0;
cpx(ind) = 1;
cpy(ind) = 0;
cpz(ind) = 0;

dist = (cpx-xg).^2 + (cpy-yg).^2 + (cpz-zg).^2;
dist = sqrt(dist);

%profile viewer

end

function [d2 grad_d2 hess_d2] = dist2(X,x0,y0,z0)

n = length(X);
x = X(1:3:n-2);
y = X(2:3:n-1);
z = X(3:3:n);

vec_d2 = (x-x0).^2 + (y-y0).^2 + (z-z0).^2;
d2 = sum(vec_d2);

grad_d2 = 2*[x-x0, y-y0, z-z0];
grad_d2 = grad_d2';
grad_d2 = grad_d2(:);

hess_d2 = 2*speye(3*length(x0));

end

function [c,ceq,grad_c,grad_ceq] = nonlinfcn(X)
% We only have equality constraint, non-equality constraints are empty.
c = [];
grad_c = [];

% level set: phi = (x-z^2)^2 + y^2 + z^2 - 1
% grad_phi = [2*(x-z^2)  2*y  4z^3-4*x*z+2*z]
n = length(X);
ceq = (X(1:3:n-2)-X(3:3:n).^2).^2 + X(2:3:n-1).^2 + X(3:3:n).^2 - 1;
% gradc_{i,j} = [\partial c_j / \partial x_i].
m = n/3;
gi = repmat((1:m)',1,3);
gj = reshape((1:n),3,m)'; 
gij = [2*(X(1:3:n-2)-X(3:3:n).^2) 2*X(2:3:n-1) 4*X(3:3:n).^3-4*X(1:3:n-2).*X(3:3:n)+2*X(3:3:n) ];
grad_ceq = sparse(gi,gj,gij,m,n);
grad_ceq = grad_ceq';
end
