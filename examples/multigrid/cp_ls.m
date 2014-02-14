function [cpx, cpy, cpz, dist] = cp_ls(xg,yg,zg,cpf_tom)


dist = zeros(size(xg(:)));
len = length(xg(:));

%[cpx, cpy, cpz] = cpf_tom(xg,yg,zg);
cpx = xg;
cpy = yg;
cpz = zg;

parfor i = 1:len
    obj = @(X) dist2(X,xg(i),yg(i),zg(i));
    con = @(X) nonlinfcn(X);
    opts = optimset('GradObj','on','GradConstr','on','Algorithm','Interior-point',...
        'TolCon',1e-15,'TolFun',1e-14,'display','off');
    X0 = [cpx(i);cpz(i);cpy(i)];
    
    [X,f]= fmincon(obj,X0,[],[],[],[],[],[],con,opts);
    cpx(i) = X(1);
    cpy(i) = X(2);
    cpz(i) = X(3);
    dist(i) = sqrt(f);
end


end

function [d2 grad_d2 hess_d2] = dist2(X,x0,y0,z0)

d2 = (X(1)-x0)^2 + (X(2)-y0)^2 + (X(3)-z0)^2;
grad_d2 = 2*[X(1)-x0; X(2)-y0; X(3)-z0];

hess_d2 = 2*speye(3);

end

function [c,ceq,grad_c,grad_ceq] = nonlinfcn(X)
% We only have equality constraint, non-equality constraints are empty.
c = [];
grad_c = [];

% level set: phi = (x-z^2)^2 + y^2 + z^2 - 1
% grad_phi = [2*(x-z^2)  2*y  4z^3-4*x*z+2*z]
ceq = (X(1)-X(3)^2)^2 + X(2)^2 + X(3)^2 - 1;
% gradc_{i,j} = [\partial c_j / \partial x_i].
grad_ceq = [2*(X(1)-X(3)^2); 2*X(2); 4*X(3)^3-4*X(1)*X(3)+2*X(3) ];

end
