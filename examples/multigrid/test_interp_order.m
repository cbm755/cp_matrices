%% test impact of polynomial interpolation order to solve poisson equation on a circle

%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../../cp_matrices');

% add functions for finding the closest points
addpath('../../surfaces');


n = 6;
err = zeros(6,n);
for p = 2:6
   err(p,:) = test_interp_order_backslash(p);
end

x = 0.1./2.^(n-1:-1:0);
N = 1./x;

a = 5e-2./4.^(n-2:-1:1);
loglog(N,err(2,:),'r.-',N,err(3,:),'b^--',N,err(4,:),'g-s',N,err(5,:),'m-x',N,err(6,:),'k+-',N(3:n-1),a(1:n-3));
legend('p=2','p=3','p=4','p=5','p=6')

%loglog(N,error_inf_matlab,'^--',N(3:n-1),a(1:n-3));
s = sprintf('%.0f-nd order',2);
text(N(4),1e-3,s)

set(gca,'Fontsize',20)
xlabel('\fontsize{20} N=1/dx')
ylabel('\fontsize{20} ||error||_{\infty}')
%title('poisson equation on a circle: sin(\theta)+sin(11\theta)')
title('\fontsize{20} poisson equation on a circle: sin(\theta)+sin(11\theta)')
