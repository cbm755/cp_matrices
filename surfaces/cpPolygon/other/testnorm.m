function n=testnorm(x,y)
n = 0*x;

W = [1 -1/sqrt(2); 0 1/sqrt(2)];
W = [0 -1/sqrt(2); 1 -1/sqrt(2)];
Winv = W^(-1);

a = Winv(1,1)*x(:) + Winv(1,2)*y(:);
b = Winv(2,1)*x(:) + Winv(2,2)*y(:);

n(:) = max( abs(a), abs(b) );
%n(:) = abs(a)+abs(b);