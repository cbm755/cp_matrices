N = 2000;
dx = 2*pi/N;
th = (0:dx:2*pi)';
f = - exp( cos(th) ).*( sin(th).^2 - cos(th) );
f = f(2:N);
uexactfn = @(th) exp( cos(th) );

n = N-1;
e = ones(n,1);
M = spdiags([-e 2*e -e], -1:1, n, n);

f = f*dx.^2;
f(1) = f(1) + exp(1);
f(N-1) = f(N-1) + exp(1);


u = M \ f;

norm(u - uexactfn(th(2:N,1)), inf)
