function [pass, str] = test_vectorized_vs_loop_lagrange()
  str = 'test the vectorized vs scalar lagrange weights';

  pass = [];
  for N=1:8
    pass = [pass helper(N)];
  end
end

function pass = helper(N)
%N = 4;   % test N = degree + 1 interp
m = 5000;
n = m/10;  % this many are off-grid (non-pathelogical)

dx = rand(m,1);
xg = 2*rand(m,1)-1;
x = xg + floor(N*rand(m,1)).*dx;

% the 1.2 makes some of them extrapolations
x(1:n) = xg(1:n) + N*( 1.2*(rand(n,1))-0.1 ) .* dx(1:n);



T1 = cputime();
w1 = LagrangeWeights1D_vec(xg,x,dx,N);
T1 = cputime()-T1;
%fprintf('  vector code      (N=%d) elapsed time=%g seconds\n', N, T);


T2 = cputime();
w2 = zeros(size(w1));
for i=1:m
  w2(i,:) = LagrangeWeights1D(xg(i), x(i), dx(i), N);
end
T2 = cputime()-T2;
%fprintf('  loop code        (N=%d) elapsed time=%g seconds\n', N, T);


T3 = cputime();
w3 = zeros(size(w1));
for i=1:m
  w3(i,:) = LagrangeWeights1D_vec(xg(i), x(i), dx(i), N);
end
T3 = cputime()-T3;
fprintf('  N=%d: [vec,loop,loop-w-vec-code] elapsed times=[%5.3g,%5.3g,%5.3g] seconds\n', N, T1,T2,T3);

pass = (max(max(abs(w1-w2))) == 0) & ...
       (max(max(abs(w1-w3))) == 0);

end