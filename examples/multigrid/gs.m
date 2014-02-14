function [y] = gauss_seidel(M, x, f, cnt)

% n = length(b);
% 
% for j = 1:1:cnt
%    x(1) = (b(1) - A(1,2:n)*x(2:n,1))/A(1,1);
%    for i = 2:n-1
%        x(i) = (b(i) - A(i,1:i-1)*x(1:i-1,1) - A(i,i+1:n)*x(i+1:n,1))/A(i,i);
%    end
%    x(n) = (b(n) - A(n,1:n-1)*x(1:n-1,1))/A(n,n);
% end
% 

y = x;

for i = 1:1:cnt
    y = tril(M) \ (f - triu(M,1)*y);
end

end