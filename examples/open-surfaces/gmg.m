function [v, error_inf, res] = gmg(M, L, E, E_out_out, E_out_in, V, F, TMf2c, TMc2f, BAND, BDYG, n1, n2, start, w, Eplot, uexact, MAX)

% geometric multigrid method, if the solution still has the potential to be 
% 'improved', then invoke function 'helper_vcycle_test' to improve it

% allocate space for returning solution
v = zeros(size(V{start}));

% if the (relative) difference of v after two successive vcylces is less
% than 'tol', then stop doing vcycles
tolV = 1e-8;    
tolRes = 1e-6;

% cnt and res are used for debugging and seeing effectiveness of multigrid
% method
% cnt is used to record total number of vcycles of multigrid method
cnt = 1;

error_inf = zeros(1,MAX);
res = zeros(1,MAX);
error_inf(1) = 1;
res(1) = 1;

not_bdy = ~logical(BDYG{start});

while cnt < MAX
   v1 = helper_vcycle(M, L, E, E_out_out, E_out_in, V, F, TMf2c, TMc2f, BAND, BDYG, n1, n2, start, w);
   error_inf(cnt+1) = norm(Eplot{start}*v1-uexact,inf) / norm(uexact,inf);
   res(cnt+1) = norm((F{start}(not_bdy) - M{start}(not_bdy,:)*v1),inf) / norm(F{start}(not_bdy),inf);
   r = res(cnt+1)
   if r < tolRes
       break;
   end
   if norm(v-v1,inf)/norm(v,inf) < tolV
      break;
   end
   v = v1;
   V{start} = v;
   
   
   cnt = cnt + 1;
end 

%number_of_vcycle = cnt

end
