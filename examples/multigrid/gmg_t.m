function [v error_circ_inf_mg res] = gmg_t(M, DM, L, E, V, F, TMf2c, TMc2f, BAND, BDYG, n1, n2, start, w, uexact, Eplot, MAX)

% geometric multigrid method, if the solution still has the potential to be 
% 'improved', then invoke function 'helper_vcycle_test' to improve it

% allocate space for returning solution
v = zeros(size(V{start}));

% if the (relative) difference of v after two successive vcylces is less
% than 'tol', then stop doing vcycles
tolV = 1e-6;    
tolRes = 1e-10;    

% cnt and res are used for debugging and seeing effectiveness of multigrid
% method
% cnt is used to record total number of vcycles of multigrid method
cnt = 0;

error_circ_inf_mg = zeros(1,MAX);
res = zeros(1,MAX);

while cnt < MAX
   v1 = helper_vcycle_t(M, DM, L, E, V, F, TMf2c, TMc2f, BAND, BDYG, n1, n2, start, w);
   circplot_mg = Eplot{start}*v1;
   error = uexact - circplot_mg;
   error_circ_inf_mg(cnt+1) = norm(error,inf);
   if norm((v-v1))/norm(v) < tolV
      break;
   end
   v = v1;
   V{start} = v;
   res(cnt+1) = norm(Eplot{start}*(F{start} - M{start}*v),inf);
   r = res(cnt+1)
   if r < tolRes
       break;
   end
   
   cnt = cnt + 1;
end 

%number_of_vcycle = cnt

end
