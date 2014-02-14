function [v, error_circ_inf_mg, res, n_vc, err_normal_v, err_normal_f] = gmg_debug1(M, L, E, V, F, TMf2c, TMc2f, BAND, BDYG, n1, n2, start, w, uexact, Eplot, MAX)

% geometric multigrid method, if the solution still has the potential to be 
% 'improved', then invoke function 'helper_vcycle_test' to improve it

% allocate space for returning solution
v = zeros(size(V{start}));

% if the (relative) difference of v after two successive vcylces is less
% than 'tol', then stop doing vcycles
tolV = 1e-8;    
tolRes = 1e-14;

% cnt and res are used for debugging and seeing effectiveness of multigrid
% method
% cnt is used to record total number of vcycles of multigrid method
cnt = 1;

error_circ_inf_mg = zeros(1,MAX);
res = zeros(1,MAX);
error_circ_inf_mg(1) = 1;
res(1) = 1;

n_level = length(BAND);
N = 2*(n_level-start);
err_normal_v = zeros(N,MAX);
err_normal_f = zeros(N,MAX);

while cnt < MAX
   [v1,err_normal_v(:,cnt),err_normal_f(:,cnt)] = helper_vcycle_debug(M, L, E, V, F, TMf2c, TMc2f, BAND, BDYG, n1, n2, start, w);
   circplot_mg = Eplot{start}*v1;
   error = uexact - circplot_mg;
   error_circ_inf_mg(cnt+1) = norm(error,inf) / norm(uexact,inf);
   %res(cnt+1) = norm(Eplot{start}*(F{start} - M{start}*v),inf);
   res(cnt+1) = norm((F{start} - M{start}*v1),inf) / norm(F{start},inf);
   r = res(cnt+1)
   if r < tolRes
       break;
   end
   if norm((v-v1))/norm(v) < tolV
      break;
   end
   v = v1;
   V{start} = v;
   
   
   cnt = cnt + 1;
end 

n_vc = cnt

end
