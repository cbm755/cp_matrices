function [v error_circ_inf_mg, res] = gmg_2b(M, L, E, Ei, D, V, F, TMf2c, TMc2f, BAND, BDYG, n1, n2, start, w, uexact, Eplot, MAX)

% geometric multigrid method, if the solution still has the potential to be 
% 'improved', then invoke function 'helper_vcycle_test' to improve it

% allocate space for returning solution
v = zeros(length(V{start}),1);

% if the (relative) difference of v after two successive vcylces is less
% than 'tol', then stop doing vcycles
tol = 1e-10;    

% cnt and res are used for debugging and seeing effectiveness of multigrid
% method
% cnt is used to record total number of vcycles of multigrid method
cnt = 0;
% res is used to measure the norm of residual after each vcycle

error_circ_inf_mg = zeros(1,MAX);
res = zeros(1,MAX);

while cnt < MAX
   v1 = helper_vcycle_2b(M, L, E, Ei, D, V, F, TMf2c, TMc2f, BAND, BDYG, n1, n2, start, w);
   circplot_mg = Eplot{start}*v1;
   error_circ_inf_mg(cnt+1) = max(abs( uexact - circplot_mg ));
   if norm(v-v1)/norm(v) < tol
      break;
   end
   v = v1;
   V{start} = v;
%   D = diag(diag(L{1}));
%   M_bar = E{1}*(L{1}-D) + D;
   res(cnt+1) = norm(Eplot{start}*(F{start} - M{start}*v),inf);
   res(cnt+1)
   cnt = cnt +1;
end 

%number_of_vcycle = cnt

end
