function [v, error_circ_inf_mg, res] = gmg_test(L, a_Ebar, a_Edouble, a_Etriple, E_out_out, E_out_in, Ecp_Omega_S, Ecp_f2c_Omega, Ecp_f2c_S, Ecp_f2c_Omega_S, V, F, FonS, TMf2c, TMc2f, BAND, BDYG, n1, n2, start, w, uexact, MAX)

% geometric multigrid method, if the solution still has the potential to be 
% 'improved', then invoke function 'helper_vcycle_test' to improve it

% allocate space for returning solution
v = zeros(size(V{start}));

% if the (relative) difference of v after two successive vcylces is less
% than 'tol', then stop doing vcycles
tolV = 1e-6;    
tolRes = 1e-7;

% cnt and res are used for debugging and seeing effectiveness of multigrid
% method
% cnt is used to record total number of vcycles of multigrid method
cnt = 1;

error_circ_inf_mg = zeros(1,MAX);
res = zeros(1,MAX);
error_circ_inf_mg(1) = 1;
res(1) = 1;

not_bdy = ~BDYG{start};

while cnt < MAX
   v1 = helper_vcycle_test(L, a_Ebar, a_Edouble, a_Etriple, E_out_out, E_out_in, Ecp_Omega_S, Ecp_f2c_Omega, Ecp_f2c_S, Ecp_f2c_Omega_S, V, F, FonS, TMf2c, TMc2f, BAND, BDYG, n1, n2, start, w);
   error = uexact{start} - v1;
   error_circ_inf_mg(cnt+1) = norm(error(not_bdy),inf) / norm(uexact{start}(not_bdy),inf);
   res(cnt+1) = norm((F{start}(not_bdy) - L{start}(not_bdy,:)*v1),inf) / norm(F{start}(not_bdy),inf);
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

%number_of_vcycle = cnt

end
