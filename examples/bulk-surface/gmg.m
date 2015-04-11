function [u, v, error_inf_u, error_inf_v] = ...
    gmg(Au, Lu, Eu, Mu_coarsest, Lv, Eoo_v, Eoi_v, Bvu, Buv, U, V, G, F, TMf2c, TMc2f, TMf2c_S, TMc2f_S, ...
            a_band, a_band_S, a_bdyg, n1, n2, start, w, Eplot, uexact, vexact, MAX)

% geometric multigrid method, if the solution still has the potential to be 
% 'improved', then invoke function 'helper_vcycle_test' to improve it

% allocate space for returning solution
u = zeros(size(G{start}));
v = zeros(size(V{start}));

% if the (relative) difference of v after two successive vcylces is less
% than 'tol', then stop doing vcycles
tolV = 1e-6;    

% cnt and res are used for debugging and seeing effectiveness of multigrid
% method
% cnt is used to record total number of vcycles of multigrid method
cnt = 1;

error_inf_u = zeros(1,MAX);
error_inf_v = zeros(1,MAX);
res_u = zeros(1,MAX);
res_v = zeros(1,MAX);
error_inf_u(1) = 1;
error_inf_v(1) = 1;
res_u(1) = 1;
res_v(1) = 1;

not_bdy = ~a_bdyg{start};

while cnt < MAX
   [u1,v1] = helper_vcycle(Au, Lu, Eu, Mu_coarsest, Lv, Eoo_v, Eoi_v, Bvu, Buv, U, V, G, F, ...
                TMf2c, TMc2f, TMf2c_S, TMc2f_S, a_band, a_band_S, a_bdyg, n1, n2, start, w);
   error_u = uexact{start} - Eplot{start}*u1;
   error_v = vexact{start} - v1;
   error_inf_u(cnt+1) = norm(error_u,inf) / norm(uexact{start},inf);
   error_inf_v(cnt+1) = norm(error_v(not_bdy),inf) / norm(vexact{start}(not_bdy),inf);
   res_u(cnt+1) = norm((G{start} - Lu{start}*u1),inf) / norm(G{start},inf);
   res_v(cnt+1) = norm((F{start}(not_bdy) - Lv{start}(not_bdy,:)*v1),inf) / norm(F{start}(not_bdy),inf);
   r = res_v(cnt+1)

   if norm((v(not_bdy)-v1(not_bdy)),inf)/norm(v(not_bdy),inf) < tolV && norm(u-u1,inf)/norm(u,inf) < tolV
      break;
   end
   v = v1;
   u = u1;
   V{start} = v;
   U{start} = u;
   
   
   cnt = cnt + 1;
end 

%number_of_vcycle = cnt

end
