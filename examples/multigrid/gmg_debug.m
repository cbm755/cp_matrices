function [v_cell, error_circ_inf_mg, res1, res2, err_matlab, number_of_vcycle] = ...
    gmg_debug(M, L, E, V, F, TMf2c, TMc2f, BAND, BDYG, n1, n2, start, w, uexact, u_matlab, Eplot, ind_pt, uexact_pt, MAX)

% geometric multigrid method, if the solution still has the potential to be 
% 'improved', then invoke function 'helper_vcycle_test' to improve it

% allocate space for returning solution
v_cell = cell(MAX,1);
v = zeros(size(V{start}));

% if the (relative) difference of v after two successive vcylces is less
% than 'tol', then stop doing vcycles
tol = 1e-6;    

% cnt and res are used for debugging and seeing effectiveness of multigrid
% method
% cnt is used to record total number of vcycles of multigrid method
cnt = 1;

dx = 0.003125*2^(start-1);
error_circ_inf_mg = zeros(1,MAX);
res1 = zeros(1,MAX);
res2 = zeros(1,MAX);
err_matlab = zeros(1,MAX);
u_mat = u_matlab{start};

while cnt < MAX
   if isa(Eplot,'cell')
       res1(cnt) = norm(Eplot{start}*(F{start} - M{start}*V{start}),inf);
       %res(cnt) = norm((F{start} - M{start}*V{start}))*dx;
       res2(cnt) = norm((F{start} - M{start}*V{start}),inf);
   else
       res1(cnt) = norm(Eplot*(F{start} - M{start}*V{start}),inf);
       res2(cnt) = norm((F{start} - M{start}*V{start}),inf);
   end
   r = res1(cnt)

   if isa(Eplot,'cell')
       circplot_mg = Eplot{start}*V{start};
   else
       circplot_mg = Eplot*V{start};
   end
   error = circplot_mg - uexact;
   
   %error = error - (V{start}(ind_pt{start}) - uexact_pt );
   error_circ_inf_mg(cnt) = norm(error,inf) / norm(uexact,inf);
   
   
   %err_matlab(cnt) = norm(V{start} - u_mat)*dx;   
   err_matlab(cnt) = norm((V{start} - u_mat), inf);
   
   %v1 = helper_vcycle_M(M, L, E, V, F, TMf2c, TMc2f, BAND, BDYG, n1, n2, start, w);
   v1 = helper_vcycle(M, L, E, V, F, TMf2c, TMc2f, BAND, BDYG, n1, n2, start, w);
   %v1 = helper_vcycle_test(M, L, E, V, F, TMf2c, TMc2f, BAND, BDYG, n1, n2, start, w);
   
   %v1 = gauss_seidel(M{start}, V{start}, F{start});
   %v1 = weighted_jacobi(L{start},V{start},F{start},1);
   %v1 = E{start}*v1;
   
%    figure
%    x = linspace(0, 2*pi, size(uexact,1));
%    y = Eplot{start}*v1;
%    plot(x,y,'b-', x,uexact,'r');
%    xlim([0 2*pi])

%    figure
%    x = linspace(0, 2*pi, size(uexact,1));
%    plot(x,error);
%    xlim([0 2*pi])
%    

   if norm((v-v1))/norm(v) < tol
      break;
   end
%   if cnt>=2 && res1(cnt)/res1(cnt-1) > 0.8
%     break;
%   end
   v = v1;
   V{start} = v1;
   v_cell{cnt} = v1;   
   
   cnt = cnt + 1;
end 

number_of_vcycle = cnt-1;

end
