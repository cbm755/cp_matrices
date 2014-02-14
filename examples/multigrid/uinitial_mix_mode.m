function [u_init] = uinitial_mix_mode(th, n_mode, random_coefficient)
u_init = 0;

for i = 1:1:n_mode
    u_init = u_init + random_coefficient(i)*sin(i*th);
end
end