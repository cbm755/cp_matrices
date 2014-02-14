function [rhs] = rhsfn_mix_mode(th, n_mode)
rhs = 0;
for i = 1:1:n_mode
    rhs = rhs - i*sin(i*th);
end
end
