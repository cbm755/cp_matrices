function [uexact] = uexactfn_mix_mode(th, n_mode)
uexact = 0;
for i = 1:1:n_mode
    uexact = uexact + sin(i*th)/i;
end
end