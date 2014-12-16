function [v] = fmg(M, L, E, V, F, TMf2c, TMc2f, BAND, BOUNDARY, n1, n2, nc, coarsest_level, w)

n_level = length(BAND);
for i = 1:1:n_level-1
    F{i+1} = TMf2c{i}*F{i};
end

V{n_level} = M{n_level} \ F{n_level};

for start = n_level-1:-1:coarsest_level
   V{start} = TMc2f{start}*V{start+1};
   for cnt = 1:1:nc
       for i = start+1:1:n_level
           V{i} = zeros(size(BAND{i}));
       end
       V{start} = helper_vcycle(M, L, E, V, F, TMf2c, TMc2f, BAND, BOUNDARY, n1, n2, start, w);
   end
   res = norm( F{start} - M{start}*V{start}, inf );
end

v = V{coarsest_level};

end
