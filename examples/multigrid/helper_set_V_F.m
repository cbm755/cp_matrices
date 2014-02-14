function [V,F] = helper_set_V_F(BAND, uold)
V = cell(size(BAND));
F = cell(size(BAND));
V{1} = uold;
F{1} = uold;
for i = 2:1:length(BAND)
    V{i} = zeros(length(BAND{i}),1);
    %F{i} = zeros(length(BAND{i}),1);
end

end
