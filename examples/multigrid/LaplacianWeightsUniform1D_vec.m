function [ww] = LaplacianWeightsUniform1D_vec(dx, a)

ww = zeros(length(a),4);
ww(:,1) = dx-a;
ww(:,2) = 3*a-2*dx;
ww(:,3) = dx-3*a;
ww(:,4) = a;

ww = ww / dx^3;

end