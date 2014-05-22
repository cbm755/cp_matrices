function [ww, dd] = LaplacianWeightsNonUniform1D_vec(dx, a)

flag1 = a<dx/2;
flag2 = a>=dx/2;

ww = zeros(length(a),4);
dd = zeros(size(a));

ww(flag1,1) = (6*dx-4*a(flag1))./( 6*dx^2*(dx+a(flag1)) );
ww(flag2,1) = 4*(dx-a(flag2))./( 3*dx^2*(dx+a(flag2)) );

ww(flag2,2) = (4*a(flag2)-2*dx)./(2*dx^2*a(flag2));

ww(flag1,3) = (2*dx-4*a(flag1))./( 2*dx^2*(dx-a(flag1)) );

ww(flag1,4) = 4*a(flag1)./( 3*dx^2*(2*dx-a(flag1)) );
ww(flag2,4) = (2*dx+4*a(flag2))./( 6*dx^2*(2*dx-a(flag2)) );

dd(flag1) = (4*dx-6*a(flag1))./( (dx+a(flag1)).*(dx-a(flag1)).*(2*dx-a(flag1)) );
dd(flag2) = (6*a(flag2)-2*dx)./( (dx+a(flag2)).*a(flag2).*(2*dx-a(flag2)) );
 
end