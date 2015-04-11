function [f] = laplace_beltrami_ls2d( phi, u, has_t )

% Given an expression of the level set representation for a surface in R^3: 'phi',
% and the expression of a  function 'u' w.r.t to cartesian coordinates, compute the
% expression of 'f = \Delta_S u'.

if nargin < 3
    has_t = false
end

syms st sx sy;

sphi = phi( sx, sy );
grad_sphi = jacobian( sphi, [sx sy] );
snormal = grad_sphi / sqrt( grad_sphi * transpose(grad_sphi) );

if has_t
    su = u(st, sx, sy);
else
    su = u( sx, sy );
end
grad_su = jacobian( su, [sx sy] );
gradS_su = simplify( grad_su - ( snormal*transpose(grad_su) )*snormal );


sf = divergence(gradS_su,[sx sy]) - ...
    ( snormal*gradient(gradS_su(1),[sx sy]) )*snormal(1) - ...
    ( snormal*gradient(gradS_su(2),[sx sy]) )*snormal(2);

sf = simplify(sf);

if has_t
    f = matlabFunction(sf, 'vars', [st sx sy]);
else
    f = matlabFunction(sf, 'vars', [sx sy]);
end

end