 function [f] = laplace_beltrami_ls3d( phi, u )

% Given an expression of the level set representation for a surface in R^3: 'phi',
% and the expression of a  function 'u' w.r.t to cartesian coordinates, compute the
% expression of 'f = \Delta_S u'.

syms sx sy sz;

sphi = phi( sx, sy, sz );
grad_sphi = jacobian( sphi, [sx sy sz] );
snormal = grad_sphi / sqrt( grad_sphi * transpose(grad_sphi) );

su = u( sx, sy, sz );
grad_su = jacobian( su, [sx sy sz] );
gradS_su = simplify( grad_su - ( snormal*transpose(grad_su) )*snormal );


sf = divergence(gradS_su,[sx sy sz]) - ...
    ( snormal*gradient(gradS_su(1),[sx sy sz]) )*snormal(1) - ...
    ( snormal*gradient(gradS_su(2),[sx sy sz]) )*snormal(2) - ...
    ( snormal*gradient(gradS_su(3),[sx sy sz]) )*snormal(3);

sf = simplify(sf);

f = matlabFunction(sf, 'vars', [sx sy sz]);

end