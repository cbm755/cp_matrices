function t1 = fd_weno5_1d(v, WENOEPS)
%FD_WENO5_1D  Fifth-order WENO finite difference in 1D
%   Input is a matrix v with 5 columns and we perform WENO on each
%   row.
%
%   Example: if f is the positive flux than we reconstruct
%   f_{j+1/2} as:
%
%     v = [f_{j-2}, f_{j-1}, f_{j}, f_{j+1}, f_{j+2}]
%     f_jph = fd_weno5_1d( v )
%
%   For the negative flux g, we can reconstruct g_{j+1/2} using the
%   same routine with the inputs in the opposite order:
%
%     v = [g_{j+3}, g_{j+2}, g_{j+1}, g_{j}, g_{j-1}]
%     g_jph = fd_weno5_1d( v )
%
%   Can also use fd_weno5_1d(v, 1e-10) to use a different value of
%   the WENO epsilon parameter (defaults to 1e-6).
%
%   Follows Osher-Fedkiw, pg 33--37

  if nargin < 2
    WENOEPS = 1e-6;
  end

  % TODO: could let user for the ENO choices too
  ForceStencil = -10;

  if (ForceStencil == -10)
    % smoothness indicators: S0, S1, S2
    w0 = 13.0/12.0 * ( v(:,1+0) -2.0*v(:,1+1) + v(:,1+2) ).^2  ...
         +    0.25 * ( v(:,1+0) -4.0*v(:,1+1) + 3.0*v(:,1+2) ).^2;
    w1 = 13.0/12.0 * ( v(:,1+1) -2.0*v(:,1+2) + v(:,1+3) ).^2  ...
         +    0.25 * ( v(:,1+1) - v(:,1+3) ).^2;
    w2 = 13.0/12.0 * ( v(:,1+2) -2.0*v(:,1+3) + v(:,1+4) ).^2  ...
         + 0.25 * ( 3.0*v(:,1+2) -4.0*v(:,1+3) +  v(:,1+4) ).^2;

    % alpha values
    % TODO: O&F specifies a different WENOEPS, pg 36
    w0 = 0.1 ./ ( WENOEPS + w0 ).^2;
    w1 = 0.6 ./ ( WENOEPS + w1 ).^2;
    w2 = 0.3 ./ ( WENOEPS + w2 ).^2;

    % the weights
    t1 = w0 + w1 + w2;
    w0 = w0 ./ t1;
    w1 = w1 ./ t1;
    w2 = w2 ./ t1;
  elseif (ForceStencil == 10)  % force the smooth region values
    w0 = 0.1;
    w1 = 0.6;
    w2 = 0.3;
  end

  t1 = w0 .* ( 2.0*v(:,1) - 7.0*v(:,2) + 11.0*v(:,3) ) / 6.0  ...
     + w1 .* (    -v(:,2) + 5.0*v(:,3) +  2.0*v(:,4) ) / 6.0  ...
     + w2 .* ( 2.0*v(:,3) + 5.0*v(:,4) -      v(:,5) ) / 6.0;

