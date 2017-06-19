function [cpx, cpy, cpz, dist, varargout] = cpSurfOfRevolution(x, y, z, cpf, whichaxis, cpfdata)
%CPSURFOFREVOLUTION  Rotate a curve to make a surface
%   [cpx cpy cpz dist] = cpSurfOfRevolution(x, y, z, @cpEllipse)
%      Takes a curve in 2D (here an ellipse) and spins it around the
%      x-axis to create an ellipsoid.  In general, pass a function
%      handle as the fourth argument.
%
%   [cpx cpy cpz dist] = cpSurfOfRevolution(x, y, z, cpf, 'y')
%      Rotates around 'y' axis instead.
%
%   [cpx cpy cpz dist] = cpSurfOfRevolution(x, y, z, cpf, 'x', P)
%      Passes the cell array P as the extra arguments of the
%      function handle cpf.
%
%   [cpx cpy cpz dist bdy] = cpSurfOfRevolution(...)
%      Any extra outputs (such as bdy) are requested from cpf.
%
%   One needs to be a little careful with this.  The curve should
%   either be symmetric about the axis (or defined in the upper half
%   plane only).  This code converts the y-z plane into (r,th) polar
%   coors and passes those to cpf.  In particular r >= 0.


  % defaults
  if (nargin < 4)
    cpf = @cpEllipse;
  end
  if (nargin < 5)
    whichaxis = 'x';
  end
  if (nargin < 6)
    cpfdata = {};
  end

  % 5th output is (probably) bdy.  Anyway, we just pass on whatever
  % we get from the 2D cp-function.
  if (nargout >= 5)
    P = cell(1,nargout - 4);
  else
    P = {};
  end

  switch whichaxis
    case 'x'
      [th,r,zz] = cart2pol(y,z,x);
      [cpzz, cpr, dist, P{:}] = cpf(zz, r, cpfdata{:});
      [cpy,cpz,cpx] = pol2cart(th, cpr, cpzz);
    case 'y'
      [th,r,zz] = cart2pol(z,x,y);
      [cpr, cpzz, dist, P{:}] = cpf(r, zz, cpfdata{:});
      [cpz,cpx,cpy] = pol2cart(th, cpr, cpzz);
    otherwise
      error('invalid axis input: can only rotate around x or y axis');
  end

  varargout = P;

