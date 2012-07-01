function [cpx cpy cpz dist] = cpSurfOfRevolution(x, y, z, cpf, whichaxis, cpfdata)
%CPSURFOFREVOLUTION  Rotate a curve to make a surface
%   [cpx cpy cpz dist] = cpSurfOfRevolution(x, y, z, @cpEllipse)
%      Takes a curve in 2D (here an ellipse) and spins it around the
%      x-axis to create an ellipsoid.  In general, pass a function
%      handle as the third argument.
%
%   [cpx cpy cpz dist] = cpSurfOfRevolution(x, y, z, cpf, 'y')
%      Rotates around 'y' axis instead.
%
%   [cpx cpy cpz dist] = cpSurfOfRevolution(x, y, z, cpf, 'x', P)
%      Passes the cell area P as the extra arguments of the
%      function handle cpf.
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

  switch whichaxis
    case 'x'
      [th,r,zz] = cart2pol(y,z,x);
      if (isempty(cpfdata))
        [cpzz cpr dist] = cpf(zz, r);
      else
        [cpzz cpr dist] = cpf(zz, r, cpfdata{:});
      end
      [cpy,cpz,cpx] = pol2cart(th, cpr, cpzz);

    case 'y'
      [th,r,zz] = cart2pol(z,x,y);

      if (isempty(cpfdata))
        [cpr cpzz dist] = cpf(r, zz);
      else
        [cpr cpzz dist] = cpf(r, zz, cpfdata{:});
      end

      [cpz,cpx,cpy] = pol2cart(th, cpr, cpzz);
    otherwise
      error('axis not implemented');
  end
