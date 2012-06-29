function [cpx cpy cpz dist] = cpSurfOfRevolution(x, y, z, cpf, whichaxis, cpfdata)
%CPSURFOFREVOLUTION  Rotate a curve to make a surface
%
%   TODO: whaxis, cen, cellarray of options for cpf


  % defaults
  if (nargin < 4)
    cpf = @cpCircle;
  end
  if (nargin < 5)
    whichaxis = 'z';
  end
  if (nargin < 6)
    cpfdata = {};
  end

  switch whichaxis
    case 'x'
      [th,r,zz] = cart2pol(y,z,x);
    case 'y'
      [th,r,zz] = cart2pol(z,x,y);
    case 'z'
      [th,r,zz] = cart2pol(x,y,z);
    otherwise
      error('axis not implemented');
  end

  if (isempty(cpfdata))
    [cpzz cpr dist] = cpf(zz, r);
  else
    [cpzz cpr dist] = cpf(zz, r, cpfdata{:});
  end

  switch whichaxis
    case 'x'
      [cpy,cpz,cpx] = pol2cart(th, cpr, cpzz);
    case 'y'
      [cpz,cpx,cpy] = pol2cart(th, cpr, cpzz);
    case 'z'
      [cpx,cpy,cpz] = pol2cart(th, cpr, cpzz);
    otherwise
      error('axis not implemented');
  end


