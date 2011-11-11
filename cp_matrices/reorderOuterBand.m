function [cp2, P, T]= reorderOuterBand(cp)
%REORDEROUTERBAND   Put the inner band before the outer band
%   The outerband contains the inner band and it sometimes might be
%   desirable to order the outerband so that the first n items are
%   the innerband.  Then innerband == outerband(1:n).
%
%   [cp2, P, T]= reorderOuterBand(cp) Supports the cpgrid structure.
%   'cp' should contain at least cp.iband and cp.oband.  'cp2' will be
%   a new cpgrid structure.  All things with some dimension with
%   length(oband) will be permuted automatically in cp2.  'P' and 'T'
%   are optional outputs for the permutation.  'P' is the permutation
%   matrix and 'T' is the same as a permutation index (i.e., v =
%   P*u reorders vector u.  v = u(T) does the same thing).

  if(1==0)
    cp.x1d = x1d;
    cp.y1d = y1d;
    cp.z1d = z1d;
    cp.xout = xout;
    cp.yout = yout;
    cp.zout = zout;
    cp.x = x;
    cp.y = y;
    cp.z = z;
    cp.cpx = cpx;
    cp.cpy = cpy;
    cp.cpz = cpz;
    cp.cpxout = cpxout;
    cp.cpyout = cpyout;
    cp.cpzout = cpzout;
    cp.L = L;
    cp.E = E;
    cp.R = R;
    cp.iband = iband;
    cp.oband = oband;
    cp.obandInParent = obandInParent;
    cp.ibandInParent = ibandInParent;  % R * obandInParent
  end

  iband = cp.iband;
  oband = cp.oband;

  T = zeros(size(oband));

  osz = length(oband);
  isz = length(iband);
  c = 0;
  for i=1:length(oband)
    I = find(oband(i) == iband);
    if (isempty(I))
      c = c + 1;
      j = isz + c;
    else
      j = I;
    end
    %disp( {i I j} )
    T(j) = i;
  end
  osz = length(oband);
  P = sparse(1:osz,T,ones(size(oband)), osz, osz, osz);

  % this will happen with the code below
  %cp2 = [];
  cp2.iband = iband;

  %cp2.oband = oband(T);

  names = fieldnames(cp);
  for i=1:length(names)
    fname = names{i};
    f = getfield(cp, fname);
    sz = size(f);
    if ( (sz(1) == osz) & (sz(2) == osz) )
      cp2 = setfield(cp2, fname, f(T,T));
    elseif (sz(1) == osz)
      cp2 = setfield(cp2, fname, f(T,:));
    elseif (sz(2) == osz)
      cp2 = setfield(cp2, fname, f(:,T));
    else
      cp2 = setfield(cp2, fname, f);
    end
  end

  %keyboard