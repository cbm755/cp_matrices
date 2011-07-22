function [cpx,cpy,cpz, dist, bdy1] = cpFishBowl(x,y,z, R, xc,yc,zc)
%CPHEMISPHERE  Closest point function for a hemisphere.
%   The hemisphere consists of those points with z >= 0.
%   [cpx,cpy,cpz, dist, bdy1] = cpHemisphere(x,y,z)
%      A unit hemisphere with "center" (sphere center) at the origin.
%      "bdy1" is non-zero for points on the boundary.
%   [cpx,cpy,cpz, dist, bdy1] = cpHemisphere(x,y,z, R)
%      A radius R hemisphere centered at the origin.
%   [cpx,cpy,cpz, dist, bdy1] = cpHemisphere(x,y,z, R, xc,yc,zc)
%      A radius R hemisphere centered at (xc,yc,zc).
%
%   Code is vectorized: any size/shape for x should work.


  % default radius of 1
  if (nargin < 4)
    R = 1;
  end
  if (nargin == 5) | (nargin == 6)
    error('must specify all of (xc,yc,zc)');
  end
  % default center is the origin
  if (nargin < 7)
    xc = 0; yc = 0; zc = 0;
  end

  % shift to the origin
  x = x - xc;
  y = y - yc;
  z = z - zc;

  [cpx, cpy, cpz] = cpSphere(x, y, z, R);

  % should not include the bdy1 itself (hence not z <= 0)
  
 
   
  dist = sqrt( (x-cpx).^2 + (y-cpy).^2 + (z-cpz).^2 );
 

 
  
  % for points above zt = 0.4
   
  
  zt = 0.4;
  
  [th1, phi, rr] = cart2sph(x, y, z);
  bdy2 = (sin(phi) > zt) & (z>0);
  cpx0 = R*cos(zt)*cos(th1);
  cpy0 = R*cos(zt)*sin(th1);
   
  ind = find(bdy2==1);
  
  
  n_ind = length(ind);
  
  for ii = 1:n_ind  %n_ind
  
    temp_dist = sqrt( (x(ind(ii))-cpx0(ind(ii))).^2 + (y(ind(ii))-cpy0(ind(ii))).^2 + (z(ind(ii))-zt).^2 );
      dist(ind(ii)) = temp_dist; 
      cpx(ind(ii)) = cpx0(ind(ii));
      cpy(ind(ii)) = cpy0(ind(ii));
      cpz(ind(ii)) = zt;
      
  end
  

  
  
  bdy1 = bdy2;
  
  
 
  
  
  % shift back to center
  cpx = cpx + xc;
  cpy = cpy + yc;
  cpz = cpz + zc;
