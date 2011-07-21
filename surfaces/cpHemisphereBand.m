function [cpx,cpy,cpz, dist, bdy1] = cpHemisphereBand(x,y,z, R, xc,yc,zc)
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
  
  if 0
  
  bdy1 = (z < 0);

  
  
  
  % points with z < 0 map to the z=0 plane circle of radius R
  [th, r, zp] = cart2pol(x, y, z);
  cpth = th;
  cpr = R*ones(size(th));
  cpzp = zeros(size(th));
  [cpx2, cpy2, cpz2] = pol2cart(cpth, cpr, cpzp);

  cpx(bdy1) = cpx2(bdy1);
  cpy(bdy1) = cpy2(bdy1);
  cpz(bdy1) = cpz2(bdy1);
 
   
  end
   
  dist = sqrt( (x-cpx).^2 + (y-cpy).^2 + (z-cpz).^2 );
 

  
 % for points below zt = -0.25
   
  zt = 0.25 ;
  [th1, phi, rr] = cart2sph(x, y, z);
  bdy1 = (sin(phi) < -zt) & (z<0);
  cpx0 = R*cos(zt)*cos(th1);
  cpy0 = R*cos(zt)*sin(th1);
   
  ind = find(bdy1==1);
  
  
  n_ind = length(ind);
  
  for ii = 1:n_ind  %n_ind
  
    temp_dist = sqrt( (x(ind(ii))-cpx0(ind(ii))).^2 + (y(ind(ii))-cpy0(ind(ii))).^2 + (z(ind(ii))+zt).^2 );
      dist(ind(ii)) = temp_dist; 
      cpx(ind(ii)) = cpx0(ind(ii));
      cpy(ind(ii)) = cpy0(ind(ii));
      cpz(ind(ii)) = -zt;
      
  end
  
 
  
  % for points above zt = 0.25
   
  zt = 0.25 ;
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
  

  
  
  bdy1 = (bdy1|bdy2);
  
  
 
  
  
  % shift back to center
  cpx = cpx + xc;
  cpy = cpy + yc;
  cpz = cpz + zc;
