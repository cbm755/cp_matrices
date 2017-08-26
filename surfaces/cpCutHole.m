function [cpx,cpy,cpz,dist,bdy,SA_final] = cpCutHole(x,y,z,cpf,hole_cen,SA_target,tol,dx,bw,given_cp_data)
%cpCutHole  Cut holes in surfaces described by cpfun
%
% x,y,z coordinate matrices of any same dimension
%
% hole_cen = [hx hy hz] gives the center of the hole.  If this is not on
% the surface given by cpf, it will be projected onto the surface (by
% calling cpf).
%
% Multiple holes are supported: pass each center as a row of "hole_cen"
% and pass a vector to "SA_target".
%
% TODO: the grid in x y z must already be banded with bandwidth bw*dx.
% This is a restriction based on how we estimate the surface area
% integral.  It should be fixed using Kublic and Tsai or similar.

  if (nargin >= 10)
    cpx = given_cp_data{1};
    cpy = given_cp_data{2};
    cpz = given_cp_data{3};
  else
    [cpx, cpy, cpz, tilde] = cpf(x, y, z);
    %[cpx, cpy, cpz, tilde, tilde] = cpf(x, y, z);
  end
  %toc
  % TODO: if nargout is 5, get bdy from cpf instead
  %nargout(cpf)
  %tic
  bdy = zeros(size(cpx));
  [num_holes, temp] = size(hole_cen);
  assert (temp == 3)

  % below assumes cpx, cpy, cpz are long vectors
  savesz = size(cpx);
  cpx = cpx(:);
  cpy = cpy(:);
  cpz = cpz(:);
  bdy = bdy(:);

  for k=1:num_holes
    xh = hole_cen(k, 1)
    yh = hole_cen(k, 2)
    zh = hole_cen(k, 3)
    % project onto surface
    [xh yh zh] = cpf(xh, yh, zh)

    SA_wanted_k = SA_target(k)
    R_bar = sqrt(SA_wanted_k/pi)
    RR = linspace(0.75*R_bar, 1.05*R_bar, 100);
    for i = 1:length(RR)
      I = ((cpx-xh).^2+(cpy-yh).^2+(cpz-zh).^2 < RR(i)^2);
      assert(nnz(I) > 0, 'empty sphere/surface intersection!')
      SASA(i) = sum(I(:))*(dx)^2/(2*bw);
    end
    SASA./SA_wanted_k
    [tilde, i] = min(abs(SASA./SA_wanted_k-1))
    R_H = RR(i)
    assert (i > 1 && i < length(RR), 'should be in the middle')
    assert (tilde < 0.1)

    %else
    %    I = ((cpx-xh).^2+(cpy-yh).^2+(cpz-zh).^2 < Rh^2);
    %    SA_final = 0;
    %end

    % index of all points whose cp are within the sphere
    I = ((cpx-xh).^2 + (cpy-yh).^2 < R_H^2);


    % l for loop, points to project onto boundary of hole on surface
    xl = cpx(I);
    yl = cpy(I);
    zl = cpz(I);
    
    
    sdist = 1;
    
    tol2 = 0.01;
    
    while  abs(sdist) > tol | abs(dist1) > tol
       
        % project onto sphere "hole"
        [xl2,yl2,zl2,sdist] = cpSphere(xl,yl,zl,R_H,[xh,yh,zh]);
                        
        % project onto surface
        [xll,yll,zll,dist1] = cpf(xl2, yl2, zl2);
         xl=xll;
         yl=yll;
         zl=zll;

        % avoid projecting back and forth infinite times
        dist2=sqrt((xl-xll).^2 + (yl-yll).^2 + (zl-zll).^2);
        dist3=sqrt((xl-xl2).^2 + (yl-yl2).^2 + (zl-zl2).^2);
        if  dist2  < tol2 & dist3 > R_H/2
           xll = xll+10^(-5)*rand; % add perturbation
           yll = yll+10^(-5)*rand;
           zll = zll+10^(-5)*rand;

        end
        
        xl=xll;
        yl=yll;
        zl=zll; 
          
    end
       
    cpx(I) = xl;
    cpy(I) = yl;
    cpz(I) = zl;
    
    bdy = bdy | I; 
end


  cpx = reshape(cpx, savesz);
  cpy = reshape(cpy, savesz);
  cpz = reshape(cpz, savesz);
  bdy = reshape(bdy, savesz);
dist = sqrt((cpx-x).^2 + (cpy-y).^2 + (cpz-z).^2);
end
