function [cpx,cpy,cpz,dist,bdy,SA_final] = cpCutHole(x,y,z,cpf,hole_cen,SA_target,tol,dx,bw,cpx0,cpy0,cpz0)
%cpCutHole  Cut holes in surfaces described by cpfun
%
% x,y,z coordinate matrices of any same dimension
%
% hole_cen = [hx hy hz] gives the center of the hole.  If this is not on
% the surface given by cpf, it will be projected onto the surface (by
% calling cpf.
%
% Multiple holes are supported: pass each center as a row of "hole_cen"
% and pass a vector to "SA_target"

  warning('WIP')
  % TODO: if nargout is 5, get bdy from cpf instead
  %nargout(cpf)
  %tic
  % TODO: at least for triangulated surface, need a faster approach
  % TODO: we could optionally pass cpx0 etc as previously computed cp data.
  [cpx, cpy, cpz, tilde] = cpf(x, y, z);
  %toc
  %[cpx, cpy, cpz, tilde, tilde] = cpf(x, y, z);
  bdy = zeros(size(cpx));
  [num_holes, temp] = size(hole_cen);
  assert (temp == 3)

  for k=1:num_holes
    hole = hole_cen(k,:);

    xh = hole(1);
    yh = hole(2);
    zh = hole(3);


        
        SA_wanted_k = SA_target(k);
        R_bar = sqrt(SA_wanted_k/pi);
        RR = linspace(0.8*R_bar,1.2*R_bar, 50);
        for i = 1:length(RR)          
            I = ((cpx-xh).^2+(cpy-yh).^2+(cpz-zh).^2 < RR(i)^2);
            SASA(i) = sum(I(:))*(dx)^2/bw;
        end
          [tilde,II] = min(abs(SASA./SA_wanted_k-1));
          R_H = RR(II);
          SA_final = SASA(II);
              
    %else
    %    I = ((cpx-xh).^2+(cpy-yh).^2+(cpz-zh).^2 < Rh^2);
    %    SA_final = 0;
    %end
    
    
    
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
        %[xll,yll,zll,dist1,tilde] = cpVase(xl2, yl2, zl2, lim, ab, cen);
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


dist = sqrt((cpx-x).^2 + (cpy-y).^2 + (cpz-z).^2);
end
