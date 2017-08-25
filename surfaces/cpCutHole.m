function [cpx,cpy,cpz,dist,bdy,SA_final] = cpVase_w_holesSA(x, y, z, ...
    lim,ab,cen,holes,tol,SA_wanted,dx,bw)
%cpf, holes,...

% x,y,z coordinate matrices of any same dimension
% lim = [angle from xy-plane, angle from xy-plane], a range of angle.
%   e.g. lim = [-pi/2, pi/2] makes a full ellipsoid


% just project onto vase
%cpf = @(x,y,z) cpVase(x, y, z, lim, ab, cen);
%[cpx,cpy,cpz,~,~] = cpf;

[cpx,cpy,cpz,~,~] = cpVase(x, y, z, lim, ab, cen);
bdy = zeros(size(cpx));
size_holes=size(holes);
num_holes=size_holes(1);

for k=1:num_holes
    hole = holes(k,:);
    
    xh = hole(1);
    yh = hole(2);
    zh = hole(3);
    Rh = hole(4);
    
    if exist('SA_wanted','var')
        
        SA_wanted_k = SA_wanted(k);
        R_bar = sqrt(SA_wanted_k/pi);
        RR = linspace(0.8*R_bar,1.2*R_bar, 50);
        for i = 1:length(RR)          
            I = ((cpx-xh).^2+(cpy-yh).^2+(cpz-zh).^2 < RR(i)^2);
            SASA(i) = sum(I(:))*(dx)^2/bw;
        end
          [~,II] = min(abs(SASA./SA_wanted_k-1));
          R_H = RR(II);
          SA_final = SASA(II);
              
    else
        I = ((cpx-xh).^2+(cpy-yh).^2+(cpz-zh).^2 < Rh^2);
        SA_final = 0;
    end
    
    
    
    % l for loop, points to project onto boundary of hole on surface
    xl = cpx(I);
    yl = cpy(I);
    zl = cpz(I);
    
    
    sdist = 1;
    
    tol2 = 0.01;
    
    while  abs(sdist) > tol | abs(dist1) > tol
       
        % project onto sphere "hole"
        [xl2,yl2,zl2,sdist] = cpSphere(xl,yl,zl,Rh,[xh,yh,zh]);
                        
        % project onto surface
        %[xll,yll,zll,dist1] = cpFromTriSlow(xl2, yl2, zl2, Faces, Vertices);
        [xll,yll,zll,dist1,~] = cpVase(xl2, yl2, zl2, lim, ab, cen);
         xl=xll;
         yl=yll;
         zl=zll;
       
        % avoid projecting back and forth infinite times
        dist2=sqrt((xl-xll).^2 + (yl-yll).^2 + (zl-zll).^2);
        dist3=sqrt((xl-xl2).^2 + (yl-yl2).^2 + (zl-zl2).^2);
        if  dist2  < tol2 & dist3 > Rh/2
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
