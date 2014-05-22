function [] = test_positivity()
dx = 1;
M = 100;
a = rand(M,1);   a = a*dx/2;
b = rand(M,1);   b = b*dx/2;
c = rand(M,1);   c = c*dx/2;

% a = sort(a);
% b = dx/4;
% c = dx/8;


gamma1 = 48*dx^4 ./ ( (dx+a).*(2*dx-a).*(dx+b).*(2*dx-b).*(dx+c).*(2*dx-c) );

X = (dx+a).*a.*(dx-a).*(2*dx-a);
Y = (dx+b).*b.*(dx-b).*(2*dx-b);
Z = (dx+c).*c.*(dx-c).*(2*dx-c);

d = ( a.*(dx-b).*(dx-c) + b.*(dx-a).*(dx-c) + c.*(dx-a).*(dx-b) )...
    - ( (dx-2*a).*(dx-b).*(dx-c) + (dx-2*b).*(dx-a).*(dx-c) + (dx-2*c).*(dx-a).*(dx-b) )...
    .* gamma1/(4*dx^2)./(1./X + 1./Y + 1./Z);


d1 = ( a.*(dx-b).*(dx-c) + b.*(dx-a).*(dx-c) + c.*(dx-a).*(dx-b) ) .* ( X.*Y + Y.*Z + Z.*X )...
    - 12*dx^2*(dx-a).*(dx-b).*(dx-c).*...
              ( a.*(dx-2*a).*b.*(dx-b).*c.*(dx-c) + ...
                a.*(dx-a).*b.*(dx-2*b).*c.*(dx-c) + ...
                a.*(dx-a).*b.*(dx-b).*c.*(dx-2*c) );

d2 = ( a.*(dx-b).*(dx-c) + b.*(dx-a).*(dx-c) + c.*(dx-a).*(dx-b) ) ...
    .* (X.*Y + Y.*Z + Z.*X + 12*dx^2*(a.*b.*c).*(dx-a).*(dx-b).*(dx-c))...         
       - 36*dx^2*a.*b.*c.*((dx-a).*(dx-b).*(dx-c)).^2;
            
   
d3 = 3* (a.*(dx-b).*(dx-c).*Y.*Z + b.*(dx-a).*(dx-c).*Z.*X + c.*(dx-a).*(dx-b).*X.*Y) ...
    - 12*dx^2*(dx-a).*(dx-b).*(dx-c).*...
              ( a.*(dx-2*a).*b.*(dx-b).*c.*(dx-c) + ...
                a.*(dx-a).*b.*(dx-2*b).*c.*(dx-c) + ...
                a.*(dx-a).*b.*(dx-b).*c.*(dx-2*c) );
  
            
d4 =  3* (a.*(dx-b).*(dx-c).*Y.*Z + b.*(dx-a).*(dx-c).*Z.*X + c.*(dx-a).*(dx-b).*X.*Y) ...
      - ( a.*(dx-b).*(dx-c) + b.*(dx-a).*(dx-c) + c.*(dx-a).*(dx-b) ) .* ( X.*Y + Y.*Z + Z.*X );
           
  
d5 = a.*(dx-b).*(dx-c).*Y.*Z - 4*dx^2*(dx-a).*(dx-b).*(dx-c).*a.*(dx-2*a).*b.*(dx-b).*c.*(dx-c);
d6 = b.*(dx-a).*(dx-c).*Z.*X - 4*dx^2*(dx-a).*(dx-b).*(dx-c).*a.*(dx-a).*b.*(dx-2*b).*c.*(dx-c);
d7 = c.*(dx-a).*(dx-b).*X.*Y - 4*dx^2*(dx-a).*(dx-b).*(dx-c).*a.*(dx-a).*b.*(dx-b).*c.*(dx-2*c);

flag = d<0;
flag1 = d1 < 0;
flag2 = d2 < 0;
flag3 = d3 < 0;
flag4 = d4 < 0;
flag5 = d5 < 0;
flag6 = d6 < 0;
flag7 = d7 < 0;


nnz(flag)
nnz(flag1)
nnz(flag2)
nnz(flag3)
nnz(flag4)
nnz(flag5)
nnz(flag6)
nnz(flag7)

end


