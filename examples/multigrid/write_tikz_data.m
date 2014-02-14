% Write tikz data for illustration of restriction and prolongation operator
% for multigrid closest point method, the output .dat file is just the x & y
% coordinates of all the grid points

addpath('../../cp_matrices');
addpath('../../surfaces');

x0 = -2;
x1 = 2;
y0 = -2;
y1 = 2;

% grid size of the coarse grid
dxc = 0.2;
x1dc = (x0:dxc:x1)';
y1dc = x1dc;
% grid size of the fine grid
dxf = 0.1;
x1df = (x0:dxf:x1)';
y1df = x1df;

dim = 2;
p = 1;
order = 2;

bw = 1.0002*sqrt( (dim-1)*((p+1)/2)^2 + (order/2+(p+1)/2)^2  );

cpf = @cpSemicircle;
[xxc yyc] = meshgrid(x1dc, y1dc);
[cpxc, cpyc, distc] = cpf(xxc,yyc);
bandc = find(abs(distc) <= bw*dxc);
xc = xxc(bandc); yc = yyc(bandc);
flag = (xc > 0) & (yc > 0);
xc = xc(flag);  yc = yc(flag);

[xxf yyf] = meshgrid(x1df, y1df);
[cpxf, cpyf, distf] = cpf(xxf,yyf);
bandf = find(abs(distf) <= bw*dxf);
xf = xxf(bandf); yf = yyf(bandf);
flag = (xf > 0) & (yf > 0);
xf = xf(flag);  yf = yf(flag);

plot(xc, yc, 'bo');
hold on, plot(xf, yf, 'r.');
axis equal

file = fopen('coarse.dat','w');
formatSpec = '%1.1f %1.1f \n';
fprintf(file, formatSpec,[xc'; yc']);
% for i = 1:1:length(xc)
%     fprintf(file, '%1.1f %1.1f\n',xc(i),yc(i));
% end
fclose(file);

file = fopen('fine.dat','w');
formatSpec = '%1.1f %1.1f \n';
fprintf(file, formatSpec,[xf'; yf']);
% for i = 1:1:length(xf)
%     fprintf(file, '%1.1f %1.1f\n',xf(i),yf(i));
% end
fclose(file);
