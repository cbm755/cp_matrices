%% Load tri data for viz
addpath('../readply')

[F,V] = plyread('bumpy_torus_scaled.ply','tri');

% vertices of the triangles, for plotting
Xp = V(:,1);
Yp = V(:,2);
Zp = V(:,3);
up = sin(8*Zp);  % use interp to get this from 3D grid values
figure(1); clf;
trisurf(F, Xp, Yp, Zp, up)
shading interp
axis equal
hold on;
pause(0)


%% Load CP band data
dx = 0.1;
p = 3;
grid_data = load( ['bumpy_torus_griddata_p' num2str(p) '_dx' num2str(dx) '.txt' ]);


numband = size(grid_data,1);
ig = grid_data(:,1)+1; % +1 to adjust C to Matlab indexes
jg = grid_data(:,2)+1;
kg = grid_data(:,3)+1;
dg = grid_data(:,4);
cpxg = grid_data(:,5);
cpyg = grid_data(:,6);
cpzg = grid_data(:,7);
xg = grid_data(:,8);
yg = grid_data(:,9);
zg = grid_data(:,10);

%% Put the banded data into 3D arrays
x = -2:dx:2;
y = x;   z = x;
[xx,yy,zz] = meshgrid(x,y,z);
% shouldn't matter what the cp is if its outside the band, so set them
% all to zero then loop over the band
cpx = zeros(size(xx));
cpy = zeros(size(xx));
cpz = zeros(size(xx));


for c=1:numband
  i=ig(c);  j=jg(c);  k=kg(c);
  cpx(j,i,k) = cpxg(c); % i,j swapped b/c of meshgrid
  cpy(j,i,k) = cpyg(c);
  cpz(j,i,k) = cpzg(c);


  check = norm([xx(j,i,k)-xg(c) yy(j,i,k)-yg(c) zz(j,i,k)-zg(c)]);
  if check > 10*eps
    c, check
    error('detecting indexing issue')
  end
end

%% plot some of the x-CP vectors
for c=round(rand(1,500)*numband)
  plot3([xg(c) cpxg(c)], [yg(c) cpyg(c)], [zg(c) cpzg(c)], 'k-');
end
pause(0)
