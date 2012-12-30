function [pass, str] = test_normals_pig()
  str = 'normals test 3d: build oriented normals for the cp pig';


  addpath('../../surfaces/readply');
  addpath('../../surfaces/tri');
  addpath('../../surfaces/tri2cp');

  %% Construct a grid in the embedding space
  %dx = 0.1;  % fails, too coarse
  dx = 0.05;
  x1d = (-2:dx:2)';
  y1d = x1d;
  z1d = x1d;
  nx=length(x1d);
  ny=length(y1d);
  nz=length(z1d);


  PlyFile = 'pig_loop2.ply';
  %PlyFile = 'annies_pig.ply';
  %PlyFile = 'bumpy_torus_scaled.ply';
  disp( ['reading triangulation from "' PlyFile '"'] );
  [Faces, Vertices] = plyread(PlyFile, 'tri');

  disp('converting to closest point representation');
  [IJK,DIST,CP,XYZ] = tri2cp(Faces, Vertices, dx, x1d(1), 5, 1);
  i = IJK(:,1);
  j = IJK(:,2);
  k = IJK(:,3);
  dist = DIST;
  cpxg = CP(:,1);
  cpyg = CP(:,2);
  cpzg = CP(:,3);
  xg = XYZ(:,1);
  yg = XYZ(:,2);
  zg = XYZ(:,3);

  % here is one place where meshgrid comes in: note ordering here.
  band = sub2ind([ny,nx,nz], j,i,k);

  p = 3;  % TODO: also try 1,5
  E = interp3_matrix(x1d,y1d,z1d, cpxg, cpyg, cpzg, p, band);

  make_plot = 0;
  if make_plot
    figure(11); clf;
    trisurf(Faces, Vertices(:,1), Vertices(:,2), Vertices(:,3))
    axis equal
    axis tight
    shading flat
    alpha(0.7);
    hold on;
  end

  dpt = (xg - 0).^2 + (yg - 0.4).^2 + (zg - 0).^2;
  [temp,seedin] = min(dpt);
  if temp > 2*dx
    error('should be closer');
  end
  [temp,seedout] = max(zg);

  pass = [~isempty(seedin) ~isempty(seedout)];


  inside = orientation_fill(xg,yg,zg,dist,dx,seedin,1);
  outside = orientation_fill(xg,yg,zg,dist,dx,seedout,1);
  [in2,out2,unknown] = orientation_stage2(xg,yg,zg,cpxg,cpyg,cpzg,dist,dx,E,inside,outside,1);
  pass = [pass isempty(unknown)];

  [insideg,sdist,unknown] = orientation_from_cp(xg,yg,zg,cpxg,cpyg,cpzg,...
                                                dist,dx,E, ...
                                                seedin, seedout, 1);

  pass = [pass isempty(unknown)];

  [n1,n2,n3] = normals_from_cp(xg,yg,zg, cpxg,cpyg,cpzg,dist, dx, 1);

  [n1,n2,n3] = normals_from_cp(xg,yg,zg, cpxg,cpyg,cpzg,sdist, dx, 1);


  if make_plot  % draw a hairy pig
  %for i = 1:length(xg)
  for ii=2000
    i = round(length(xg)*rand)
    plot3([cpxg(i) cpxg(i)+n1(i)/20], [cpyg(i) cpyg(i)+n2(i)/20], ...
          [cpzg(i) cpzg(i)+n3(i)/20], 'k-')
    end
  end

