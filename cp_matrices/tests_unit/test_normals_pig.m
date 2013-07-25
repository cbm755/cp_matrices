function [pass, str] = test_normals_pig()
  str = 'normals test 3d: build oriented normals for the cp pig';

  if isoctave()
    warning('no tri2cp yet for octave, failing');
    pass = 0;
    return
  end

  verb = 1;
  make_plot = 0;
  wh = 1;  % also supports genus3 but its quite thin and needs a
           % very high resolution


  %% Construct a grid in the embedding space
  %dx = 0.1;  % fails, too coarse
  %dx = 0.025;  % has trouble around the tail
  dx = 0.05;
  x1d = (-2:dx:2)';
  y1d = x1d;
  z1d = x1d;
  nx=length(x1d);
  ny=length(y1d);
  nz=length(z1d);

  if wh == 1
    if isoctave()
      load data_pig  % plyready() doesn't work in Octave
    else % matlab
      PlyFile = 'pig_loop2.ply';
      disp( ['reading triangulation from "' PlyFile '"'] );
      [Faces, Vertices] = plyread(PlyFile, 'tri');
    end
  else
    PlyFile = 'genus3.off';
    disp( ['reading triangulation from "' PlyFile '"'] );
    [Faces, Vertices] = offread(PlyFile);
  end

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


  if make_plot
    figure(11); clf;
    trisurf(Faces, Vertices(:,1), Vertices(:,2), Vertices(:,3))
    axis equal
    axis tight
    shading flat
    alpha(0.7);
    xlabel('x'); ylabel('y'); zlabel('z');
    hold on;
  end

  if wh == 1
    nrpt = [0 0.4 0];
  else
    nrpt = [0.4 -0.25 -0.2];
  end
  dpt = (xg - nrpt(1)).^2 + (yg - nrpt(2)).^2 + (zg - nrpt(3)).^2;
  [temp,seedin] = min(dpt);
  if temp > 2*dx
    error('should be closer');
  end

  [temp,seedout] = max(zg);

  pass = [~isempty(seedin) ~isempty(seedout)];


  inside = orientation_fill(xg,yg,zg,dist,dx,seedin,verb);

  if make_plot
    figure(12); clf;
    H = plot3(xg(inside),yg(inside),zg(inside),'k.');
    axis equal
    xlabel('x'); ylabel('y'); zlabel('z');
  end

  outside = orientation_fill(xg,yg,zg,dist,dx,seedout,verb);
  tic
  [in2,out2,unknown] = orientation_stage2(...
      xg,yg,zg,cpxg,cpyg,cpzg,dist,dx,E,inside,outside,verb);
  toc
  pass = [pass isempty(unknown)];

  % testing another approach
  % TODO: lots of different things we could do, all unreliable
  % (probably) in some specific geometry.  Clean up, pick one, etc.
  tic
  [in3,out3,un3] = orientation_stage2_global1(...
      xg,yg,zg,cpxg,cpyg,cpzg,dist,dx,E,inside,outside,verb);
  toc
  tic
  [in4,out4,un4] = orientation_stage2_global2(...
      xg,yg,zg,cpxg,cpyg,cpzg,dist,dx,E,inside,outside,verb);
  toc
  pass = [pass isempty(un3) isempty(un4)]



  [insideg,sdist,unknown] = orientation_from_cp(...
      xg,yg,zg,cpxg,cpyg,cpzg,dist,dx,E,seedin, seedout, verb);

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

