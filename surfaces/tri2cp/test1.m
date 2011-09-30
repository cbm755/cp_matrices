addpath('../../cp_matrices');
addpath('../tri');
addpath('../readply');


dx = 0.025;

% ply file contains the triangles
disp('reading plyread');
%PlyFile = 'bunny.ply';
PlyFile = 'pig_loop2.ply';
%PlyFile = 'annies_pig.ply';
[Faces, Vertices] = plyread(PlyFile, 'tri');

disp('running tri2cp');
[IJK,DIST,CP,XYZ] = tri2cp2(Faces, Vertices, dx, -2);
%[IJK,DIST,CP,XYZ] = tri2cp(Faces, Vertices, dx, -2);
