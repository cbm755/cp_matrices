addpath('../../cp_matrices');
addpath('../tri');
addpath('../readply');


dx = 0.05;

% ply file contains the triangles
disp('reading plyread');
%PlyFile = 'bunny.ply';
PlyFile = 'pig_loop2.ply';
%PlyFile = 'annies_pig.ply';
[Faces, Vertices] = plyread(PlyFile, 'tri');

%mex tri2cp_helper.c

disp('running tri2cp');
[IJK,DIST,CP,XYZ,WHICH_FACES] = tri2cp(Faces, Vertices, dx, -2);
%[IJK2,DIST2,CP2,XYZ2] = tri2cp_old(Faces, Vertices, dx, -2);

%max(IJK-IJK2)
%max(DIST-DIST2)
%max(CP-CP2)
%max(XYZ-XYZ2)
