function [Faces, Vertices] = offread(fname)
%OFFREAD Read a triangulation from an ascii .off file
%   [Faces, Vertices] = offread('genus3.off');
%
%   Note: this may not be a standard off file.

  f = fopen(fname);

  s = fgetl(f);
  assert(strcmpi(s, 'OFF'));

  [temp] = fscanf(f, '%d %d', 2);
  numv = temp(1);
  numf = temp(2);
  s = fgetl(f);  % get the rest of the numv numf 0 line

  s = fgetl(f);  % a blank line?
  assert(strcmp(s, ''));

  Vertices = fscanf(f, '%g', [3, numv]);
  Faces = fscanf(f, '%g', [4, numf]);

  fclose(f);

  %  3  ix iy iz
  if any(Faces(1,:) ~= 3)
    error('non-triangle in .off file')
  end

  Vertices = transpose(Vertices);

  % discard triangles columns and c-based indexing to matlab
  Faces = transpose(Faces(2:4,:)) + 1;

  % debugging plot
  if (1==0)
    xp = Vertices(:,1);
    yp = Vertices(:,2);
    zp = Vertices(:,3);

    up = sin(10*zp).*cos(7*xp).*sin(6*yp);

    figure(1); clf;
    trisurf(Faces,xp,yp,zp, up);
    camlight left
  end


