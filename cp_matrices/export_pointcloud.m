%i=4;
%eigenvec = V(:,I(i));
%u = eigenvec
RGB = scalar2rgb(u, [], [], true);

normals = [cpxg - xg, cpyg - yg, cpzg - zg];
%normals = normals ./ sqrt(normals(:,1).^2 + normals(:,2).^2 + normals(:,3).^2);
len = sqrt(normals(:,1).^2 + normals(:,2).^2 + normals(:,3).^2);
min(len)
normals(:,1) = normals(:,1) ./ len;
normals(:,2) = normals(:,2) ./ len;
normals(:,3) = normals(:,3) ./ len;

%ASC=[20*[xg yg zg]  RGB  normals];
ASC=[15*[cpxg cpyg cpzg]  RGB  normals];


f = fopen('pig_cps_w_normals.asc','w');
fprintf(f,'%f %f %f %d %d %d %d %d %d\n', ASC');
fclose(f);
