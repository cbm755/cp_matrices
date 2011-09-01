function [Axb,Axf,Ayb,Ayf,Azb,Azf] = average_upw1_3d_matrices(x,y,z, band1, band2, use_ndgrid)
% AVERAGE_UPW1_3D_MATRICES builds a matrix which performs arithmetic
% averaging
%   [AXB,AXF,AYB,AYF,AZB,AZF] = AVERAGE_UPW1_3D_MATRICES(X,Y,Z,BAND) returns 
%   the arithmetic averaging in the backwards x-direction, forwards
%   x-direction, and so on.
%
%   [AXB,AXF,AYB,AYF,AZB,AZF] = AVERAGE_UPW1_3D_MATRICES(X,Y,Z,BAND1,BAND2) 
%   performs the same with an optional outer band.
%
%   [AXB,AXF,AYB,AYF,AZB,AZF] = AVERAGE_UPW1_3D_MATRICES(X,Y,Z,BAND1,BAND2,USE_NDGRID)
%   with USE_NDGRID set as 'true' uses ndgrid ordering instead of meshgrid
%   ordering.
%
%   For example, averaging a quantity C in the forward x-direction would be 
%   ( C_ijk + C_(i+1)jk ) / 2.

  if (nargin <= 5)
    use_ndgrid = false;
  end
  if (nargin <= 4)
    band2 = band1;
  end

  % input checking
  [temp1, temp2] = size(x);
  if ~(  (ndims(x) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('x must be a vector, not e.g., meshgrid output');
  end
  [temp1, temp2] = size(y);
  if ~(  (ndims(y) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('y must be a vector, not e.g., meshgrid output');
  end
  [temp1, temp2] = size(z);
  if ~(  (ndims(z) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('z must be a vector, not e.g., meshgrid output');
  end

  weights = [1  1] / 2;
  PTS = [-1   0   0; ...
          0   0   0];
  Axb = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);

  weights = [1  1] / 2;
  PTS = [ 0   0   0; ...
          1   0   0];
  Axf = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);

  weights = [1  1] / 2;
  PTS = [ 0  -1   0; ...
          0   0   0];
  Ayb = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);

  weights = [1  1] / 2;
  PTS = [ 0   0   0; ...
          0   1   0];
  Ayf = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);

  weights = [1  1] / 2;
  PTS = [ 0   0  -1; ...
          0   0   0];
  Azb = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);

  weights = [1  1] / 2;
  PTS = [ 0   0   0; ...
          0   0   1];
  Azf = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);
