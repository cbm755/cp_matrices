function [cpx,cpy,cpz,dist] = cpFromTriSlow(xx, yy, zz, faces, vertices)
%
%   A global search: for a given grid point, search all triangles one
%   at a time.
%
%   Note the various subroutines here are just copied from the tri2cp code.
%
% TODO: this could be made much faster with a tree search on vertices to
%       find nearby triangles.
%
% TODO: not vectorized, just loops on multiple inputs.

  number_faces = size(faces, 1);

  cpx = zeros(size(xx));
  cpy = cpx;
  cpz = cpx;
  dist = cpx;
  
  for j=1:length(xx(:))
    x = xx(j);
    y = yy(j);
    z = zz(j);
    [cx, cy, cz, dd_min] = FindClosestPointToOneTri(x,y,z, 1, faces, vertices);
    for i=2:number_faces
      [t1, t2, t3, dd] = FindClosestPointToOneTri(x,y,z, i, faces, vertices);
      if (dd<dd_min)
        dd_min = dd;
        cx = t1;
        cy = t2;
        cz = t3;
      end
    end
    cpx(j) = cx;
    cpy(j) = cy;
    cpz(j) = cz;
    dist(j) = sqrt(dd_min);
  end
end


function [c1, c2, c3] = ProjectOnSegment(c1, c2, c3, p1, p2, p3, q1, q2, q3)
%
% Project the point (c1,c2,c3) onto the line segment specified by
% (p1,p2,p3) and (q1,q2,q3).
% (code by Steve Ruuth)

  cmp1 = c1-p1;
  cmp2 = c2-p2;
  cmp3 = c3-p3;
  qmp1 = q1-p1;
  qmp2 = q2-p2;
  qmp3 = q3-p3;

  lambda = (cmp1*qmp1+cmp2*qmp2+cmp3*qmp3)/(qmp1*qmp1+qmp2*qmp2+qmp3*qmp3);
  lambda_star = max(0.0, min(lambda, 1.0));

  c1 = p1+lambda_star*qmp1;
  c2 = p2+lambda_star*qmp2;
  c3 = p3+lambda_star*qmp3;
end


function [c1, c2, c3, dd] = FindClosestPointToOneTri(a1, a2, a3, face_index, face, vertex)
%/*
% * Closest point and distance from a point (a1,a2,a3) to a triangle
% * indexed by face_index.  Returns the *squared* distance and the
% * closest point in (c1,c2,c3).  Uses global vars `face' and
% * `vertex'.  (Based on code by Steve Ruuth)
% * TODO: this code may not be robust to degenerate triangles (line
% * segments and points).  More testing required.
% */

  %/* obtain the indices to the three vertices */
  index_p = face(face_index, 1);
  index_q = face(face_index, 2);
  index_r = face(face_index, 3);

  %/* translate so the p is at the origin */
  a1 = a1 - vertex(index_p, 1);
  a2 = a2 - vertex(index_p, 2);
  a3 = a3 - vertex(index_p, 3);
  q1 = vertex(index_q, 1)-vertex(index_p, 1);
  q2 = vertex(index_q, 2)-vertex(index_p, 2);
  q3 = vertex(index_q, 3)-vertex(index_p, 3);
  r1 = vertex(index_r, 1)-vertex(index_p, 1);
  r2 = vertex(index_r, 2)-vertex(index_p, 2);
  r3 = vertex(index_r, 3)-vertex(index_p, 3);

  %/* evaluate the various matrix entries */
  a11 = q1*q1+q2*q2+q3*q3;
  a12 = q1*r1+q2*r2+q3*r3;
  a22 = r1*r1+r2*r2+r3*r3;
  b1  = a1*q1+a2*q2+a3*q3;
  b2  = a1*r1+a2*r2+a3*r3;

  %/* find the inverse matrix and solve for lambda and mu */
  factor = 1/(a11*a22-a12*a12);
  i11 = a22*factor;
  i12 =-a12*factor;
  i22 = a11*factor;
  lambda = i11*b1+i12*b2;
  mu     = i12*b1+i22*b2;
  c1 = lambda*q1+mu*r1;
  c2 = lambda*q2+mu*r2;
  c3 = lambda*q3+mu*r3;

  if ((lambda<0) && (mu<0) && (lambda+mu<=1))
    c1 = 0; c2 = 0; c3 = 0;
  elseif ((lambda>=0) && (mu<0) && (lambda+mu<=1))
    [c1,c2,c3] = ProjectOnSegment(c1,c2,c3,0,0,0,q1,q2,q3);
  elseif ((lambda>=0) && (mu<0) && (lambda+mu>1))
    c1 = q1;
    c2 = q2;
    c3 = q3;
  elseif ((lambda>=0) && (mu>=0) && (lambda+mu>1))
    [c1, c2, c3] = ProjectOnSegment(c1,c2,c3,q1,q2,q3,r1,r2,r3);
  elseif ((lambda<0) && (mu>=0) && (lambda+mu>1))
    c1 = r1;
    c2 = r2;
    c3 = r3;
  elseif ((lambda<0) && (mu>=0) && (lambda+mu<=1))
    [c1, c2, c3] = ProjectOnSegment(c1,c2,c3,r1,r2,r3,0,0,0);
  elseif ((lambda>=0) && (mu>=0) && (lambda+mu<=1))
    %% TODO: do nothing, what case is this?
  else
    fprintf(1, 'Error: non-enumerated case, can this happen?\n');
    fprintf(1, 'lambda mu %g %g\n', lambda, mu);
    fprintf(1, 'factor %g\n', factor);
    fprintf(1, 'a11,a22,a12=%g,%g,%g\n', a11, a22, a12);
    fprintf(1, 'det?=%g\n', a11*a22 - a12*a12);
    error('Error in CP to tri, unanticipated case');
  end

  %/* Calculate distance */
  %/* Note: dd is dist squared! */
  %/* TODO: HORRIBLE HACK 2010-07-28, this was to deal with a ply file
  %   with degenerate triangles. */
  %/*if (isinf(factor)) {*/
  %/* TODO, just gets worse and worse, where is isinf on windows? */
  if ((factor > 1e20) || (factor < -1e20))
    dd = 10000.0;
    error('"factor" is infinite: panic!');
  else
    dd  = (a1-c1)*(a1-c1)+(a2-c2)*(a2-c2)+(a3-c3)*(a3-c3);
  end

  %/* Shift everything back */
  c1 = c1 + vertex(index_p, 1);
  c2 = c2 + vertex(index_p, 2);
  c3 = c3 + vertex(index_p, 3);

end
