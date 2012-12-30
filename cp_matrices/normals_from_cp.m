function [n1,n2,n3] = normals_from_cp(x,y,z,cpx,cpy,cpz,dist,dx,verbose)
%NORMALS_FROM_CD  Compute normals from a CP representations
%   Given a closest point representation of a codim-1 surface,
%   find normals.
%
%   [n1,n2,n3] = normals_from_cp(x,y,z,cpx,cpy,cpz,dist,dx,verbose)
%   [n1,n2] = normals_from_cp(x,y,cpx,cpy,dist,dx,verbose)
%
%   If 'dist' is signed distance (with negative inside), the code will
%   attempt to give you outward-pointing normals.  If its just
%   positive distance, the orientation will be undefined.
%
%   No accuracy is promised: the normal should be considered
%   qualitatively.
%
%   The computation is easy for points far from the surface: just n =
%   x - cp(x).  But this fails if x = cp(x) and probably unstable for
%   x close to cp(x).)  For such we points, we estimate a normal from
%   the neighbouring CPs.
%
%   TODO: for point clouds, averaging nbring normals probably
%   important.  But for circle/sphere/smooth surfaces if we
%   see zero-dist nbrs, we should just take that normal.
%   Perhaps need a toggle to select which one (via the weights).
%
%   TODO: opt structure?

  if (nargin == 9)
    dim = 3;
  elseif (nargin == 8)
    dim = 3;
    verbose = 0;
  elseif (nargin == 7)
    dim = 2;
    verbose = dist;
    dx = cpz;
    dist = cpy;
    cpy = cpx;
    cpx = z;
    clear z cpz
  elseif (nargin == 6)
    dim = 2;
    verbose = 0;
    dx = cpz;
    dist = cpy;
    cpy = cpx;
    cpx = z;
    clear z cpz
  else
    error('wrong inputs');
  end

  %% Parameters
  % how many neighbours near the CP
  nnbrs = 6;
  TOO_CLOSE_TOL = 100*eps;  % probably should match orientation_from_cp

  tooclose = find(abs(dist) <= TOO_CLOSE_TOL);
  ntc = length(tooclose);

  %% Find most normals with x - cp(x)
  % may give NaN's near surface but will be overwritten later
  n1 = x - cpx;
  n2 = y - cpy;
  if dim == 3
    n3 = z - cpz;
    len = sqrt(n1.*n1 + n2.*n2 + n3.*n3);
  else
    len = sqrt(n1.*n1 + n2.*n2);
  end
  % normalize the normals
  n1 = n1 ./ len;
  n2 = n2 ./ len;
  if dim == 3
    n3 = n3 ./ len;
  end

  %% check if we're oriented
  inside = dist < 0;
  if nnz(inside) >= 1
    isoriented = true;
  else
    isoriented = false;
  end

  if isoriented
    n1 = n1.*(-1*inside) + n1.*(~inside);
    n2 = n2.*(-1*inside) + n2.*(~inside);
    if dim == 3
      n3 = n3.*(-1*inside) + n3.*(~inside);
    end
  end

  % if it seems we have trouble finding enough neighbours, we could try
  % two passes, deferring on the first.

  for ii=1:ntc
    i = tooclose(ii);
    if (dim == 3)
      tcp = [cpx(i) cpy(i) cpz(i)];
      % could accelerate this search with a kdtree
      dd = (cpx(:)-tcp(1)).^2 + (cpy(:)-tcp(2)).^2 + (cpz(:)-tcp(3)).^2;
    elseif (dim == 2)
      tcp = [cpx(i) cpy(i)];
      dd = (cpx(:)-tcp(1)).^2 + (cpy(:)-tcp(2)).^2;
    end
    dd(i) = 1e42;  % big, so we don't find the current point
    [dd,II] = sort(dd, 'ascend');

    % want this but discarding ones that have a small dist
    %I = II(1:nnbrs);
    %d = dd(1:nnbrs);

    % a logical index: 1 if far enough away
    wh = dist(II) > TOO_CLOSE_TOL;
    % want the first nnbr nonzeros of this
    wh2 = find(wh, nnbrs);
    II = II(wh2);
    dd = dd(wh2);

    if any(dd > dx*dx)
      % want them close but not sure how much it matters
      warning('nbrs not close enough?');
      disp('dropping to keyboard');
      keyboard
    end

    ns = zeros(nnbrs, dim);
    ns(:,1) = n1(II);
    ns(:,2) = n2(II);
    if dim == 3
      ns(:,3) = n3(II);
    end
    if dim == 3
      ips = n1(II)*n1(II)' + n2(II)*n2(II)' + n3(II)*n3(II)';
    else
      ips = n1(II)*n1(II)' + n2(II)*n2(II)';
    end

    %% Weighting of the nbr normals
    % ultimately, this is not really much different from
    % "mean(ns,1)", but we might want to be more sophisticated
    % based on distance (or other properties potentially).
    % TODO: maybe users want to be able to switch weights with an option.

    % Weights fall off with 1/distance.
    % if any nbrs are (almost) the same, this selects their normals.
    dd2 = sqrt(dd);
    dd2(dd2 == 0) = eps;  % just to avoid inf
    w = 1./dd2;
    weights = w / sum(w);

    % Weight using a normal curve.
    % (scale by dx b/c reasonable to expect closer neighbours
    % with a smaller dx, recall dd is squared dist).
    %weights = exp(-dd/(dx*dx));

    if (1==0)
      % debug plotting in 2D
      figure(10); clf;
      plot(cpx(i),cpy(i),'k*');
      hold on;
      th = 0:(pi/1024):2*pi;
      plot(cos(th),sin(th),'k--')
      plot(cpx(II),cpy(II),'ro');
      for j=1:nnbrs
        I = II(j);
        H = plot([x(I) cpx(I)],[y(I) cpy(I)],'r-');
      end
      axis equal
      drawnow
    end

    if ~isoriented
      % TODO: maybe we should do a PCA of ips or something similar?  For
      % now, just assume they reasonably line up but with some facing
      % the other way.
      ipsc = abs(ips);
      % arbitrarily choose the orientation of the first nbr
      weights = weights .* (ips(1,:) > 0)' + ...
               -weights .* (ips(1,:) <= 0)';
    else
      ipsc = ips;
    end
    if all(all(ipsc > 0.6))
      tn = zeros(size(ns(1,:)));
      for j=1:nnbrs
        tn = tn + weights(j)*ns(j,:);
      end
    else
      dd'
      ns'
      ips
      warning('more sophistication needed?')
      disp('dropping to keyboard');
      keyboard
    end

    len = tn*tn';

    if len < 0.5
      warning('normal very far from normal, weight problem?');
      disp('dropping to keyboard');
      keyboard
    end

    n1(i) = tn(1) / len;
    n2(i) = tn(2) / len;
    if dim == 3
      n3(i) = tn(3) / len;
    end

    %keyboard

    if (verbose>=2)
      dd'
      ns'
      weights'
      tn
      mean(ns',2)'
    end

    if (verbose)
      fprintf('#%d dist=%9.3g, %.2g, [%.2g,%.2g]\n', ...
              i,dist(i), max(dd), min(min(ips)),max(max(ips)));
    end
    %keyboard
    %pause
  end

