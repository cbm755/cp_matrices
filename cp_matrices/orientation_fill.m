function inside = orientation_fill(xx,yy,zz,dist,dx,seeds,verbose)
%ORIENTATION_FILL  Orient a surface distance fcn by flood filling
%   inside = orientation_fill(xx,yy,zz, dist, dx, seeds, verbose)
%   This code starts at the points in 'seeds' and flood fills
%   toward the surface.  'seeds' should be list of indices into the
%   arrays xx, yy, zz, dist.  'inside' will be the same shape as xx
%   and contains ones for each points classified as inside.  Set
%   verbose to 2 to show isosurface of the progress.
%
%   This code it not guaranteed to find all inside points.  Its just
%   helper code for a larger routine (orientation_from_cp).  You
%   probably want that instead.

  inside = zeros(size(xx));
  looked = logical(inside);

  if length(seeds) ~= 1
    %warning('more than one seed might not be optimal');
  end
  inside(seeds) = 1;

  if verbose == 2
    figure(10); clf;
  end

  c = 0;
  while(1)
    c = c + 1;
    % look for points inside but not yet looked at (looked)
    I = find(inside & (~looked));

    if (isempty(I))
      break
    end
    % chould just take the first one, but instead we sort and
    % choose the furthest
    %i = I(1);
    %i = I(ceil(rand*length(I)));
    dd = abs(dist(I));
    [temp,ii] = sort(dd,'descend');
    i = I(ii(1));
    d = abs(dist(i));

    % now everything within dist is inside too
    x = xx(i);  y = yy(i);  z = zz(i);
    %[x y z d]
    % todo: add a small tolerance to prevent leakage here?
    tol = 0.5*dx;
    II = find( (xx-x).^2 + (yy-y).^2 + (zz-z).^2 < (d-tol)^2);
    inside(II) = 1;

    %looked(i) = 1;
    % the middle of the sphere can be considered done, but we need a
    % small band at the edge to seed further fill
    % TODO: -1.5 here is tunable parameter, larger values seem to
    % need less stage2, but how much slower is this stage?
    pd = d - 1.5*dx;
    if (pd > 0)
      II = find( (xx-x).^2 + (yy-y).^2 + (zz-z).^2 < pd^2 );
      if ~isempty(II)
        looked(II) = 1;
        text = ['looked at ball:' num2str(length(II))];
      else
        % none found in the ball, doubt this can happen
        looked(i) = 1;
        text = 'none found: looked at point';
      end
    else
      text = 'small rad, looked at point';
      looked(i) = 1;
    end

    if verbose == 2 && ((c < 30) || (mod(c,100) == 0))
      set(0, 'CurrentFigure', 10);  clf;
      isosurface(xx,yy,zz,inside,0.5);
      axis equal
      title([ 'inside, iteration ' num2str(c) ]);
      alpha(0.7)
      drawnow()
    end

    if (verbose >= 1)
      fprintf('iter%d, #look,in=%d,%d, %s', c, ...
              nnz(looked), nnz(inside), text);
    end
  end

