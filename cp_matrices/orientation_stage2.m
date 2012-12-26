function [in2,out2] = orientation_stage2(xx,yy,zz,cpx,cpy,cpz,dist,dx,inside,outside,verbose)
%ORIENTATION_STAGE2  Helper function for ORIENTATION_FROM_CP
%   Classify a few closest points as inside/outside supposing that
%   most are already known.
%
%   Uses a local calculation based on neighbouring CPs.  You can
%   probably invent (fine-scale) geometry that will fool this.

  noclass = find(~inside & ~outside);

  nnc = length(noclass);

  in2 = inside;
  out2 = outside;

  % TODO: maybe better to have more than one pass: could defer?
  for ii=1:nnc
    i = noclass(ii);
    tcp = [cpx(i) cpy(i) cpz(i)];
    % could accelerate this search with a kdtree
    dd = (cpx(:)-tcp(1)).^2 + (cpy(:)-tcp(2)).^2 + (cpz(:)-tcp(3)).^2;
    [dd,II] = sort(dd, 'ascend');

    % take 4 nbring closest points, TODO: lots of choice here
    nnbrs = 6;
    % Need the ones that are not the same cp! (why?)
    %III=find(dd >= 10*eps);
    %I = II(III(1:nnbrs));
    %d = dd(III(1:nnbrs));
    I = II(1:nnbrs);
    d = dd(1:nnbrs);
    if (any(d > dx))
      % want them close but not sure how much it matters
      warning('nbrs not close enough?');
      disp('dropping to keyboard');
      keyboard
    end
    A = I(find(inside(I) == 1));
    B = I(find(outside(I) == 1));

    if (abs(dist(i)) < 100*eps)
      tooclose = true;
      % also, we'll probably need to invent a normal
      % which is probably a similar procedure to below?
    else
      tooclose = false;
    end

    % nonnormalized normal :-)
    n1 = [xx(i)-cpx(i)  yy(i)-cpy(i)  zz(i)-cpz(i)];

    ips = [];  % inner products
    conf = [];  % confidence
    for jj=1:nnbrs
      j = I(jj);
      n2 = [xx(j)-cpx(j)  yy(j)-cpy(j)  zz(j)-cpz(j)];
      ip = n1 * n2' / (sqrt(n1 * n1') * sqrt(n2 * n2'));
      if inside(j)
        ips = [ips -ip];
        % todo: should scale confidence by dx somehow?
        conf = [conf exp(-abs(d(jj)))];
      elseif outside(j)
        ips = [ips ip];
        conf = [conf exp(-abs(d(jj)))];
      else
        % don't include this one
      end
    end

    if (length(ips) <= 2)
      warning('not enough classified nearby CPs')
      disp('dropping to keyboard');
      keyboard
      % maybe here should be able to defer and come back to this point
    end

    ipmean = mean(ips);
    ipstd = std(ips);

    if ~tooclose && ipmean > 0.8 && ipstd < 0.1
      % clear concensus
      out2(i) = 1;
      text = 'out (clear)  ';
    elseif ~tooclose && ipmean < -0.8 && ipstd < 0.1
      in2(i) = 1;
      text = 'in  (clear)  ';
    elseif tooclose && ipmean > 0.1 && ipstd < 0.1
      out2(i) = 1;
      text = '*close* [out]';
    elseif tooclose && ipmean < -0.1 && ipstd < 0.1
      in2(i) = 1;
      text = '*close* [in] ';
    else
      ips
      conf
      ipmean
      ipstd
      warning('not sure how to classify');
      disp('dropping to keyboard');
      keyboard
    end

    if (verbose)
      fprintf('#%d dist=%8.2g  %s (%.2g, [%d,%d]/%d, %.4g,%.2g)', ...
              i,dist(i),text, max(d),length(A),length(B),nnbrs, ...
              mean(ips),std(ips));
    end
  end

