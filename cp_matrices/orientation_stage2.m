function [inside,outside,unknown] = orientation_stage2(xx,yy,zz,cpx,cpy,cpz,dist,dx,E,inside0,outside0,verbose)
%ORIENTATION_STAGE2  Helper function for ORIENTATION_FROM_CP
%   Classify a few closest points as inside/outside supposing that
%   most are already known.
%
%   Uses a local calculation based on interpolating distances and
%   neighbouring CPs.  You can probably invent (fine-scale) geometry
%   that will fool these.
%
%   WARNING: not well tested, will probably need other cases added
%   to the classifier.

  %% Parameters
  nnbrs = 6;   % how many neighbours near the CP
  DROP_TO_KEYBOARD_OK = 1;
  % important that this stays same as in normals_from_cp.m
  TOO_CLOSE_TOL = 100*eps;

  % many other tweakable parameters related to classifying

  %xx = xx(:);
  %yy = yy(:);
  %zz = zz(:);
  %cpx = cpx(:);
  %cpy = cpy(:);
  %cpz = cpz(:);
  %dist = dist(:);

  inside = inside0;
  outside = outside0;

  nnc0 = nnz(~inside0 & ~outside0);

  unknown = zeros(size(inside));

  n1 = xx - cpx;
  n2 = yy - cpy;
  n3 = zz - cpz;
  len = sqrt(n1.*n1 + n2.*n2 + n3.*n3);
  % normalize the normals
  n1 = n1 ./ len;
  n2 = n2 ./ len;
  n3 = n3 ./ len;

  skip = zeros(size(xx));

  pass = 0;
  num_class = 42;  % dummy, nonzero value
  while(1)
    pass = pass + 1;
    if num_class == 0
      force = 1;
      if verbose >= 2
        disp('stage2: no progress last pass, forcing first');
      end
    else
      force = 0;
    end

    num_class = 0;  % reset's each pass


    noclass = find(~inside & ~outside & ~unknown);
    nnc = length(noclass);
    if (nnc == 0)
      break;
    end

    if verbose >= 1
      fprintf('stage2: starting pass %d, %d pts to do\n', pass, nnc);
    end

    for ii=1:nnc
      i = noclass(ii);

      if (abs(dist(i)) < TOO_CLOSE_TOL)
        tooclose = true;
      else
        tooclose = false;
      end

      [io1,skip1,data1] = classify_dist(...
          i, abs(dist), inside, outside, E, force);

      [io2,io3,skip2,skip3,data2] = classify_nbrs(...
          i,inside,outside,dx,nnbrs,xx,yy,zz,cpx,cpy,cpz,n1,n2,n3,force);

      io = [io1 io2 io3];

      drop_key = 0;
      if all(io == [1 1 1])
        outside(i) = 1;
        text = 'outside';
        num_class = num_class + 1;
      elseif all(io == [-1 -1 -1])
        inside(i) = 1;
        text = 'inside';
        num_class = num_class + 1;
      elseif pass>1 && all(io == [1 0 1]) || all(io == [1 1 0])
        outside(i) = 1;
        text = 'outside';
        num_class = num_class + 1;
      elseif pass>1 && all(io == [-1 0 -1]) || all(io == [-1 -1 0])
        inside(i) = 1;
        text = 'inside';
        num_class = num_class + 1;
      % I thought I needed these cases but I had a bug instead
      %elseif pass>1 && io2==1 && io3==1 ...
      %      && mean(data2.ips)>0.9 && data2.weightsum > 0.8
      %  outside(i) = 1;
      %  text = 'outside';
      %  num_class = num_class + 1;
      %elseif pass>1 && io2==-1 && io3==-1 ...
      %      && mean(data2.ips)<-0.9 && data2.weightsum < -0.8
      %  inside(i) = 1;
      %  text = 'inside';
      %  num_class = num_class + 1;
      elseif tooclose
        % probably doesn't really matter which we pick here
        if io2 == 1
          outside(i) = 1;
          text = 'close, outside';
        else  %if io2 == -1;
          inside(i) = 1;
          text = 'close, inside';
        end
        num_class = num_class + 1;
      else
        if ~force && (skip1 == 1 || skip2 == 1 || skip3 == 1)
          text = '*skipping*';
        else
          text = '*unknown*';
          if DROP_TO_KEYBOARD_OK
            drop_key = 1;
          else
            unknown(i) = 1;
          end
        end
      end

      if (verbose >= 2)
        if isempty(data1)
          L1 = [];
          L2 = [];
        else
          L1 = abs(mean(data1.cpval1 ./ data1.cpval2));
          L2 = data1.nncII;
        end
        fprintf('  p#%d:%d(%d%%) dist=%8.2g  %s\t%2d%2d%2d,%d%d%d %.0e[%d],[%d%d],%.3g\n', ...
                pass,i,round(i/nnc),dist(i),text,io1,io2,io3,skip1,skip2,skip3,...
                L1,L2,data2.numneg,data2.numpos,data2.weightsum);
      end

      if drop_key %&& DROP_TO_KEYBOARD_OK
        disp('dropping to keyboard');
        i
        dist(i)
        tooclose
        %weightsum
        data1
        data2
        %plot3(xx(I), yy(I), zz(I), 'k.')
        %plot3(cpx(II), cpy(II), cpz(II), 'ko')
        hold on;
        plot3([xx(i) cpx(i)], [yy(i) cpy(i)], [zz(i) cpz(i)], 'k-')
        plot3(xx(i), yy(i), zz(i), 'k*')
        plot3(cpx(i), cpy(i), cpz(i), 'ko')
        II = data2.II;
        for kk = 1:length(II)
          k = II(kk);
          H1 = plot3([xx(k) cpx(k)],[yy(k) cpy(k)],[zz(k) cpz(k)],'r-');
          H2 = plot3(cpx(k), cpy(k), cpz(k), 'ro');
          if inside(k)
            set(H1, 'color', [0.7, 0.7, 0.7])
            set(H2, 'color', [0.7, 0.7, 0.7])
          end
        end
        fprintf('To manually select, type either:\n');
        fprintf('  inside(i) = 1; return\n');
        fprintf('  outside(i) = 1; return\n');
        fprintf('  unknown(i) = 1; return\n');
        keyboard
      end

    end % for
  end % passes

  unknown = find(unknown);

  if ~isempty(unknown)
    fprintf('Orientation: could not classify %d (of %d non-fill-classified points)\n', ...
            length(unknown), nnc0);
    if (verbose >= 1)
      unknown
      disp('distances/dx');
      (dist(unknown)/dx)'
    end
  end

end




function [io,skip,data] = classify_dist(i,dist,inside,outside,E,force)
%

  Ei = E(i,:);
  II = find(Ei);
  ncII = ~inside & ~outside;
  ncII = II(ncII(II));
  % remove the current point
  ncII = ncII(ncII ~= i);
  %[nnz(inside(II))  nnz(outside(II))  length(ncII)]

  nncII = length(ncII);

  skip = 0;
  % If too many stencil points are unknown, then it seems this approach
  % would not be that reliable (in practice it seems fine).  It also
  % gets very slow for large values, so we use random patterns below.
  % (empirical results: 5: 0.04s, 8: 0.3s, 12: 5s).
  if nncII >= 8
    skip = 1;
    %if ~force
    %  io = 0;
    %  data = [];
    %  return
    %end
  end

  %sdist = -1*inside(II).*dist(II) + outside(II).*dist(II);
  sdist = -1*inside.*abs(dist) + outside.*abs(dist);
  sdist1 = sdist; sdist1(i) = -abs(dist(i));
  sdist2 = sdist; sdist2(i) = abs(dist(i));
  if nncII >= 8
    pat = make_random_binary_patterns(nncII,256);
  else
    pat = make_binary_patterns(nncII);
  end
  cpval1 = zeros(1,size(pat,1));
  cpval2 = zeros(1,size(pat,1));
  for s=1:size(pat,1)
    for k=1:nncII
      iii = ncII(k);
      sdist1(iii) = (-1)^pat(s,k)*abs(dist(iii));
      sdist2(iii) = (-1)^pat(s,k)*abs(dist(iii));
    end
    cpval1(s) = Ei*sdist1(:);
    cpval2(s) = Ei*sdist2(:);
  end
  data.nncII = nncII;
  data.pat = pat;
  data.cpval1 = cpval1;
  data.cpval2 = cpval2;

  if all(abs(cpval1) < 1/3*abs(cpval2))
    io = -1;
  elseif all(abs(cpval1) > 3*abs(cpval2))
    io = 1;
  elseif all(abs(cpval1) < abs(cpval2))
    io = -1;
    skip = 1;
  elseif all(abs(cpval1) > abs(cpval2))
    io = 1;
    skip = 1;
  else
    io = 0;
    skip = 1;
  end
end




function [io1,io2,skip1,skip2,data] = classify_nbrs(...
    i,inside,outside,dx,nnbrs,xx,yy,zz,cpx,cpy,cpz,n1,n2,n3,force)
%

  dim = 3;
  skip1 = 0;
  skip2 = 0;

  tcp = [cpx(i) cpy(i) cpz(i)];

  % TODO: could accelerate this search with a kdtree
  dd = (cpx(:)-tcp(1)).^2 + (cpy(:)-tcp(2)).^2 + (cpz(:)-tcp(3)).^2;
  dd(i) = 1e42;  % big, so we don't find the current point
  [dd,II] = sort(dd, 'ascend');

  wh = inside(II) | outside(II);
  wh2 = find(wh, nnbrs);

  % we allow at most two unknown neighbour
  if max(wh2) >= nnbrs+2
    %disp('maybe skip');
    skip1 = 1; skip2 = 1;
  end

  II = II(wh2);
  dd = dd(wh2)';

  if any(dd > 2*dx*dx)
    % want them close but not sure how much it matters
    % TODO: skip here?
    warning('nbrs not close enough?');
    disp('dropping to keyboard');
    keyboard
  end
  A = II(find(inside(II) == 1));
  B = II(find(outside(II) == 1));

  % Weights fall off with 1/distance^2.
  dd2 = dd;
  dd2(dd2 == 0) = eps;  % just to avoid inf
  w = 1./dd2;
  weights = w / sum(w);

  % nonnormalized normal, TODO: what if this is exactly zero? (a
  % too-close case, who cares?)
  nt = [xx(i)-cpx(i)  yy(i)-cpy(i)  zz(i)-cpz(i)];
  ntl = norm(nt,2);

  ns = zeros(nnbrs, dim);
  ns(:,1) = n1(II);
  ns(:,2) = n2(II);
  ns(:,3) = n3(II);

  ips = nt(1)*n1(II) + nt(2)*n2(II) + nt(3)*n3(II);
  % TODO: ips will be inf if ntl is exactly 0
  % it really doesn't matter how we classify in that case
  ips = nt * ns' / ntl;

  % now flip the sign of inside normals
  ips = -1*inside(II)'.*ips  +  outside(II)'.*ips;

  %ips2 = sort(ips,'descend');
  %ips2 = ips2(1:end-1);
  %ipmean = mean(ips);
  %ipstd = std(ips);

  numneg = sum(sign(ips) == -1);
  numpos = sum(sign(ips) == 1);


  if numneg >= nnbrs-1
    io1 = -1;
  elseif numpos >= nnbrs-1
    io1 = 1;
  else
    io1 = 0;
    skip1 = 1;
  end

  weightsum = weights * ips';
  if weightsum > 0.25
    io2 = 1;
  elseif weightsum < -0.25
    io2 = -1;
  else
    io2 = 0;
    skip2 = 1;
  end

  data.ips = ips;
  data.weights = weights;
  data.numneg = numneg;
  data.numpos = numpos;
  data.weightsum = weightsum;
  data.II = II;
end



function pat = make_binary_patterns(N)
% Given 3, gives [0 0 0; 0 0 1; 0 1 0; etc]
% a bit ugly but gets the job done
  if N == 0
    pat = [1];
    return
  end
  pat = zeros(2^N,N);

  NN = 2*ones(1,N);
  for s=1:2^N
    [iii{1:N}] = ind2sub(NN, s);
    for k=1:N
      pat(s,k) = iii{k} - 1;
    end
  end

end

function pat = make_random_binary_patterns(N,M)
% some (close to but less than M) patterns of length N with random 1
% entries.
  pat = zeros(M,N);

  for s=1:M
    A = rand(1,N);
    A = A > rand;
    pat(s,:) = A;
  end
  pat = unique(pat,'rows');
end
