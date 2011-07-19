function plot2d_compdomain(u,x,y,dx,dy,fignum)
%PLOT2D_COMPDOMAIN  Plot a 2D Closest Point computational domain
%   plot2d_compdomain(u,x,y,dx,dy)
%   plot2d_compdomain(u,x,y,dx,dy,cpx,cpy,fignum)

  if (nargin < 6) | isempty(fignum)
    fignum = figure()
  end

  figure(fignum); clf;
  % plot a patch for each grid point
  xpat = dx/2*[-1; 1; 1; -1];
  ypat = dy/2*[1; 1; -1; -1];
  X = repmat(x',4,1) + repmat(xpat,1,length(u));
  Y = repmat(y',4,1) + repmat(ypat,1,length(u));
  patch(X,Y,u');
  % TODO: could add the outerband as well
  %ufull = 0/0 * zeros(ny,nx);
  %ufull(innerband) = u;
  %pcolor(x1d,y1d,ufull);
  shading flat
  %caxis([-1.05 1.05]);
  hold on
  %plot(xp,yp,'k-', 'linewidth', 2);
  axis equal;  axis tight;
  colorbar;
  %title( ['embedded domain: soln at time ' num2str(t) ...
  %        ', timestep #' num2str(kt)] );
  xlabel('x'); ylabel('y');
