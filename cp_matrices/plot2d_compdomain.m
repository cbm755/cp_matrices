function plot2d_compdomain(u,x,y,dx,dy,fignum)
%PLOT2D_COMPDOMAIN  Plot a 2D Closest Point computational domain
%   plot2d_compdomain(u,x,y,dx,dy)
%   plot2d_compdomain(u,x,y,dx,dy,fignum)
%   You must have already created the figure with "figure(fignum)"
%   for the later version.
%
%   plot2d_compdomain(u, g);
%   plot2d_compdomain(u, g, fignum);
%   Alternatively, pass a cpgrid object 'g'.
%
%   TODO: could add the outerband as well
%   TODO: allow plotting a parameterization as well

  if nargin == 2 || nargin == 3
    g = x;
    if nargin == 3
      fignum = y;
    else
      fignum = [];
    end
    x = g.x;
    y = g.y;
    dx = g.dx;
    dy = g.dx;  % TODO: hardcoped for dx=dy
  end

  if nargin == 5
    fignum = []
  end

  if isempty(fignum)
    fignum = figure();
  end

  set(0, 'CurrentFigure', fignum);
  clf;
  % plot a patch for each grid point
  xpat = dx/2*[-1; 1; 1; -1];
  ypat = dy/2*[1; 1; -1; -1];
  X = repmat(x',4,1) + repmat(xpat,1,length(u));
  Y = repmat(y',4,1) + repmat(ypat,1,length(u));

  % this works on matlab
  %H = patch(X,Y,u');
  % some bugs in octave 3.4.2? anyway, below works on both:
  H = patch(X,Y,'g');
  set(H,'FaceColor','flat','FaceVertexCData',u)

  hold on
  %shading flat
  %caxis([-1.05 1.05]);
  % plot parameterization of the surface
  %plot(xp,yp,'k-', 'linewidth', 2);
  axis equal; axis tight;
  colorbar;
  xlabel('x'); ylabel('y');
