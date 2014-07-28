function gg = refine_cpgrid_bw(g, bw)
%REFINE_CPGRID_BW  CP grid refinement using bandwidth
%
%   cpgrid_new = refine_cpgrid_bw(cpgrid_old)
%
%   cpgrid_new = refine_cpgrid_bw(cpgrid_old, bw)
%     uses a bandwidth bw*dx to determine the band.  bw defaults to output of
%     rm_bandwidth() if omitted.
%
%   NOTE: there are other possible algorithms.  This one refines
%   computational band and then discards points that are too far away.
%   Perhaps not optimal in terms of how many CPs much be computed.
%   See implementation of refine_cpgrid() for example which uses a
%   sample of surface points and finds a band around it.
%
%   NOTE 2: in high dim/high co-dim, this approach is generally faster
%   (but makes larger bands by virtue of the bandwidth approach).
%
%   TODO: fix or replace the private functions to take cpgrid objects

  if nargin < 2
    bw = rm_bandwidth(g.dim);
  end

  if iscell(g.x1d)
    gg = refine_cpgridnd(g, bw);
  else
    %gg = refine_cpgrid_bw_2d(g, bw);
    if isfield(g,'bdy')
      bdy = g.bdy;
    else
      bdy = [];
    end
    if g.dim == 2
    [band,x,y,cpx,cpy,dist,bdy,dx,x1d,y1d] = ...
        refine_grid2d(g.cpfun,g.dx,g.x1d,g.y1d,bw,g.band, ...
                      g.dist,bdy, false);
    gg.dim = g.dim;
    gg.dx = dx;
    gg.x1d = x1d;
    gg.y1d = y1d;
    gg.cpfun = g.cpfun;
    gg.band = band;
    gg.x = x;
    gg.y = y;
    gg.cpx = cpx;
    gg.cpy = cpy;
    gg.dist = dist;
    gg.bdy = bdy;
    elseif g.dim == 3
    %gg = refine_cpgrid_bw_3d(g, bw);
    [band,x,y,z,cpx,cpy,cpz,dist,bdy,dx,x1d,y1d,z1d] = ...
        refine_grid3d(g.cpfun,g.dx,g.x1d,g.y1d,g.z1d,bw,g.band, ...
                      g.dist,bdy, false);
    gg.dim = g.dim;
    gg.dx = dx;
    gg.x1d = x1d;  gg.y1d = y1d;  gg.z1d = z1d;
    gg.cpfun = g.cpfun;
    gg.band = band;
    gg.x = x;  gg.y = y;  gg.z = z;
    gg.cpx = cpx;  gg.cpy = cpy;  gg.cpz = cpz;
    gg.dist = dist;
    gg.bdy = bdy;
    else
      error('this grid object cannot be handled');
    end
  end

