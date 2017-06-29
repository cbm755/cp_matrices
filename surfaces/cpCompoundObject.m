function [cpx, cpy, dist] = cpCompoundObject(x,y, cpfs)
%CPCOMPOUNDOBJECT  Closest Point function for a combination of objects
%   [cpx,cpy,dist] = cpCompoundObject(x, y, cp_fcn_handles)
%      where 'cp_fcn_handles' is a cell array of function handles to
%      cp functions, each of which takes exactly two arguments (x,y).
%
%   TODO: currently works for closed objects.  Open objects will
%   not have their boundaries identified correctly.
%
%   TODO: this could probably be written in a dimensional
%   independent form.


  [cpx, cpy, dist, bdy] = cpfs{1}(x, y);
  %figure(2); clf;
  %porcupine_plot2d(x,y,cpx,cpy,bdy,2);
  %figure(3); clf; pcolor(x,y,dist); axis equal
  %pause(.1);
  for i=2:length(cpfs)
    [cpx2, cpy2, dist2, bdy2] = cpfs{i}(x, y);
    %figure(2); clf;
    %porcupine_plot2d(x,y,cpx2,cpy2,bdy2,2);
    %figure(3); clf; pcolor(x,y,dist2); axis equal
    %pause(.1);
    wh = (dist2 < dist);
    dist = (wh) .* dist2  +  (~wh) .* dist;
    cpx  = (wh) .* cpx2   +  (~wh) .* cpx;
    cpy  = (wh) .* cpy2   +  (~wh) .* cpy;
    % TODO: do something with bdy for open?
  end
