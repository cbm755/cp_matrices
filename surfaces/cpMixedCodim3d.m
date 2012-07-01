function [cpx cpy cpz dist bdy ls] = cpMixedCodim3d(x,y,z)
%CPMIXEDCODIM3D  A mixed codimension test case

  [cpx cpy cpz dist bdy ls] = ...
      cpSurfOfRevolution(x,y,z, @cpMixedCodim2d, 'y', ...
                         {1 0.95 0.6});

