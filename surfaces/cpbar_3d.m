function [cpx,cpy,cpz, dist, bdy] = cpbar_3d(x,y,z, f, varargin)
%CPBAR_3D  Calculate the cpbar fcn used for boundary conditions.
%   "cpbar" refers to the double closest point projection described
%   in [Macdonald, Brandman, Ruuth, 2011] which is used to
%   implement higher-order accurate boundary conditions.
%
%   "f" is a handle to a cp function which takes x,y,z and whatever is
%   in varargin as inputs.  It must return "cpx,cpy,cpz, dist, bdy".
%   See cpHemisphere.m for example.
%
%   Notes:
%     * dist refers to the distance to the original closest point,
%       not the cpbar point.  This means is safe to use this for
%       banding
%     * bdy is non-zero for points where the original closest point
%       hit a boundary
%
%   Code is vectorized: any size/shape for x should work, provided
%   function handle f is vectorized as well.


[cpx,cpy,cpz, dist, bdy] = f(x,y,z, varargin{:});

% DJS: original implementation.  Changed to match the python code
%cpynew = cpy(Reg);
%cpxnew = cpx(Reg);
%cpznew = cpz(Reg);
%x1 = x(Reg);
%y1 = y(Reg);
%z1 = z(Reg);
%cpx(Reg) = (2*cpxnew - x1);
%cpy(Reg) =(2*cpynew - y1);
%cpz(Reg) = (2*cpznew - z1);
%[cpx,cpy,cpz,Regtemp,dist3Dtemp] = cphemi(cpx,cpy,cpz,R);


% project into interior
x2 = (2*cpx(bdy) - x(bdy));
y2 = (2*cpy(bdy) - y(bdy));
z2 = (2*cpz(bdy) - z(bdy));

% recompute just those ones (it should be ok to do all of them, but
% this should be faster).
[cpx(bdy),cpy(bdy),cpz(bdy), disttemp, bdytemp] = f(x2,y2,z2, varargin{:});
