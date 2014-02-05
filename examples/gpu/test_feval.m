%% Tests for function evaluations f & g in the R-D equations, accelerated on a GPU

cpf = @cpSphere;
paramf = @paramSphere;
%cpf = @cpEllipsoid;
%paramf = @paramEllipsoid;

runcpu = 1
rungpu = 1

loaddata = 1
makeplots = 0;

if (loaddata == 1)
  dx = 0.1/4;      % grid size

  % make vectors of x, y, z positions of the grid
  x1d = (-2.0:dx:2.0)';
  y1d = x1d;
  z1d = x1d;
  nx = length(x1d);
  ny = length(y1d);
  nz = length(z1d);

  % meshgrid is only needed for finding the closest points, not afterwards
  [x y z] = meshgrid(x1d, y1d, z1d);

  [cpx, cpy, cpz, dist] = cpf(x,y,z);
  cpx = cpx(:); cpy = cpy(:); cpz = cpz(:);

  %% Banding: do calculation in a narrow band around the surface
  dim = 3;  % dimension
  p = 3;    % interpolation order
  % "band" is a vector of the indices of the points in the computation
  % band.  The formula for bw is found in [Ruuth & Merriman 2008] and
  % the 1.0001 is a safety factor.
  bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
  band = find(abs(dist) <= bw*dx);

  % store closest points in the band, discarding
  cpx = cpx(band); cpy = cpy(band); cpz = cpz(band);
  x = x(band); y = y(band); z = z(band);
  dist = dist(band);

% parameters and functions for Gray--Scott
FF = 0.054;  kk = 0.063;  nuu = 1/(3/dx)^2;  nuv = nuu/3;
f = @(u,v) (-u.*v.*v  +  FF*(1-u));
g = @(u,v) ( u.*v.*v  -  (FF+kk)*v);

%% initial conditions - small perturbation from steady state
pert = 0.5*exp(-(10*(z-.1)).^2) + 0.5*rand(size(x));
u0 = 1 - pert;  v0 = 0.5*pert;
u = u0;  v = v0;

Tf = 1000;
dt = .2 * (1/max(nuu,nuv)) * dx^2
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps

%% GPU stuff
N = length(u);
blocksz = 512;   % block size for GPU
numblocks = ceil(N / blocksz);

gpu_fg = parallel.gpu.CUDAKernel('kernel_fg.ptx', ...
                               'kernel_fg.cu');
gpu_fg.ThreadBlockSize = [blocksz,1,1];
gpu_fg.GridSize = [numblocks,1];

% upload the arrays
u_d = gpuArray(u);
v_d = gpuArray(v);
unew_d = gpuArray.zeros(size(u_d));
vnew_d = gpuArray.zeros(size(u_d));
t1 = gpuArray.zeros(size(u_d));
rhsu_d = gpuArray.zeros(size(u_d));
rhsv_d = gpuArray.zeros(size(u_d));


if rungpu
tic
for kt = 1:numtimesteps
  
    %unew_d = f(u_d,v_d);
    %vnew_d = g(u_d,v_d);
    [unew_d, vnew_d] = feval(gpu_fg, unew_d, vnew_d, u, v, N);
    t = kt*dt;

end
GPUtime = toc
end

if runcpu
tic
for kt = 1:numtimesteps

unew = f(u,v);
vnew = g(u,v);
t = kt*dt;

CPUtime = toc
end
end

end
factor = CPUtime / GPUtime

