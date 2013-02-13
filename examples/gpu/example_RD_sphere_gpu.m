%% Reaction-diffusion, accelerated on a GPU

cpf = @cpSphere;
paramf = @paramSphere;
%cpf = @cpEllipsoid;
%paramf = @paramEllipsoid;

runcpu = 0
rungpu = 1

loaddata = 1
makeplots = 0;

if (loaddata == 1)
  dx = 0.2;      % grid size

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


  %% discrete operators
  disp('building laplacian and interp matrices');
  L = laplacian_3d_matrix(x1d,y1d,z1d, 2, band,band);
  [Li,Lj,Ls] = laplacian_3d_matrix_tempcomp(x1d,y1d,z1d,2,band,band);
  % some are zero because their stencil is outside the band
  I = find(Lj == 0);
  Lj(I) = 1;
  Ls(I) = 0;  % these are not necessarily zero b/c of how we build L
  L2 = sparse(Li, Lj, Ls, size(L,1), size(L,2));
  L - L2

  % TODO: could modify laplacian to return the component form
  %[Li,Lj,Ls] = find(L);
  % this won't work b/c of zeros
  %Li = reshape(Li, length(u), 7);
  %Lj = reshape(Lj, length(u), 7);
  %Ls = reshape(Ls, length(u), 7);
  if (1==0)
  Li = repmat((1:length(x))', 1, 2*dim+1);
  Lj = zeros(length(x), 2*dim+1);
  Ls = zeros(length(x), 2*dim+1);
  tic
  for i=1:length(x)
    [I,J,V] = find(L(i,:));
    n = length(I);
    if n < 2*dim + 1
      J = [J ones(1,2*dim+1-n)];
      V = [V zeros(1,2*dim+1-n)];
    end
    Lj(i,:) = J;
    Ls(i,:) = V;
  end
  toc
  L2 = sparse(Li, Lj, Ls, size(L,1), size(L,2));
  L - L2
  end
  E = interp3_matrix(x1d,y1d,z1d, cpx, cpy, cpz, p, band);
  [Ei,Ej,Es] = interp3_matrix(x1d,y1d,z1d, cpx, cpy, cpz, p, band);
  % TODO: a bit silly, interp3_matrix just straightened them out...
  M = (p+1)^dim;
  Ei = reshape(Ei, length(x), M);
  Ej = reshape(Ej, length(x), M);
  Es = reshape(Es, length(x), M);
  E2 = sparse(Ei,Ej,Es, size(E,1), size(E,2));
  E - E2
  %[Ei,Ej,Es] = interp3_matrix(x1d,y1d,z1d, cpx, cpy, cpz, p);
  I = speye(size(E));

  %% plotting grid
  [xp,yp,zp] = paramf(128);
  % Eplot is a matrix which interpolations data onto the plotting grid
  Eplot = interp3_matrix(x1d, y1d, z1d, xp(:), yp(:), zp(:), p, band);
end

% u_t = f(u,g) + nuu*Lap u
% v_t = g(u,g) + nuv*Lap u

% parameters and functions for Gray--Scott
FF = 0.054;  kk = 0.063;  nuu = 1/(3/dx)^2;  nuv = nuu/3;
f = @(u,v) (-u.*v.*v  +  FF*(1-u));
g = @(u,v) ( u.*v.*v  -  (FF+kk)*v);


%% initial conditions - small perturbation from steady state
pert = 0.5*exp(-(10*(z-.1)).^2) + 0.5*rand(size(x));
u0 = 1 - pert;  v0 = 0.5*pert;
u = u0;  v = v0;

Tf = 200;
dt = .2 * (1/max(nuu,nuv)) * dx^2
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps

if makeplots
  figure(1); clf;
  sphplot = Eplot*u;
  sphplot = reshape(sphplot, size(xp));
  Hplot = surf(xp, yp, zp, sphplot);
  xlabel('x'); ylabel('y'); zlabel('z');
  axis equal
  view(-10, 60)
  axis off;
  shading interp
  camlight left
  colorbar
end

%% Method-of-lines approach
% See [vonGlehn/Macdonald/Maerz 2013]
%lambda = 8*max(nuu,nuv)/(dx^2);
%Au = nuu*(E*L) - lambda*(I-E);
%Av = nuv*(E*L) - lambda*(I-E);

%% GPU stuff
N = length(u);
blocksz = 512;   % block size for GPU
numblocks = ceil(N / blocksz);

spMV = parallel.gpu.CUDAKernel('kernel_spmatvec2.ptx', ...
                               'kernel_spmatvec2.cu');
spMV.ThreadBlockSize = [blocksz,1,1];
spMV.GridSize = [numblocks,1];

% convert the indices to int32 (8 bytes, same a cuda's int)
Ei = int32(Ei);
Ej = int32(Ej);
Li = int32(Li);
Lj = int32(Lj);

% upload the arrays
u_d = gpuArray(u);
v_d = gpuArray(v);
Ej_d = gpuArray(Ej);
Es_d = gpuArray(Es);
Lj_d = gpuArray(Lj);
Ls_d = gpuArray(Ls);
unew_d = gpuArray.zeros(size(u_d));
vnew_d = gpuArray.zeros(size(u_d));
t1 = gpuArray.zeros(size(u_d));
rhsu_d = gpuArray.zeros(size(u_d));
rhsv_d = gpuArray.zeros(size(u_d));

starttime = cputime();
for kt = 1:numtimesteps
  %% MOL: explicit Euler timestepping
  %unew = u + dt*( E*f(u,v) + Au*u );
  %vnew = v + dt*( E*g(u,v) + Av*v );
  %u = unew;
  %v = vnew;
  %% MOL: without precomputing matrices
  %rhsu = nuu*(L*u) + f(u,v);
  %rhsv = nuv*(L*v) + g(u,v);
  %unew = u + dt*( E*rhsu - lambda*(u - E*u) );
  %vnew = v + dt*( E*rhsv - lambda*(v - E*v) );
  %u = unew;
  %v = vnew;

  %% Ruuth-Merriman

  % TODO: have one call do multiple RHS?

  if rungpu
    tic
    t1 = feval(spMV, t1, Ls_d, Lj_d, u_d, N, 2*dim+1);
    % TODO: there was a command to run the function on each element
    % instead---check if that is faster
    rhsu_d = f(u_d,v_d) + nuu*t1;
    t1 = feval(spMV, t1, Ls_d, Lj_d, v_d, N, 2*dim+1);
    rhsv_d = g(u_d,v_d) + nuv*t1;
    unew_d = u_d + dt*rhsu_d;
    vnew_d = v_d + dt*rhsv_d;
    u_d = feval(spMV, u_d, Es_d, Ej_d, unew_d, N, M);
    v_d = feval(spMV, v_d, Es_d, Ej_d, vnew_d, N, M);
    GPUtime = toc;
  end

  if runcpu
    tic
    rhsu = nuu*(L*u) + f(u,v);
    rhsv = nuv*(L*v) + g(u,v);
    unew = u + dt*rhsu;
    vnew = v + dt*rhsv;
    u = E*unew;
    v = E*vnew;
    CPUtime = toc;
  end

  if runcpu && rungpu
    %[kt norm(u-u_d)  norm(v-v_d) GPUtime CPUtime]
  end

  t = kt*dt;

  if makeplots && ((mod(kt,20)==0) || (kt<=10) || (kt==numtimesteps))
    disp([kt t]);
    sphplot = Eplot*u;
    sphplot = reshape(sphplot, size(xp));
    set(0, 'CurrentFigure', 1);
    set(Hplot, 'CData', sphplot);
    title( ['u at time ' num2str(t) ', kt= ' num2str(kt)] );
    drawnow;
  end
end
tstime = cputime() - starttime

