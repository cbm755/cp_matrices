%% Test of my sparse matrix-vector multiply on GPU
% This uses Matlab's Parallel Computing Toolkit.

% TODO: needs to be inside a function to avoid a memcpy off of GPU?

N = 100000;     % vector length
M = 128;        % nnz per row
blocksz = 512;  % block size for GPU

numblocks = ceil(N / 10);

kern = parallel.gpu.CUDAKernel('kernel_spmatvec.ptx', ...
                               'kernel_spmatvec.cu');

kern.ThreadBlockSize = [blocksz,1,1];
kern.GridSize = [numblocks,1];

%% build some random matrices
% (easy to do this on the GPU really)
x = rand(N,1);

Ai = (1:N)';
Ai = repmat(Ai, 1, M);
Aj = ceil(rand(N,M)*N);
if M == 1
  As = (1:N)'*100
elseif M == 2
  As = [(1:N)'*100 (1:N)'*1000];
else
  As = rand(N,M);
end
A = sparse(Ai,Aj,As,N,N);

% CPU approach with Matlab's sparse matrix
tic
y2 = A*x;
cputime = toc

% CPU loop approach
if (1==0)
  y = zeros(size(x));
  for i=1:N
    y(i) = As(i,:) * x(Aj(i,:));
  end
  norm(y-y2)
end


%% upload to GPU
tic
xd = gpuArray(x);
Ajd = gpuArray(Aj);
Asd = gpuArray(As);
ultime = toc



y3 = gpuArray(zeros(size(x)));
tic
y3 = feval(kern, y3, Asd, Ajd, xd, N, M);
GPUtime1 = toc
tic
y3 = feval(kern, y3, Asd, Ajd, xd, N, M);
GPUtime2 = toc

%[y2 y3 y2-y3 Aj As x]
disp('should be small:')
norm(y2-y3)
disp(sprintf('N=%d, block=%d, blocksz=%d',N,numblocks,blocksz));
factor = cputime / GPUtime2

