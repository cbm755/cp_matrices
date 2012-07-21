function mu = logarithmicInfNorm(A)
%LOGARITHMICINFNORM   Compute the logarithmic norm (inf norm) of a matrix
%

  Ad = diag(A);

  A_offdiag = A - diag(diag(A));

  s = real(Ad) + sum(abs(A_offdiag),2);
  %minmax(s)
  %s = real(Ad) + sum((A_offdiag),2)
  mu = max(full(s));

  %keyboard