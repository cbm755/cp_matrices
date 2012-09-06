TODO:

modify interpn_matrix to do the banding in the Ei,Ej,Es form then contruct the matrix: that uses less memory.


there is an experimental interpn_matrix_lowmem(): not clear it uses less memory but its certainly slower.  In theory it doesn't construct the Ei,Ej,Es and thus has potential savings in RAM.


possible unit tests:
  highmem gives exactly the same matrix as interp_matrix


Results in 4D
-------------
dx=0.2, lambda=2*dim/dx^2, err=0.014045
dx=0.1, lambda=2*dim/dx^2, err=0.0068736
dx=0.05, lambda=2*dim/dx^2, err=0.003416
dx=0.025, lambda=2*dim/dx^2, err=0.0017013


dx=0.2, lambda=-2*dim/dx^2, err=0.11624
dx=0.1, lambda=-2*dim/dx^2, err=0.0068148
dx=0.05, lambda=-2*dim/dx^2, err=0.0030842
dx=0.025, lambda=-2*dim/dx^2, err=0.0015363

lambda=-4/dx^2, 0.0012, unstable, 0.0010566, <running>

Results in 5D
-------------
timings, etc on cyclops.maths

           dx     Vector length   time
cpgrid2   0.125       455206
cpgrid3   0.0625     1725958
cpgrid4   0.03125    6793374      0.4s
cpgrid5   0.015625  26993730      1.5s

Memory use:
dx       cpgrid structure
0.125       44 MiB
0.0625     166 MiB
0.03125    649 MiB
0.015625  2591 MiB

dx      storage for L   time
0.125        76 MiB      2.7s
0.0625      282 MiB      12s
0.03125    1136 MiB      55s
0.015625   4512 MiB     245s

%dx     storage for Ei2
%0.125      885 MiB
%0.0625    3355 MiB

dx     storage for E1   for E2     time
0.125     226 MiB      1668 MiB     7s, 51s
0.0625    876 MiB      6515 Mib    30s, 225s
0.03125  3472 MiB     25910 MiB   146s, 1100s

dx        cost of L*u, E1*u, E2*u
0.125      0.02s,  0.07s, 0.5s
0.0625     0.089s, 0.6s, 3.68s
0.03125    0.038s, 3.6s, 19.6s


