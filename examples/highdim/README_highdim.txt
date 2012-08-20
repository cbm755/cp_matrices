TODO:

modify interpn_matrix to do the banding in the Ei,Ej,Es form then contruct the matrix: that uses less memory.


there is an experimental interpn_matrix_lowmem(): not clear it uses less memory but its certainly slower.  In theory it doesn't construct the Ei,Ej,Es and thus has potential savings in RAM.


possible unit tests:
  highmem gives exactly the same matrix as interp_matrix



           dx     Vector length:
cpgrid2   0.125    139345
cpgrid3   0.0625   455206
cpgrid4   0.03125  1725958

Memory use:

dx       cpgrid structure
0.125      44 MiB
0.0625    166 MiB
0.03125   649 MiB


dx     storage for L   time
0.125     75 MiB       0.6s
0.0625   288 MiB       6.5s
0.03125   machine started swapping (170GiB RAM)

%dx     storage for Ei2
%0.125      885 MiB
%0.0625    3355 MiB
%0.03125

dx     storage for E1   for E2     time
0.125     226 MiB      1668 MiB     7s, 51s 
0.0625    876 MiB      6515 Mib    30s, 225s
0.03125                  26 GiB       , 1100s



dx        cost of E1*u, E2*u
0.125          , 0.5s
0.0625     0.6s, 3.68s
0.03125        , 19.6s
