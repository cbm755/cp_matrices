function helper_tri2cp(dx, relpt, bw, F, V, gridhashsz, debuglvl)
%HELPER_TRI2CP  Helper function for TRI2CP
%   Don't call this directly, see TRI2CP instead

warning('Don''t call this directly, see TRI2CP instead');

%mex helper_tri2cp.c
mex -O CFLAGS='\$CFLAGS -Wall -std\=c99' helper_tri2cp.c
