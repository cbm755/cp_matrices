CP_MATRICES

These are matrices used for implementing the Closest Point Method.
See the various example_*.m files, maybe starting with
example_heat_circle.m

These codes assume there is an underlying grid based on meshgrid (or
ndgrid).  But its important to note that one doesn't actually need to
create that grid...

tests/ contains a set of unit tests (feel free to add more!).
run_all_tests.m in the parent directory runs all of them.

../surfaces/ contains various CP representations and tools

TODO: move run_all_tests.m to tests/ and clean up to run only files
that start with the word "test"

TODO: add documentation or clean up an example on how to use these.

