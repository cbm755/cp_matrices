How to run
==========

1. Start IPython (arguments are for visualization)

    ipython --gui=wx

   If you're using an old IPython version, use this instead

    ipython -wthread

2. Run any example:

    run ex_heat_hemisphere.py

    run ex_heat_sphere.py

    run im_heat_hemisphere.py

    # This crashes for me unless using --gui=qt or running it from a

    # terminal ($ python ex_heat_circle.py)

    run ex_heat_circle.py

The old examples can be ran like this:

    run -i cp/cpm_load

Then execute a particular PDE problem

    run -i cp/cpm_patterns

Also cpm_heat_ex, cpm_heat_im, cpm_eigen
