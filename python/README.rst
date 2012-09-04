Installation
============

You need NumPy >= 1.6 and scipy. You need a C compiler to build an
extension and Cython >=0.16. Matplotlib (>=1) and
mayavi are used for visualization, and IPython to interact with the
examples.

The easiest way to get everything is to get EPD, or Ubuntu 12.04.

Checkout the code from the repo, and run::

    $ python setup.py build_ext --inplace

to build the Cython extension, and you're good to go.

How to run
==========

1. Start IPython (arguments are for visualization)::

    ipython --gui=wx

   If you're using an old IPython version, use this instead::

    ipython -wthread

2. Run the example using the new functions::

    run examples/example_heat_triangulated_sphere.py

   Or any of the cleaned-up examples::

    run examples/ex_heat_hemisphere.py
    run examples/ex_heat_sphere.py
    run examples/im_heat_hemisphere.py
    # This crashes for me unless using --gui=qt or running it from a
    # terminal ($ python ex_heat_circle.py)
    run examples/ex_heat_circle.py

The old examples can be ran like this::

    run -i examples/cpm_load

Then execute a particular PDE problem::

    run -i examples/cpm_patterns

Also: cpm_heat_ex, cpm_heat_im, cpm_eigen
