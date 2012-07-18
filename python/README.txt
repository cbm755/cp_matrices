# How to run this:
# start the python interpreter (arguments are for vizualization)
ipython -wthread -pylab
# In recent IPython versions, -wthread doesn't exist. Use --gui=wx
ipython --gui=wx
# then within python, type the following (-i is important)
run -i cpm_load
# then execute a particular PDE problem
run -i cpm_patterns
# (also cpm_heat_ex, cpm_heat_im, cpm_eigen)

