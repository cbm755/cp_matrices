import numpy as np

from cp import surfaces

sur = surfaces.sur(interp_degree=3, levels=5)

sur.grid()
sur.find_stencils()

D = sur.buildDiffMatrix()
E = sur.buildExtensionMatrix()
EPlot = sur.buildEPlotMatrix()

u0 = np.random.rand(D.shape[0])
