"""
Laplacian Eigenvalues and Eigenfunctions of a surface.  Run cpm_load.py first.
"""
import numpy as np
try:
    from mayavi import mlab
except ImportError:
    from enthought.mayavi import mlab
from time import time
import scipy.linalg
import scipy.sparse.linalg as splinalg
#import scipy.sparse.linalg.eigen.arpack as arpack


if (1==0):
    Ldense,Xdense = scipy.linalg.eig(M.todense())
    Labs = abs(Ldense)
    I = np.argsort(Labs)
    II = I[0:13]
    Ldensef = Ldense[II]
    for i in II:
        eval = Ldense[i]
        evec = Xdense[:,i]
        print eval



print "computing eigenvalues..."
st = time()
#Evals,Evecs = arpack.eigen(-M, k=32, which="SM")
Evals,Evecs = splinalg.eigs(-M, k=32, which="SM")
print "found eigenvalues in time = " + str(time() - st)
# gives very poor results (as expected, our matrix is not symmetric)
#L,X = arpack.eigen_symmetric(M, k=50, which="SM")

#for i,lam in enumerate(L):
    #print i,lam
    #evec = X[:,i]



howManyEvals = len(Evals)
I = np.argsort(Evals)
Evals2 = Evals[I]

Evals_ex = np.array([0, 2,2, 6,6,6, 12,12,12,12, 20,20,20,20,20, 30,30,30,30,30,30,42])
errs = Evals2[0:22].real - Evals_ex
A = "Dx.append(%g); Err.append([" % dx
B = ",".join(map(str, errs))
print A + B + "])"

idx = np.r_[0,1,3,6,10,15,21]
errs2 = Evals2[idx].real - Evals_ex[idx]
A = "Dx.append(%g); Err.append([" % dx
B = ",".join(map(str, errs2))
print A + B + "])"

n = 13



# setup viz
fig = mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1), size=(640,640))
mlab.clf()
# build a pipeline so we can explicitly change the data later
src = mlab.pipeline.grid_source(x, y, z, scalars=z)
normals = mlab.pipeline.poly_data_normals(src)
surf = mlab.pipeline.surface(normals)
# TODO: size is not fontsize
#mlab.title('u: initial conditions', size=0.2)
mlab.show()
pause(1)


#for ii in range(0,n):
#for ii in range(0, howManyEvals):
for ii in range(9, 15+1):
    i = I[ii]
    eval = Evals[i]
    print i,eval
    evec = Evecs[:,i]
    emax = evec.real.max()
    emin = evec.real.min()
    absmax = max(-emin, emax)
    print (-absmax,absmax), (emin,emax)
    # Eplot projects a soln onto a surface grid x,y,z
    #evec_plot = Eplot*E*evec
    src.data.point_data.scalars = Eplot*E*(evec.real)
    src.data.point_data.scalars.name = 'scalars'
    src.data.modified()
    #mlab.title(str(ii) + ' ew=' + str(eval), size=0.2)
    mlab.show()
    mlab.savefig('hemisphere_dx' + str(dx) + '_n' + str(ii) + '_' + str(eval) + '.png')

    if (1==0):
        mlab.clf()
        s = mlab.mesh(x, y, z, scalars=real(evec_plot.reshape(x.shape)), vmin=-absmax, vmax=absmax)
        mlab.title(str(ii) + ' ew=' + str(eval), size=0.2)
        mlab.show()
        mlab.savefig('hemisphere_lev' + str(maxlev) + '_n' + str(ii) + '_' + str(eval) + '.png')


if (1==0):
    Lr = Evals.real
    Ls = np.flipud(np.sort(Lr))
    #Ls.reverse()
    Lf = Ls[0:13]
    Lex = []
    for i in range(1,14):
        Lex.append(-(i)*(i))
    #Lex = [0]
    #for i in range(0,5+1):
    #    Lex.append(-(i+1)*(i+1))
    #    Lex.append(-(i+1)*(i+1))
