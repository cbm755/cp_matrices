"""
Laplacian Eigenvalues and Eigenfunctions of a surface.  Run cpm_load.py first.
"""
from enthought.mayavi import mlab
from time import time
import scipy.linalg
#import scipy.sparse.linalg.eigen
import scipy.sparse.linalg.eigen.arpack as arpack


if (1==0):
    Ldense,Xdense = scipy.linalg.eig(M.todense())
    Labs = abs(Ldense)
    I = numpy.argsort(Labs)
    II = I[0:13]
    Ldensef = Ldense[II]
    for i in II:
        eval = Ldense[i]
        evec = Xdense[:,i]
        print eval



print "computing eigenvalues..."
st = time()
Evals,Evecs = arpack.eigen(M, k=64, which="SM")
print "found eigenvalues in time = " + str(time() - st)
# gives very poor results (as expected, our matrix is not symmetric)
#L,X = arpack.eigen_symmetric(M, k=50, which="SM")

#for i,lam in enumerate(L):
    #print i,lam
    #evec = X[:,i]



howManyEvals = len(Evals)
I = np.argsort(-Evals)

n = 33


f = mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1), size=(640,640))

#for ii in range(0, howManyEvals):
for ii in range(n, n+1):
    i = I[ii]
    eval = Evals[i]
    print i,eval
    evec = Evecs[:,i]
    emax = evec.real.max()
    emin = evec.real.min()
    absmax = max(-emin, emax)
    print (-absmax,absmax), (emin,emax)
    # Eplot projects a soln onto a surface grid x,y,z
    evec_plot = Eplot*E*evec
    mlab.clf()
    s = mlab.mesh(x, y, z, scalars=real(evec_plot.reshape(x.shape)), vmin=-absmax, vmax=absmax)
    mlab.title(str(ii) + ' ew=' + str(eval), size=0.2)
    mlab.show()
    mlab.savefig('hemisphere' + str(ii) + '_' + str(eval) + '.png')


if (1==0):
    Lr = Evals.real
    Ls = numpy.flipud(numpy.sort(Lr))
    #Ls.reverse()
    Lf = Ls[0:13]
    Lex = []
    for i in range(1,14):
        Lex.append(-(i)*(i))
    #Lex = [0]
    #for i in range(0,5+1):
    #    Lex.append(-(i+1)*(i+1))
    #    Lex.append(-(i+1)*(i+1))
