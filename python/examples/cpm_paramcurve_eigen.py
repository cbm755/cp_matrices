"""
Laplacian Eigenvalues and Eigenfunctions of a surface.  Run cpm_load.py first.
"""
#from enthought.mayavi import mlab
from time import time
import scipy.linalg
#import scipy.sparse.linalg.eigen
import scipy.sparse.linalg.eigen.arpack as arpack


if (1==0):
    Ldense,Xdense = scipy.linalg.eig(M.todense())
    #Labs = abs(Ldense)
    I = numpy.argsort(Ldense)
    II = I[0:13]
    Ldensef = Ldense[II]
    for i in II:
        eval = Ldense[i]
        evec = Xdense[:,i]
        print eval



print "computing eigenvalues..."
st = time()
fac = (1/dx)**2
#Evals,Evecs = arpack.eigen(-M, k=16, which="SM", tol=1e-15, maxiter=200000, ncv=128)
Evals,Evecs = arpack.eigen(-M, k=16, which="SM", maxiter=200000, ncv=128)
print "found eigenvalues in time = " + str(time() - st)
# gives very poor results (as expected, our matrix is not symmetric)
#L,X = arpack.eigen_symmetric(M, k=50, which="SM")

#for i,lam in enumerate(L):
    #print i,lam
    #evec = X[:,i]


pylab.figure(1)
pylab.clf()
pylab.plot(real(Evals),imag(Evals),'rx')
pylab.title('spectrum')


howManyEvals = len(Evals)
I = np.argsort(Evals)
Evals2 = Evals[I]

n = 33


#f = mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1), size=(640,640))

figure(2)
clf()

for ii in range(0, 5):
#for ii in range(0, howManyEvals):
#for ii in range(n, n+1):
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
    plot(x, evec_plot, '-', color=numpy.random.rand(4,))
    #mlab.clf()
    #s = mlab.mesh(x, y, z, scalars=real(evec_plot.reshape(x.shape)), vmin=-absmax, vmax=absmax)
    #mlab.title(str(ii) + ' ew=' + str(eval), size=0.2)
    #mlab.show()
    #mlab.savefig('hemisphere' + str(ii) + '_' + str(eval) + '.png')

# using maple to evalute the length of the curve cos(x) from x=0..3*Pi/2
#length = 5.73029668354156802685714312327
# x = 1/4..4
length = 4.510147677480894444218572846480212229452198844910047382329043518
# circular arc, pi/6 to 3*pi/2
#length = pi

if (1==1):
    # cos, dirichlet
    Evals_ex = arange(1,10+1)**2 * ((0.5*2*pi/length)**2)
    print dx, Evals_ex - Evals2[0:10].real

    # cos, neumann
    Evals_ex = arange(0,9+1)**2 * ((0.5*2*pi/length)**2)
    print dx, Evals_ex - Evals2[0:10].real
    #print dx, dx**2*Evals_ex - Evals2[0:10].real

    # circular arc, pi/6 to 3*pi/2
    Evals_ex = arange(1,10+1)**2 * ((2*pi/length)**2)
    #Evals_ex = arange(1,5+1)**2
    print dx, Evals_ex - Evals2[0:10].real
    #print dx, dx**2*Evals_ex - Evals2[0:10].real

    # circle
    Evals_ex = arange(0,9)**2
    print dx, Evals_ex - Evals2[ix_([0,1,3,5,7,9,11,13,15])].real

    # semicircle
    Evals_ex = arange(1,10)**2
    print dx, Evals_ex - Evals2[0:9].real

    # egg-shape (just use circle, its same)
    #ex = a([0,1,1,2,2,3,3,4,4,5,5,6,6])**2
    #print dx, ex - Evals2[0:13].real
show()
