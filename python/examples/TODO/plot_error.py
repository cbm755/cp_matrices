'''
Created on Aug 15, 2012

@author: nullas
'''
from matplotlib import pylab as pl
import pickle

if __name__ == '__main__':
    ferror = open('error.pickle','r')
    fdx = open('dx.pickle','r')
    error = pickle.load( ferror)
    dx = pickle.load(fdx)
    pl.loglog(dx,error)
    pl.xlabel('dx')
    pl.ylabel('max error')
    pl.title('error')
    pl.show()