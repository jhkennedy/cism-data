import numpy as np
import scipy
from scipy.spatial import cKDTree

x = [1,2,3]
y = [1,2,3]

xx, yy =  scipy.meshgrid(x, y, indexing='ij')

xy = np.ndarray( (len(xx.ravel() ), 2) )

xy[:,0] = xx.ravel()
xy[:,1] = yy.ravel()

print(xy)

tree = cKDTree( xy )

qp = (4., 4.)
d, indx = tree.query( qp, k=4) 

print('x neighbors: ')
print(xy[indx,0])

print('y neighbors: ')
print(xy[indx,1])

print('Distances: ')
print(d)
