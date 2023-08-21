import numpy as np
import math

A = np.asarray([[1,0,0],[0,2,0],[0,0,3]]);
B= np.asarray([[1,2,4],[2,13,23],[4,23,77]]);


#Householder
print(np.dot(A, B))