import numpy as np
import math

A = np.array([[1,2,0,1],[1,0,3,1],[1,0,3,2],[1,2,0,2]])
B = np.array([[1,0,0],[3,1,0],[2,-2,1]])



#Householder
print(np.linalg.inv(B))