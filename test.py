import numpy as np
import math

A = np.array([[1,2,0,1],[1,0,3,1],[1,0,3,2],[1,2,0,2]])
B = np.array([[3,5,1],[0,8/3,13/3],[0,0,9/8]])



#Householder
print(np.linalg.det(A))