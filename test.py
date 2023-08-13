import numpy as np

A = np.array([[1,0,0],[2/3,1,0],[1/3,1/8,1]])
B = np.array([[3,5,1],[0,8/3,13/3],[0,0,9/8]])

print(np.dot(A,B))