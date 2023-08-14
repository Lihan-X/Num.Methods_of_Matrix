import numpy as np
import math

A = np.array([[1,2,0,1],[1,0,3,1],[1,0,3,2],[1,2,0,2]])
B = np.array([[3,5,1],[0,8/3,13/3],[0,0,9/8]])



#Householder
u1 = 0.5*np.array([[-1],[1],[1],[1]])
H1=np.identity(4)-2*np.dot(u1,np.transpose(u1))
print(H1)
A1=np.dot(H1,A)
print(A1)
A1=np.array([[ 0,  0, -1.],[ 0,  0,  0.],[ 2, -3,  0]])
u2=math.sqrt(2)*0.5*np.array([[-1],[0],[1]])
H2=np.identity(3)-2*np.dot(u2,np.transpose(u2))
print(H2)
A2=np.dot(H2,A1)
print(A2)