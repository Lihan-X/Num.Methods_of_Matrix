import numpy as np
import math

A = np.asarray([[4,-1,1],[9,-8,9],[11,-11,12]]);
B= np.asarray([[1,2,4],[2,13,23],[4,23,77]]);

A = A*np.identity(3)
#Householder
Q, R = np.linalg.qr(A)
print(Q, R)