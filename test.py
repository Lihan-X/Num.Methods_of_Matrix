import numpy as np
import math

A = np.asarray([[4,-1,1],[9,-8,9],[11,-11,12]])
B= np.asarray([[0,3,1],[0,4,-2],[2,1,2]])

A = A*np.identity(3)

R = np.array([[-14.77, 13.34, -14.7], [0,2.82,-3.15],[0,0,-0.24]])
Q = np.array([[-0.271, 0.926, 0.264], [-0.61,0.047,-0.791],[-0.745,-0.375,-0.552]])

Q, R=np.linalg.qr(B)
print(Q, R)
