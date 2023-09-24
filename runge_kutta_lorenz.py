import numpy as np
import matplotlib.pyplot as plt
import math

def dp_singlestep(func, tk, xk, dt, order=4):
    f1 = func(tk, xk)
    f2 = func(tk+dt/5, xk+dt/5*f1)
    f3 = func(tk+dt*3/10, xk+dt*(3/40*f1+9/40*f2))
    f4 = func(tk+dt*4/5, xk+dt*(44/45*f1-56/15*f2+32/9*f3))
    f5 = func(tk+dt*8/9, xk+dt*(19372/6561*f1-25360/2187*f2+64448/6561*f3-212/729*f4))
    f6 = func(tk+dt, xk+dt*(9017/3168*f1-355/33*f2+46732/5247*f3+49/176*f4-5103/18656*f5))
    f7 = func(tk+dt, xk+dt*(35/384*f1-0*f2+500/1113*f3+125/192*f4-2187/6784*f5+11/84*f6))

    if order == 5:
        xout = 35/384*f1-0*f2+500/1113*f3+125/192*f4-2187/6784*f5+11/84*f6
    else:
        xout = 5179/57600*f1-0*f2+7571/16695*f3+393/640*f4-92097/339200*f5+187/2100*f6+1/40*f7
    return xout

def lorenz(t,x):
    sigma = 10
    beta = 8/3
    rho = 28
    dx = [sigma*(x[1]-x[0]),
                     x[0]*(rho-x[2])-x[1],
                     x[0]*x[1]-beta*x[2]]
    return np.asarray(dx)

def versuch(t, x):
    return -200*t*x*x

x_exact = []
t_exact = []
i = -3
while i <= 0:
    x_exact.append(1/(1+100*i*i))
    t_exact.append(i)
    i += 0.01

plt.show()




