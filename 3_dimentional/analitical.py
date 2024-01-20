import math 
import numpy as np 
from matplotlib import cm
import matplotlib.pyplot as plt

def solution(x, y, t) :
    return t * x * y + np.cos(np.pi * math.sqrt(2) * t) * np.sin(np.pi * x / 2) * np.sin(np.pi * y / 2)



def get_analitical()  :
    M = 30
    N = 60  
    vec_x = np.linspace(0.0, 2.0, M + 1)
    vec_t = np.linspace(0.0, 1.0, N + 1)
    X, T = np.meshgrid(vec_x, vec_t)
    Z = solution(X, T)
    return X, T, Z

M = 30
N = 60
t = 1
X, Y = np.meshgrid(np.linspace(0.0, 2.0, M + 1), np.linspace(0.0, 1.0, N + 1))
Z = solution(X, Y, t)
print(X.shape)
print(Y.shape)
print(Z.shape)

plt.figure("t = " + str(t))
ax1 = plt.axes(projection='3d')
ax1.plot_surface(X, Y, Z, rstride=1, cstride=1,
                cmap=cm.coolwarm, edgecolor='none')


ax1.set_xlim(0.0, 2.0)
ax1.set_ylim(0.0, 1.0)
ax1.set_xlabel('X-axis')
ax1.set_ylabel('Y-axis')
ax1.set_zlabel('U-axis')


plt.show()

