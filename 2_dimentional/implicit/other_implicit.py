import math 
import numpy as np 
from matplotlib import cm
import matplotlib.pyplot as plt

def solution(x, t) :
    return (0.5 / np.pi) * np.sin(2 * np.pi * t) *  np.sin(2 * np.pi * x) 
    #return np.sin(x) * (t - np.sin(t))


def matrix(v, M, n, r2, d) : # Ax = b

    A = np.zeros((M-1, M-1))
    np.fill_diagonal(A, -(2.0 * r2 * d + 1))
    np.fill_diagonal(A[1:], r2 * d)
    np.fill_diagonal(A[:,1:], r2 * d)


    b = np.empty(M - 1)
    i = 0
    for m in range (1, M) : 
        b[i] = -r2 * (1 - 2 * d)*v[n][m+1] - (2 - 2 * r2 * (1 - 2 * d)) * v[n][m] - r2 * (1- 2 * d) * v[n][m - 1] \
        - r2 * d * v[n-1][m+1] + 2 * r2 * d * v[n-1][m] - r2 * d * v[n-1][m - 1] + v[n - 1][m] #- (tao ** 2) * f(x_arr[m], t_arr[n])
        i+=1
    sol = np.linalg.solve(A, b)
    print(sol)
    return sol

def explicit_solution(T, l, M, N) :

    if M <= 0 or N <= 0 :
        raise ValueError("Invalid step arguments")
    h_x = l / M 
    h_t = T / N
    x_arr = np.linspace(0.0, l, M + 1)
    t_arr = np.linspace(0.0, T, N + 1)

    r = float(np.power(float(h_t / h_x), 2))
    #if r < 0 :
    if h_t >= h_x :
        #raise ValueError("Stability condition not satisfied")
        print("nonono" + str(r))
        
    v = np.zeros((N + 1, M + 1), dtype=float)
    v[1] =h_t * np.sin( 2 * np.pi * x_arr)


    # 2rd shertic sksac
    for n in range(1, N) :  
        for m in range(1, M) : 
            v[n + 1][m] = r * v[n][m - 1] + 2.0 * (1.0 - r) * v[n][m] + r * v[n][m + 1] - v[n - 1][m]
    return v

def implicit_solution(T, l, M, N) :

    if M <= 0 or N <= 0 :
        print("Negative value")
        exit()
    h_x = l / M 
    h_t = T / N
    x_arr = np.linspace(0.0, l, M + 1)
    t_arr = np.linspace(0.0, T, N + 1)

    if h_t >= h_x :
        print("unstable scheme")
        exit()

    r = float(np.power(float(h_t / h_x), 2))
    v = np.zeros((N + 1, M + 1), dtype=float)
    # v[0] = np.sin(x_arr * 3)
    v[1] =h_t * np.sin( 2 * np.pi * x_arr)

    for n in range(1, N) :
        # v[n + 1][0] = 0  
        # v[n + 1][M] = 0 
        v[n + 1, 1 : M] = matrix(v, M, n, r, 0.3)  
   
    return v


T = 1.0
M = 30
N = 50
l = 1.0
v1 = implicit_solution(T, l, M, N)
v2 = explicit_solution(T, l, M, N)
     
vec_x = np.linspace(0.0, 1.0, M + 1)
vec_t = np.linspace(0.0, 1.0, N + 1)
X, T = np.meshgrid(vec_x, vec_t)

Z = solution(X, T)

plt.figure("explicit")
ax1 = plt.axes(projection='3d')
ax1.plot_surface(X, T, Z, rstride=1, cstride=1,
                cmap=cm.coolwarm, edgecolor='none')

plt.figure("implicit")
ax2 = plt.axes(projection='3d')
ax2.plot_surface(X, T, v1, rstride=1, cstride=1,
                cmap=cm.coolwarm, edgecolor='none')

ax1.set_xlim(0.0, 1.0)
ax1.set_ylim(0.0, 1.0)
ax1.set_xlabel('X-axis')
ax1.set_ylabel('T-axis')
ax1.set_zlabel('U-axis')

ax2.set_xlim(0.0, 1.0)
ax2.set_ylim(0.0, 1.0)
ax2.set_xlabel('X-axis')
ax2.set_ylabel('T-axis')
ax2.set_zlabel('V-axis')
plt.show()


