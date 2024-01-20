import math 
import numpy as np 
from matplotlib import cm
import matplotlib.pyplot as plt

def solution(x, t) :
    return 2.0 * t + (np.power(x, 2) - 2) * np.sin(t) + np.cos(math.pi * t) * np.cos(math.pi * x)


def f(x, t) :
    return -np.power(x, 2) * np.sin(t)


# def matrix(v, M, n, r2, d, tao, x_arr, t_arr) : # Ax = b

#     A = np.zeros((M-1, M-1))
#     np.fill_diagonal(A, -(2.0 * r2 * d + 1))
#     np.fill_diagonal(A[1:], r2 * d)
#     np.fill_diagonal(A[:,1:], r2 * d)

#     b = np.empty(M - 1)
#     i = 0
#     for m in range (1, M) : 
#         b[i] = -r2 * (1 - 2 * d)*v[n][m+1] - (2 - 2 * r2 * (1 - 2 * d)) * v[n][m] - r2 * (1- 2 * d) * v[n][m-1] - r2 * d * v[n-1][m+1] + \
#                 2 * r2 * d * v[n-1][m] - r2 * d * v[n-1][m-1] + v[n-1][m] - f(x_arr[m], t_arr[n]) * np.power(tao,2)
#         i+=1
    
#     sol = np.linalg.solve(A, b)
#     # print(sol)
#     return sol

def matrix(v, M, n, r2, d, h_x, h_t, x_arr, t_arr) : # Ax = b
    A = np.zeros((M+1, M+1))
    np.fill_diagonal(A, -(2.0 * r2 * d + 1))
    np.fill_diagonal(A[1:], r2 * d)
    np.fill_diagonal(A[:,1:], r2 * d)
    A[0][0] = 1
    A[0][1] = -1
    A[M][M-1] = 1
    A[M][M] = -1
    b = np.empty(M+1)
    i = 1
    for m in range (1, M) : 
        b[i] = -r2 * (1 - 2 * d)*v[n][m+1] - (2 - 2 * r2 * (1 - 2 * d)) * v[n][m] - r2 * (1- 2 * d) * v[n][m-1] \
                - r2 * d * v[n-1][m+1] + (1 + 2 * r2 * d) * v[n-1][m] - r2 * d * v[n-1][m-1] \
                - f(x_arr[m], t_arr[n]) * np.power(h_t,2)
        i+=1
    b[0] = 0
    b[M] = -2 * h_x * np.sin(t_arr[n+1])  
    
    return np.linalg.solve(A, b)


def scheme_solution(T, l, M, N) :

    if M <= 0 or N <= 0 :
        raise ValueError("Invalid step arguments")
    h_x = l / M
    h_t = T / N
    sigma = 0.5
    r = np.power(h_t / h_x, 2)
    if r * (1-4*sigma) > 1:
        raise ValueError("Stability condition not satisfied")
 

    x_arr = np.linspace(0.0, l, M + 1)
    t_arr = np.linspace(0.0, T, N + 1)
    v = np.zeros((N + 1, M + 1), dtype=float)
 
    v[0] = np.cos(x_arr * math.pi)
    # v[1] = v[0] + h_t * np.power(x_arr, 2)
    v[1] = v[0] + h_t * (np.power(x_arr, 2) + (h_t/2)  * (-np.pi * np.sin(np.pi * x_arr) + f(x_arr, 0)))

    for n in range(1, N) :
        v[n+1] = matrix(v, M, n, r, 0.5, h_x, h_t, x_arr, t_arr) 
    return v


def getter(M, N) :
    X, T = np.meshgrid(np.linspace(0.0, 1.0, M + 1), np.linspace(0.0, 1.0, N+ 1))
    Z = scheme_solution(1.0, 1.0, M, N)
    return X, T, Z

T = 1.0
M = 20
N = 40
l = 1.0
v = scheme_solution(T, l, M, N)
X, T = np.meshgrid(np.linspace(0.0, 1.0, M + 1), np.linspace(0.0, 1.0, N + 1))

Z = solution(X, T)


# plt.figure(1)
# ax1 = plt.axes(projection='3d')
# ax1.plot_surface(X, T, Z, rstride=1, cstride=1,
#                 cmap=cm.coolwarm, edgecolor='none')

# plt.figure(2)
# ax2 = plt.axes(projection='3d')
# ax2.plot_surface(X, T, v, rstride=1, cstride=1,
#                 cmap=cm.coolwarm, edgecolor='none')

# ax1.set_xlim(0.0, 1.0)
# ax1.set_ylim(0.0, 1.0)
# ax1.set_zlim(-1.0, 2.2)

# ax1.set_xlabel('X-axis')
# ax1.set_ylabel('T-axis')
# ax1.set_zlabel('U-axis')

# ax2.set_xlim(0.0, 1.0)
# ax2.set_ylim(0.0, 1.0)
# ax2.set_zlim(-1.0, 2.2)

# ax2.set_xlabel('X-axis')
# ax2.set_ylabel('T-axis')
# ax2.set_zlabel('V-axis')

# plt.show()


