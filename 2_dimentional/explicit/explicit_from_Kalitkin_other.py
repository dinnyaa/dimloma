import math 
import numpy as np 
from matplotlib import cm
import matplotlib.pyplot as plt

def solution(x, t) :
    return 2.0 * t + (np.power(x, 2) - 2) * np.sin(t) + np.cos(math.pi * t) * np.cos(math.pi * x)

def f(x, t) :
    return -np.power(x, 2) * np.sin(t)



def change_edges(l, x_arr, h) :
    x_arr = np.insert(x_arr, 0, -h / 2)
    x_arr = np.insert(x_arr, 2, h / 2)
    x_arr = np.insert(x_arr, -1, l - h / 2)
    x_arr = np.append(x_arr, l + h / 2)
    return x_arr
    
def scheme_solution(T, l, M, N) :

    if M <= 0 or N <= 0 :
        raise ValueError("Invalid step arguments")
    h_x = l / M 
    h_t = T / N
    if h_t >= h_x :
        raise ValueError("Stability condition not satisfied")
    
    r = float(np.power(float(h_t / h_x), 2))
    x_arr = np.linspace(0.0, l, M + 1)
    t_arr = np.linspace(0.0, T, N + 1)

    x_arr = change_edges(l, x_arr, h_x)
    M = M + 4
    v = np.zeros((N + 1, M + 1), dtype=float)
    v[0] = np.cos(x_arr * math.pi)
    v[1] = v[0] + h_t * (np.power(x_arr, 2) + (h_t/2)  * (-np.pi * np.sin(np.pi * x_arr) + f(x_arr, 0)))

    for n in range(1, N) :  
        for m in range(1, M) :
            if m == 2 or m == M - 3 :
                rr = float(np.power(float(h_t / (1.5 * h_x)), 2))
                v[n+1][m] = rr * v[n][m-1] + 2.0 * (1.0 - rr) * v[n][m] + rr * v[n][m+1] - v[n-1][m] + np.power(h_t, 2) * f(x_arr[m], t_arr[n])
            else :
                v[n+1][m] = r * v[n][m-1] + 2.0 * (1.0 - r) * v[n][m] + r * v[n][m+1] - v[n-1][m] + np.power(h_t, 2) * f(x_arr[m], t_arr[n])
        v[n+1][0] = v[n+1][1]
        v[n+1][M] = v[n+1][M-1] + 2 * h_x * np.sin(t_arr[n+1])

    return v

def getter(M,N) :
    v = scheme_solution(1.0, 1.0, M, N)
    vec_x = np.linspace(0.0, 1.0, M + 1)
    vec_t = np.linspace(0.0, 1.0, N + 1)
    X, T = np.meshgrid(vec_x, vec_t)
    return X, T, v


# T = 1.0
# M = 30
# N = 60
# l = 1.0
# v = scheme_solution(T, l, M, N)
# X, Y = np.meshgrid(np.linspace(0.0, 1.0, M + 1), np.linspace(0.0, 1.0, N + 1))
# plt.figure("explicit")
# ax1 = plt.axes(projection='3d')
# ax1.plot_surface(X, Y, v, rstride=1, cstride=1,
#                 cmap=cm.coolwarm, edgecolor='none')

# ax1.set_xlim(0.0, 1.0)
# ax1.set_ylim(0.0, 1.0)
# ax1.set_xlabel('X-axis')
# ax1.set_ylabel('T-axis')
# ax1.set_zlabel('V-axis')


# plt.show()
