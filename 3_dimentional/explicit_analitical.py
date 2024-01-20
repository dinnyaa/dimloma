import math 
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm


def analitical_solution(x, y, t) :
    return t * x * y + np.cos(np.pi * math.sqrt(2) * t) * np.sin((np.pi * x) / 2) * np.sin((np.pi * y) / 2)

def initial1(x,y) :
    return np.sin((np.pi * x) / 2) * np.sin((np.pi * y) / 2)

def initial2(x,y) :
    return x * y

def boundry1(y,t) :
    return 0

def boundry2(y,t) :
    return 2 * y * t

def boundry3(x,t) :
    return 0

def boundry4(x,t) :
    return t * x

def scheme1_solution(T, X, Y, K, M, N) :

    a = 2.0
    if M <= 0 or N <= 0 or K <= 0:
        raise ValueError("Invalid step arguments")
    h_x = X / M 
    h_y = Y / N
    h_t = T / K
    if h_t >= 1.0 / math.sqrt((a/h_x) ** 2 + (a/h_y) ** 2) :
       raise ValueError("Stability condition not satisfied")
    x_arr = np.linspace(0.0, X, M + 1)
    y_arr = np.linspace(0.0, Y, N + 1)
    t_arr = np.linspace(0.0, T, K + 1)
    
    v = np.zeros((K + 1, N + 1, M + 1), dtype=float)
    mat_x, mat_y = np.meshgrid(x_arr, y_arr)
    v[0] = initial1(mat_x, mat_y)
    v[1] = v[0] + h_t * initial2(mat_x, mat_y)


    mat_y, mat_t = np.meshgrid(y_arr, t_arr[2:])
    v[2:, :, 0] = boundry1(mat_y, mat_t)
    v[2:, :, M] = boundry2(mat_y, mat_t)

    for k in range(1, K) :
        for n in range (1, N) :
            for m in range (1, M) :
                v[k+1][n][m] = 2 * v[k][n][m] - v[k-1][n][m] \
                    + ((a * h_t / h_x) ** 2) * (v[k][n][m+1] - 2 * v[k][n][m] + v[k][n][m-1]) \
                        + ((a * h_t / h_y) ** 2) * (v[k][n+1][m] - 2 * v[k][n][m] + v[k][n-1][m])
        v[k+1, 0, :] = boundry3(x_arr, t_arr[k+1])
        v[k+1, N, :] = v[k+1, N-1, :] + h_y * boundry4(x_arr, t_arr[k+1])
       # v[k+1, N, :] =  (4*v[k+1, N-1, :] -  v[k+1, N-2, :] + 2 *  h_y * boundry4(x_arr, t_arr[k+1])) / 3.0
    return v

def scheme2_solution(T, X, Y, K, M, N) :

    a = 2.0
    if M <= 0 or N <= 0 or K <= 0:
        raise ValueError("Invalid step arguments")
    h_x = X / M 
    h_y = Y / N
    h_t = T / K
    if h_t >= 1.0 / (math.sqrt((a/h_x) ** 2 + (a/h_y) ** 2)):
        raise ValueError("Stability condition not satisfied")
       print("nonono")
    x_arr = np.linspace(0.0, X, M + 1)
    y_arr = np.linspace(0.0, Y, N + 1)
    t_arr = np.linspace(0.0, T, K + 1)
    
    v = np.zeros((K + 1, N + 1, M + 1), dtype=float)
    mat_x, mat_y = np.meshgrid(x_arr, y_arr)
    v[0] = initial1(mat_x, mat_y)
    tmp = -(a ** 2) * ((np.power(np.pi, 2)/2.0) * initial1(mat_x, mat_y))
    v[1] = v[0] + h_t * (initial2(mat_x, mat_y) + (h_t / 2.0) * tmp)
    #v[1] = v[0] + h_t * initial2(mat_x, mat_y)

    

    mat_y, mat_t = np.meshgrid(y_arr, t_arr[2:])
    v[2:, :, 0] = boundry1(mat_y, mat_t)
    v[2:, :, M] = boundry2(mat_y, mat_t)

    for k in range(1, K) :
        for n in range (1, N) :
            for m in range (1, M) :
                v[k+1][n][m] = 2 * v[k][n][m] - v[k-1][n][m] \
                    + ((a * h_t / h_x) ** 2) * (v[k][n][m+1] - 2 * v[k][n][m] + v[k][n][m-1]) \
                        + ((a * h_t / h_y) ** 2) * (v[k][n+1][m] - 2 * v[k][n][m] + v[k][n-1][m]) +np.power(h_t, 2) * f(x_arr[m], y_arr[n], t_arr[k])
        v[k+1, 0, :] = boundry3(x_arr, t_arr[k+1])
        v[k+1, N, :] = v[k+1, N-1, :] + h_y * boundry4(x_arr, t_arr[k+1])
       # v[k+1, N, :] =  (4*v[k+1, N-1, :] -  v[k+1, N-2, :] + 2 *  h_y * boundry4(x_arr, t_arr[k+1])) / 3.0
    return v




T = 1.0
M = 20
N = 30
K = 80
point = 40
t = np.linspace(0.0, T, K + 1)[point]
print(t)
X, Y = np.meshgrid(np.linspace(0.0, 2.0, M + 1), np.linspace(0.0, 1.0, N + 1))
V1 = scheme1_solution(1.0, 2.0, 1.0, K,  M, N)
V2 = scheme2_solution(1.0, 2.0, 1.0, K,  M, N)

U = analitical_solution(X, Y, t)


plt.figure("anal")
ax1 = plt.axes(projection='3d')
ax1.plot_surface(X, Y, U, rstride=1, cstride=1,
                cmap=cm.coolwarm, edgecolor='none')

ax1.set_xlim(0.0, 2.0)
ax1.set_ylim(0.0, 1.0)
ax1.set_xlabel('X-axis')
ax1.set_ylabel('Y-axis')
ax1.set_zlabel('V-axis')


plt.figure("motarkum")
ax2 = plt.axes(projection='3d')
ax2.plot_surface(X, Y, V1[point], rstride=1, cstride=1,
                cmap=cm.coolwarm, edgecolor='none')

ax2.set_xlim(0.0, 2.0)
ax2.set_ylim(0.0, 1.0)
ax2.set_xlabel('X-axis')
ax2.set_ylabel('Y-axis')
ax2.set_zlabel('U-axis')



plt.figure("sxalanq")
ax3 = plt.axes(projection='3d')
ax3.plot_surface(X, Y, U-V1[point], rstride=1, cstride=1,
                cmap=cm.coolwarm, edgecolor='none')

ax3.set_xlim(0.0, 2.0)
ax3.set_ylim(0.0, 1.0)
ax3.set_xlabel('X-axis')
ax3.set_ylabel('Y-axis')
ax3.set_zlabel('S-axis')
plt.show()


