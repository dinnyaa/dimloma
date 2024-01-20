import math 
import numpy as np 
from matplotlib import cm
import matplotlib.pyplot as plt

def solution(x, t) :
    return 2.0 * t + (np.power(x, 2) - 2) * np.sin(t) + np.cos(math.pi * t) * np.cos(math.pi * x)

def f(x, t) :
    return -np.power(x, 2) * np.sin(t)

# սկզբնական պայմանների սահմանում
def initial1(x) :
    return np.cos(x * math.pi)

def initial2(x) :
    return np.power(x, 2)

# եզրային պայմանների սահմանում
def boundry1(t) :
    return 0

def boundry2(t) :
    return 2 * np.sin(t)

def explicit_solution(T, l, M, N) :

    # մուտքային պայմանների ստուգում
    if M <= 0 or N <= 0 :
        raise ValueError("Invalid step arguments")
    
    # քայլերի հաշվարկը
    h_x = l / M 
    h_t = T / N
    r = np.power(h_t / h_x, 2)

    #կայունության պայմանի ստուգում
    if r > 1 :
        raise ValueError("Stability condition not satisfied")
    
    # ցանցի կառուցում
    x_arr = np.linspace(0.0, l, M + 1)
    t_arr = np.linspace(0.0, T, N + 1)

    # լուծումների մատրից
    v = np.zeros((N + 1, M + 1), dtype=float)

    #սկզբանական պայմանների կիրառում 
    v[0] = initial1(x_arr)
    v[1] = v[0] + h_t * initial2(x_arr)

    # բացահայտ սխեմայի կիրառում
    for n in range(1, N) :  
        for m in range(1, M) : 
            v[n+1][m] = r * v[n][m-1] + 2.0 * (1.0 - r) * v[n][m] + r * v[n][m+1] - v[n-1][m] + np.power(h_t, 2) * f(x_arr[m], t_arr[n])
        # եզրային պայմանների կիրառում
        v[n+1][0] = v[n+1][1] + h_x * boundry1(t_arr[n+1])
        v[n+1][M] = v[n+1][M-1] + h_x *  boundry2(t_arr[n+1])
    
    return v


def getter(M, N) :
    X, T = np.meshgrid(np.linspace(0.0, 1.0, M + 1), np.linspace(0.0, 1.0, N+ 1))
    Z = explicit_solution(1.0, 1.0, M, N)
    return X, T, Z


T = 1.0
M = 30
N = 60
l = 1.0
v = explicit_solution(T, l, M, N)

# #լուծման վիզուալիզացիա
X, T = np.meshgrid(np.linspace(0.0, 1.0, M + 1), np.linspace(0.0, 1.0, N + 1))
plt.figure("explicit basic")
ax1 = plt.axes(projection='3d')
ax1.plot_surface(X, T, v, rstride=1, cstride=1,
                cmap=cm.coolwarm, edgecolor='none')

plt.figure("analititicqal")
ax2 = plt.axes(projection='3d')
ax2.plot_surface(X, T, solution(X,T), rstride=1, cstride=1,
                cmap=cm.coolwarm, edgecolor='none')

plt.show()


