
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib
from matplotlib import cm
import yaml
from matplotlib.animation import FuncAnimation
import os
import imageio

def fun(x):
    a = 2**2
    b = 1**2
    val = 1/(1+np.exp(-((x[0]**2)/a+(x[1]**2)/b-1)))
    #val2 = x[0]**2+x[1]**2
    return val

def fun2(x):
    a = 2**2
    b = 1**2
    val = (x[0]**2)/a+(x[1]**2)/b
    return val

nx = 50
ny = 50
x = np.linspace(-4, 4, nx)
y = np.linspace(-4, 4, ny)
X, Y = np.meshgrid(x, y)
xre = X.reshape(-1,)
yre = Y.reshape(-1,)

Zvals = np.zeros((nx*ny,))
Zvals2 = np.zeros((nx*ny,))
for idx in range(nx*ny):
    Zvals[idx]= fun(np.array([xre[idx], yre[idx]]))
    Zvals2[idx]= fun2(np.array([xre[idx], yre[idx]]))

Z = Zvals.reshape(nx,ny)
Z2 = Zvals2.reshape(nx,ny)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')
ax.scatter3D(0,0,s = 80, c = 'r')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

fig2 = plt.figure()
ax2 = plt.axes(projection='3d')
ax2.plot_surface(X, Y, Z2, cmap='viridis')
ax2.scatter3D(0,0,s = 80, c = 'r')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_zlabel('z')

plt.show()
