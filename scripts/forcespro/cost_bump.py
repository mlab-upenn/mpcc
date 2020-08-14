
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
    offset_rear = -1
    offset_front = 1.5
    mu_rear = np.array([-1,0])
    mu_front = np.array([1,0])
    sigmainv_rear = np.array([[0.5,0],[0,2]])
    sigmainv_front = np.array([[0.3,0],[0,0.3]])

    #front gauss value
    front = np.exp(-(x-mu_rear)@sigmainv_rear@(x-mu_rear))
    back = np.exp(-(x-mu_front)@sigmainv_front@(x-mu_front))
    return front - 0.4*back

nx = 50
ny = 50
x = np.linspace(-4, 4, nx)
y = np.linspace(-4, 4, ny)
X, Y = np.meshgrid(x, y)
xre = X.reshape(-1,)
yre = Y.reshape(-1,)

Zvals = np.zeros((nx*ny,))
for idx in range(nx*ny):
    Zvals[idx]= fun(np.array([xre[idx], yre[idx]]))
Z = Zvals.reshape(nx,ny)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')
ax.scatter3D(0,0,s = 80, c = 'r')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()
