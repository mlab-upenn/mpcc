import numpy as np
import matplotlib.pyplot as plt

class plotter():

    def __init__(self, table):
        self.fig, self.ax = plt.subplots(1)
        self.r = 0.3
        #downsample
        self.downsampling = 8
        self.coords = table[::self.downsampling, 2:4]
        self.phis = table[::self.downsampling, 4]
        self.tvals = table[::self.downsampling, 1]
        self.svals = table[::self.downsampling, 0]
        self.cos_phi = table[::self.downsampling, 5]
        self.sin_phi = table[::self.downsampling, 6]
        self.gvals = table[::self.downsampling, 7]
        self.normaldir = np.vstack((-self.sin_phi, self.cos_phi)).T

    def plot_track(self):
        self.ax.plot(self.coords[:,0],self.coords[:,1], color = 'k')
        self.ax.plot(self.coords[:,0] + self.r * self.normaldir[:,0] ,\
        self.coords[:,1]+ self.r * self.normaldir[:,1], linestyle = '--', color = 'k')
        self.ax.plot(self.coords[:,0] - self.r * self.normaldir[:,0] ,\
        self.coords[:,1]- self.r * self.normaldir[:,1], linestyle = '--', color = 'k')
        plt.show(block=False)
        plt.pause(0.001) # Pause for interval seconds.
        input("hit [enter] to continue.")


    def plot_traj(self, waypoints):
        self.ax.scatter(waypoints[:,0], waypoints[:,1], s = 10, color = 'r')
        self.fig.canvas.draw()
        #plt.show(block=False)
        #plt.pause(0.001) # Pause for interval seconds.
        #input("hit [enter] to continue.")
        #plt.close('all') # all open plots are correctly closed after each run

def plot_pajecka(modelparams)
