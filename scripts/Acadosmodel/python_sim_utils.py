import numpy as np
import matplotlib.pyplot as plt
import yaml

class plotter():

    def __init__(self, table, smax):
        self.fig, self.ax = plt.subplots(1, figsize=(20,20))
        self.r = 0.3
        self.smax = smax
        maxidx = np.floor(smax * 100).astype(np.int32)
        #downsample
        self.downsampling = 8
        self.coords_full = table[:, 2:4]
        self.coords = table[:maxidx:self.downsampling, 2:4]
        self.phis = table[:maxidx:self.downsampling, 4]
        self.tvals = table[:maxidx:self.downsampling, 1]
        self.svals = table[:maxidx:self.downsampling, 0]
        self.cos_phi = table[:maxidx:self.downsampling, 5]
        self.sin_phi = table[:maxidx:self.downsampling, 6]
        self.gvals = table[:maxidx:self.downsampling, 7]
        self.normaldir = np.vstack((-self.sin_phi, self.cos_phi)).T

    def plot_track(self):
        print
        self.ax.plot(self.coords[:,0],self.coords[:,1], color = 'k')
        self.ax.scatter(self.coords[::4,0],self.coords[::4,1], color = 'r')
        self.ax.plot(self.coords[:,0] + self.r * self.normaldir[:,0] ,\
        self.coords[:,1]+ self.r * self.normaldir[:,1], linestyle = '--', color = 'k')
        self.ax.plot(self.coords[:,0] - self.r * self.normaldir[:,0] ,\
        self.coords[:,1]- self.r * self.normaldir[:,1], linestyle = '--', color = 'k')
        plt.show(block=False)
        #plt.pause(0.001) # Pause for interval seconds.
        #input("hit [enter] to continue.")


    def plot_traj(self, waypoints):
        self.ax.scatter(waypoints[:,0], waypoints[:,1], s = 10, color = 'r')
        self.fig.canvas.draw()
        #plt.show(block=False)
        #plt.pause(0.001) # Pause for interval seconds.
        #input("hit [enter] to continue.")
        #plt.close('all') # all open plots are correctly closed after each run

    def plot_horizon(self, thetavals, xval):
        idxp = 100 * thetavals
        idxp= idxp.astype(np.int32)
        x_theta_vals = self.coords_full[idxp,0]
        y_theta_vals = self.coords_full[idxp,1]
        self.theta_horz = self.ax.scatter(x_theta_vals, y_theta_vals, marker = 'x', color = 'g')
        self.pos_horz = self.ax.scatter(xval[:,0], xval[:,1], marker = 'D', color = 'b')
        self.fig.canvas.draw()

    def clear_horizion(self):
        self.theta_horz.remove()
        self.pos_horz.remove()
        self.fig.canvas.draw()

def plot_pajecka(modelparams):
    #loadparameters
    with open(modelparams) as file:
        params = yaml.load(file, Loader= yaml.FullLoader)

    m = params['m'] #[kg]
    lf = params['lf'] #[m]
    lr = params['lr'] #[m]
    Iz = params['Iz'] #[kg*m^3]

    #pajecka and motor coefficients
    Bf = params['Bf']
    Br = params['Br']
    Cf = params['Cf']
    Cr = params['Cr']
    Cm1 = params['Cm1']
    Cm2 = params['Cm2']
    Cr = params['Cr']
    Cd = params['Cd']
    Df = params['Df']
    Dr = params['Dr']

    vy_vals = np.linspace(-4, 4, 100)
    vx = 2
    omega = 0
    delta = 0

    ffy_vals = Df*np.sin(Cf*np.arctan(-Bf*np.arctan(vy_vals/vx)))
    fry_vals = Dr*np.sin(Cr*np.arctan(Br*np.arctan(-vy_vals/vx)))

    fig, (ax1, ax2) = plt.subplots(nrows = 2, ncols = 1)
    ax1.plot(vy_vals,ffy_vals)
    ax1.set_xlabel("lateral velocity vy [m/s]")
    ax1.set_ylabel("lateral tire force Ffy [N]")
    ax2.plot(vy_vals, fry_vals)
    ax2.set_xlabel("lateral velocity vy [m/s]")
    ax2.set_ylabel("lateral tire force Fry [N]")
