import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import yaml

class plotter():

    def __init__(self, table, smax):
        #trackplot
        self.fig, self.ax = plt.subplots(1, figsize=(20,20))
        #input and state plot
        self.fig2, (self.ax1, self.ax2) = plt.subplots(nrows = 2, ncols = 1, figsize=(10,10))

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
        self.ax.set_xlim([-1,5])
        self.ax.set_ylim([-0.5,3])
        #plt.show(block=False)
        #plt.pause(0.001) # Pause for interval seconds.
        #input("hit [enter] to continue.")


    def plot_traj(self, xvals):

        heatmap = self.ax.scatter(xvals[10:,0], xvals[10:,1], s = 10, c=xvals[10:,3], cmap=cm.rainbow, edgecolor='none', marker='o')
        cbar = self.fig.colorbar(heatmap, fraction=0.035)
        cbar.set_label("velocity in [m/s]")
        self.fig.canvas.draw()
        #plt.show(block=False)
        #plt.pause(0.001) # Pause for interval seconds.
        #input("hit [enter] to continue.")
        #plt.close('all') # all open plots are correctly closed after each run

    def plot_horizon(self, thetavals, xval):
        self.connect_horz = []
        idxp = 100 * thetavals
        idxp= idxp.astype(np.int32)
        x_theta_vals = self.coords_full[idxp,0]
        y_theta_vals = self.coords_full[idxp,1]
        self.theta_horz = self.ax.scatter(x_theta_vals, y_theta_vals, marker = 'x', color = 'g')
        self.pos_horz = self.ax.scatter(xval[:,0], xval[:,1], marker = 'D', color = 'b')
        for idx in range(len(x_theta_vals)):
            connector = self.ax.plot([x_theta_vals[idx], xval[idx,0]],\
                                                [ y_theta_vals[idx],  xval[idx,1]],\
                                                 linestyle = '--', color = 'gray')
            self.connect_horz.append(connector)
        self.fig.canvas.draw()

    def plot_input_state_traj(self, xval, uval):
        #state trajectories
        N = len(xval)
        time = np.arange(N)
        self.xplot = self.ax1.step(time, xval[:,0], where='post')
        self.yplot = self.ax1.step(time, xval[:,1], where='post')
        self.phiplot = self.ax1.step(time, xval[:,2], where='post')
        self.vxplot = self.ax1.step(time, xval[:,3], where='post')
        self.vyplot = self.ax1.step(time, xval[:,4], where='post')
        self.omegaplot = self.ax1.step(time, xval[:,5], where='post')
        self.thetaplot = self.ax1.step(time, xval[:,6], where='post')
        self.ax1.legend(['x', 'y', 'phi','vx','theta','d','delta'])
        self.ax1.set_xlabel("time [t/Ts]")
        max = np.max(xval[:,:7])
        min = np.min(xval[:,:7])
        self.ax1.set_ylim([min-0.1,max+0.1])

        self.ddotplot = self.ax2.step(time, uval[:,0], where='post')
        self.deltadotplot = self.ax2.step(time, uval[:,1], where='post')
        self.thetadotplot = self.ax2.step(time, uval[:,2], where='post')
        self.ax2.legend(['ddot', 'deltadot', 'thetadot'])
        self.ax2.set_xlabel("time [t/Ts]")
        max = np.max(uval)
        min = np.min(uval)
        self.ax2.set_ylim([min-0.1,max+0.1])
        self.fig2.canvas.draw()

    def clear_horizion(self):
        self.theta_horz.remove()
        self.pos_horz.remove()
        for idx in range(len(self.connect_horz)):
            connector = self.connect_horz[idx]
            connector[0].remove()
        self.fig.canvas.draw()

    def clear_input_state_traj(self):
        self.xplot[0].remove()
        self.yplot[0].remove()
        self.vxplot[0].remove()
        self.vyplot[0].remove()
        self.phiplot[0].remove()
        self.omegaplot[0].remove()
        self.thetaplot[0].remove()
        self.ddotplot[0].remove()
        self.deltadotplot[0].remove()
        self.thetadotplot[0].remove()
        self.fig2.canvas.draw()

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

def compute_objective(Ts, Qc, Ql, Q_theta, R_d, R_delta, theta, theta_hat, x, u, xt, phit):
    xt_hat = xt[:,0] + np.multiply(np.cos(phit), (theta-theta_hat))
    yt_hat = xt[:,1] + np.multiply(np.sin(phit), (theta-theta_hat))
    e_cont = np.multiply(np.sin(phit), (xt_hat- x[:,0])) - np.multiply(np.cos(phit), (yt_hat- x[:,1]))
    e_lag = np.multiply(np.cos(phit), (xt_hat- x[:,0])) + np.multiply(np.sin(phit), (yt_hat- x[:,1]))
    thetadot = u[:,2]
    deltadot = u[:,1]
    ddot = u[:,0]

    stages = len(x)
    objective = 0
    for stageidx in range(stages):
        objective += (e_cont[stageidx]**2)*Qc + (e_lag[stageidx]**2)*Ql - Q_theta*thetadot[stageidx] +\
                    (ddot[stageidx]**2)*R_d + (deltadot[stageidx]**2)*R_delta

    objective = Ts * objective
    return objective
