'''
// MIT License

// Copyright (c) 2020 Peter Werner

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib
from matplotlib import cm
import yaml

class plotter():

    def __init__(self, table, smax, r, lencar):
        #trackplot
        self.fig, self.ax = plt.subplots(1, figsize=(20,12))
        #input and state plot
        self.fig2, (self.ax1, self.ax2) = plt.subplots(nrows = 2, ncols = 1, figsize=(10,10))

        self.r = r
        self.lencar = lencar
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
        self.trajplots = []

    def plot_track(self):
        self.ax.plot(self.coords[:,0],self.coords[:,1], color = 'k')
        self.ax.scatter(self.coords[::4,0],self.coords[::4,1], color = 'r')
        self.ax.plot(self.coords[:,0] + self.r * self.normaldir[:,0] ,\
        self.coords[:,1]+ self.r * self.normaldir[:,1], linestyle = '--', color = 'k')
        self.ax.plot(self.coords[:,0] - self.r * self.normaldir[:,0] ,\
        self.coords[:,1]- self.r * self.normaldir[:,1], linestyle = '--', color = 'k')
        self.ax.set_xlim([-1,5])
        self.ax.set_ylim([-0.5,3.1])
        #plt.show(block=False)
        #plt.pause(2) # Pause for interval seconds.
        #input("hit [enter] to continue.")


    def plot_traj(self, xvals):

        heatmap = self.ax.scatter(xvals[:,0], xvals[:,1], s = 40, c=xvals[:,3], cmap=cm.jet, edgecolor='none', marker='o')
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
        width = self.lencar
        height = width/2


        phi0 = xval[0,2]
        #define transform
        tr = matplotlib.transforms.Affine2D().rotate_deg_around(xval[0,0], xval[0,1], phi0*180/3.14159)
        ts = self.ax.transData
        t = tr + ts

        carrect = patches.Rectangle((xval[0,0]-width/2, xval[0,1]-height/2),width, height,\
        linewidth=1,edgecolor='k',facecolor='none', transform= t)
        self.car = self.ax.add_patch(carrect)
        center = np.array([xval[0,0], xval[0,1]])
        front = center + width/2*np.array([np.cos(phi0), np.sin(phi0)])
        self.cardir = self.ax.plot([center[0], front[0]], [center[1], front[1]], linewidth = 1, color = 'r')

        self.theta_horz = self.ax.scatter(x_theta_vals, y_theta_vals, marker = 'x', color = 'g')
        self.pos_horz = self.ax.scatter(xval[:,0], xval[:,1], marker = 'D', color = 'b')
        for idx in range(len(x_theta_vals)):
            connector = self.ax.plot([x_theta_vals[idx], xval[idx,0]],\
                                                [ y_theta_vals[idx],  xval[idx,1]],\
                                                 linestyle = '--', color = 'gray')
            self.connect_horz.append(connector)
        self.fig.canvas.draw()
        plt.show(block = False)
        plt.pause(0.01)


    def plot_input_state_traj(self, zval, varnames):

        uval = zval[:,:3]
        #state trajectories
        N = len(zval)
        time = np.arange(N)
        for idx in range(zval.shape[1]-3):
            temp = self.ax1.step(time, zval[:,idx + 3], label = varnames[idx + 3], where='post')
            self.trajplots.append(temp)
        self.ax1.legend()
        self.ax1.set_title("State Trajectories")
        self.ax1.set_xlabel("time [t/Ts]")
        max = np.max(zval[:,3:])
        min = np.min(zval[:,3:])
        self.ax1.set_ylim([min-0.1,max+0.1])

        for idx in range(3):
            temp = self.ax2.step(time, uval[:,idx], label = varnames[idx], where='post')
            self.trajplots.append(temp)
        self.ax2.legend()
        self.ax2.set_title("Input Trajectories")
        self.ax2.set_xlabel("time [t/Ts]")
        max = np.max(uval)
        min = np.min(uval)
        self.ax2.set_ylim([min-0.1,max+0.1])
        self.fig2.canvas.draw()

    def plot_input_state_traj(self, xval, uval, xvarnames, uvarnames):
        #state trajectories
        N = len(xval)
        time = np.arange(N)
        for idx in range(xval.shape[1]):
            temp = self.ax1.step(time, xval[:,idx], label = xvarnames[idx], where='post')
            self.trajplots.append(temp)
        self.ax1.legend()
        self.ax1.set_title("State Trajectories")
        self.ax1.set_xlabel("time [t/Ts]")
        max = np.max(xval[:,:])
        min = np.min(xval[:,:])
        self.ax1.set_ylim([min-0.1,max+0.1])

        for idx in range(3):
            temp = self.ax2.step(time, uval[:,idx], label = uvarnames[idx], where='post')
            self.trajplots.append(temp)
        self.ax2.legend()
        self.ax2.set_title("Input Trajectories")
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
        #self.car.remove()
        self.fig.canvas.draw()

    def clear_input_state_traj(self):
        nrtraj = len(self.trajplots)
        for idxtraj in range(nrtraj):
            self.trajplots[idxtraj][0].remove()
        self.trajplots = []
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
