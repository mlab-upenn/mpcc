import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from python_sim_utils import replay_plotter, plotter
import pickle
import sys


filename = "simdata/1_agent_obs_avoid_08_28_16_00_25.pkl"

with open(filename, 'rb') as handle:
    data = pickle.load(handle)

smax = data['smax']
lencar = data['lencar']
r = data['r']
zinit = data['zinit']
zdata = data['zdata']
laptimes = data['laptimes']
track_name = "tracks/"+data['track']+"_lutab.csv"
zvars = ['ddot', 'deltadot', 'thetadot', 'posx', 'posy', 'phi', 'vx', 'vy', 'omega', 'd', 'delta', 'theta']
theta_vals = zinit[:,zvars.index('theta')]
#indexes = np.where(0<zinit<0.01)
laptimes = np.delete(laptimes, laptimes.argmin())
print("laptimes: ", laptimes)

# load array
track_lu_table = np.loadtxt(track_name, delimiter=',')
print("Experiment with N = ", zdata.shape[0]," samples loaded")
trk_plt = plotter(track_lu_table, smax, r, lencar)
trk_plt.plot_track()

#trk_plt.animate_result(zdata[36:399,:,:], 35, 1.0, 'test.gif')
#trk_plt.animate_result(zdata[36:80,:,:], 35, 1.0, 'test.gif')
trk_plt.plot_traj(zinit[36:399,3:])
plt.show()
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib
from matplotlib import cm
import yaml
from matplotlib.animation import FuncAnimation
import os
import imageio


imagenames = []
Tf = 1.5
N = 40
Ts = Tf/N

for idx in range(355):
    imagenames.append('tempdump/f_'+str(idx)+'.jpg')


with imageio.get_writer("small_final.gif", mode='I', duration = Ts) as writer:
    for filename in imagenames:
        image = imageio.imread(filename)
        writer.append_data(image)

ratio = 7/5
figwidth = 12
#statewidth = 10

fig = plt.figure(figsize=(figwidth*2, ratio*figwidth))
grid = plt.GridSpec(2, 2, hspace=0.2, wspace=0.2)
trackax = fig.add_subplot(grid[:, 0])
stateax = fig.add_subplot(grid[0, 1] )
inputax = fig.add_subplot(grid[1, 1] )

# scatter points on the main axes
main_ax.plot(x, y, 'ok', markersize=3, alpha=0.2)

# histogram on the attached axes
x_hist.hist(x, 40, histtype='stepfilled',
            orientation='vertical', color='gray')
x_hist.invert_yaxis()

y_hist.hist(y, 40, histtype='stepfilled',
            orientation='horizontal', color='gray')
y_hist.invert_xaxis()


plt.show()
'''
