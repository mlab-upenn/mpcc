import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from python_sim_utils import replay_plotter
import pickle
import sys


def main_1ag(filename):
    print("[INFO] Move slider at the bottom of the window to advance time")
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
    trk_plt = replay_plotter(track_lu_table, smax, r, lencar)
    trk_plt.plot_track()
    trk_plt.plot_horizon(zdata[0,:,zvars.index('theta')], zdata[0,:, 3:6])
    #trk_plt.plot_input_state_traj(zdata[0,:,:], zvars)
    axcolor = 'lightgoldenrodyellow'
    axtime = plt.axes([0.1, 0.05, 0.75, 0.03], facecolor=axcolor)
    stime = Slider(axtime, 'Time', 0, zdata.shape[0]-1, valinit=0, valstep = 1)

    def update(val):
        trk_plt.clear_horizion()
        #trk_plt.clear_input_state_traj()
        time = stime.val
        trk_plt.plot_horizon(zdata[time,:,zvars.index('theta')], zdata[time,:, 3:6])
        #print(zdata[time,0,:])
        #trk_plt.plot_input_state_traj(zdata[time,:,:], zvars)
        #print(time)

    stime.on_changed(update)
    plt.show()


def main_2ag(filename):
    print("[INFO] Move slider at the bottom of the window to advance time")
    with open(filename, 'rb') as handle:
        data = pickle.load(handle)

    smax = data['smax']
    lencar = data['lencar']
    r = data['r']
    agent1 = data['agent1']
    agent2 = data['agent2']
    zdata_1 = agent1['zdata']
    zdata_2 = agent2['zdata']
    track_name = "tracks/"+data['track']+"_lutab.csv"
    zvars = ['ddot', 'deltadot', 'thetadot', 'posx', 'posy', 'phi', 'vx', 'vy', 'omega', 'd', 'delta', 'theta']

    plt.plot(np.arange(len(zdata_1[:,0,0])),zdata_1[:,0,zvars.index('theta')])
    # load array
    track_lu_table = np.loadtxt(track_name, delimiter=',')
    print("Experiment with N = ", zdata_1.shape[0]," samples loaded")
    trk_plt = replay_plotter(track_lu_table, smax, r, lencar)
    trk_plt.plot_track()
    trk_plt.plot_agents(zdata_1[0,:,:], zdata_2[0,:,:])
    #trk_plt.plot_input_state_traj(zdata[0,:,:], zvars)
    axcolor = 'lightgoldenrodyellow'
    axtime = plt.axes([0.1, 0.05, 0.75, 0.03], facecolor=axcolor)
    stime = Slider(axtime, 'Time', 0, zdata_1.shape[0]-1, valinit=0, valstep = 1)

    def update(val):
        trk_plt.clear_agents()
        #trk_plt.clear_input_state_traj()
        time = stime.val
        trk_plt.plot_agents(zdata_1[time,:,:], zdata_2[time,:,:])
        #print(zdata[time,0,:])
        #trk_plt.plot_input_state_traj(zdata[time,:,:], zvars)
        #print(time)

    stime.on_changed(update)
    plt.show()


if __name__ == "__main__":
    replayfilename = sys.argv[1]
    filestart = replayfilename[8:]
    nragents = int(filestart[0])
    if nragents == 1:
        main_1ag(replayfilename)
    elif nragents == 2:
        main_2ag(replayfilename)
    else :
        inp = input("Press 1 for 1 agent or 2 for 2 agents")
