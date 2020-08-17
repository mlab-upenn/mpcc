import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from python_sim_utils import replay_plotter
import pickle
import sys
'''
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
t = np.arange(0.0, 1.0, 0.001)
a0 = 5
f0 = 3
delta_f = 5.0
s = a0 * np.sin(2 * np.pi * f0 * t)
l, = plt.plot(t, s, lw=2)
ax.margins(x=0)

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
axamp = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

sfreq = Slider(axfreq, 'Freq', 0.1, 30.0, valinit=f0, valstep=delta_f)
samp = Slider(axamp, 'Amp', 0.1, 10.0, valinit=a0)


def update(val):
    amp = samp.val
    freq = sfreq.val
    l.set_ydata(amp*np.sin(2*np.pi*freq*t))
    fig.canvas.draw_idle()


sfreq.on_changed(update)
samp.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    sfreq.reset()
    samp.reset()
button.on_clicked(reset)

rax = plt.axes([0.025, 0.5, 0.15, 0.15], facecolor=axcolor)
radio = RadioButtons(rax, ('red', 'blue', 'green'), active=0)


def colorfunc(label):
    l.set_color(label)
    fig.canvas.draw_idle()
radio.on_clicked(colorfunc)

plt.show()
'''

def main(filename):
    with open(filename, 'rb') as handle:
        data = pickle.load(handle)

    smax = data['smax']
    lencar = data['lencar']
    r = data['r']
    zinit = data['zinit']
    zdata = data['zdata']
    track_name = "tracks/"+data['track']+"_lutab.csv"
    zvars = ['ddot', 'deltadot', 'thetadot', 'posx', 'posy', 'phi', 'vx', 'vy', 'omega', 'd', 'delta', 'theta']

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



if __name__ == "__main__":
    replayfilename = sys.argv[1]
    main(replayfilename)
