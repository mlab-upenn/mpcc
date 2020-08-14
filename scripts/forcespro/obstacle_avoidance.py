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
#from generate_solver_dynamic import get_forces_solver_dynamic
#from generate_col_avoid_solver import get_col_avoid_solver
from racing_agent import racer
import forcespro.nlp
from python_sim_utils import   plotter, plot_pajecka, compute_objective
import matplotlib.pyplot as plt
import InterpolateTrack
import yaml
import sys


def main():
    np.set_printoptions(precision=4)
    np.set_printoptions(threshold=sys.maxsize)
    # model parameters
    modelparams = "parameters/modelparams.yaml"
    solverparams = "parameters/solverparams.yaml"

    #load global constant model parameters
    with open(modelparams) as file:
        params = yaml.load(file, Loader= yaml.FullLoader)
    lf = params['lf'] #[m]
    lr = params['lr'] #[m]
    lencar = (lf+lr)
    widthcar = lencar/2

    #sim parameters
    with open(solverparams) as file:
        params = yaml.load(file, Loader= yaml.FullLoader)
    Tsim = 15
    Tf = params['Tf']
    N = params['N']
    Nsim = np.int(np.floor(N/Tf*Tsim))



    track_lu_table, smax = InterpolateTrack.generatelookuptable("tracks/slider")
    r = 0.2 #trackwidth
    track = {"track_lu_table": track_lu_table,
             "smax": smax,
             "r": r}
    trk_plt = plotter(track_lu_table, smax, r, lencar)

    agent = racer(track, Tsim, "agent_1", modelparams, solverparams)

    #plot_pajecka(paramfile)
    trk_plt.plot_track()

    #starting position in track startidx = theta0[m] * 100 [pts/m]
    trackvars = ['sval', 'tval', 'xtrack', 'ytrack', 'phitrack', 'cos(phi)', 'sin(phi)', 'g_upper', 'g_lower']
    startidx = 1300
    xt0 = track_lu_table[startidx,trackvars.index('xtrack')]
    yt0 = track_lu_table[startidx,trackvars.index('ytrack')]
    phit0 = track_lu_table[startidx,trackvars.index('phitrack')]
    theta_hat0 = track_lu_table[startidx,trackvars.index('sval')]
    #initial condition
    zvars = ['ddot', 'deltadot', 'thetadot', 'posx', 'posy', 'phi', 'vx', 'vy', 'omega', 'd', 'delta', 'theta']
    xvars = ['posx', 'posy', 'phi', 'vx', 'vy', 'omega', 'd', 'delta', 'theta']
    xinit = np.array([xt0, yt0, phit0, 0.2, 0.0, 0, 0, 0, theta_hat0])

    #static obstacle
    ob_idx = 400
    phi_ob = track_lu_table[ob_idx, trackvars.index('phitrack')]
    x_ob = track_lu_table[ob_idx,trackvars.index('xtrack')] - 0.5 * r * np.sin(phi_ob)
    y_ob = track_lu_table[ob_idx,trackvars.index('ytrack')] + 0.5 * r * np.cos(phi_ob)
    l_ob = lencar*2
    w_ob = l_ob*2
    obstacleinfo = {"phi_ob": phi_ob,
                    "x_ob": x_ob,
                    "y_ob": y_ob,
                    "l_ob": l_ob,
                    "w_ob": w_ob}

    trk_plt.plot_static_obstacle(x_ob, y_ob, phi_ob, l_ob, w_ob)
    z_current = agent.initialize_trajectory(xinit, obstacleinfo, startidx)
    trk_plt.plot_horizon(z_current[:,zvars.index('theta')], z_current[:, 3:6])
    plt.pause(0.1)
    trk_plt.clear_horizion()
    trk_plt.clear_input_state_traj()
    input("start")
    ##########################SIMULATION#######################################
    for simidx in range(Nsim):
        z_current = agent.update(obstacleinfo)
        #plotting result
        trk_plt.plot_horizon(z_current[:,zvars.index('theta')], z_current[:, 3:6])
        trk_plt.plot_input_state_traj(z_current, zvars)

        #plt.pause(0.01)
        #input("hit [enter] to continue.")
        #plt.pause(0.1)
        trk_plt.clear_horizion()
        trk_plt.clear_input_state_traj()

    ###############################/SIMULATION##################################

    zinit_vals, z_data = agent.return_sim_data()
    #trk_plt.animate_result(z_data, N,Tf, 'test.gif')
    trk_plt.plot_traj(zinit_vals[:,3:])
    plt.show()
    #np.savetxt("full_sol_x_log.csv", )
        #print(simidx)
    return 0

if __name__ == "__main__":

    main()
