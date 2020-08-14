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
from racing_agent_obs import racer
import forcespro.nlp
from python_sim_utils import   plotter, plot_pajecka, compute_objective
import matplotlib.pyplot as plt
import Bezier
import yaml
import pickle
import sys


def main():
    np.set_printoptions(precision=4)
    np.set_printoptions(threshold=sys.maxsize)

    #state variable ordering
    trackvars = ['sval', 'tval', 'xtrack', 'ytrack', 'phitrack', 'cos(phi)', 'sin(phi)', 'g_upper', 'g_lower']
    zvars = ['ddot', 'deltadot', 'thetadot', 'posx', 'posy', 'phi', 'vx', 'vy', 'omega', 'd', 'delta', 'theta']
    xvars = ['posx', 'posy', 'phi', 'vx', 'vy', 'omega', 'd', 'delta', 'theta']

    # parameters for both agents
    modelparams_1 = "parameters/modelparams_1.yaml"
    solverparams_1 = "parameters/solverparams_1.yaml"
    modelparams_2 = "parameters/modelparams_2.yaml"
    solverparams_2 = "parameters/solverparams_2.yaml"

    #sim parameters
    with open(solverparams_1) as file:
        params = yaml.load(file, Loader= yaml.FullLoader)
    Tsim = 30
    Tf = params['Tf']
    N = params['N']
    Nsim = np.int(np.floor(N/Tf*Tsim))

    track_lu_table, smax = Bezier.generatelookuptable("tracks/sample_track")
    r = 0.2 #trackwidth
    track = {"track_lu_table": track_lu_table,
             "smax": smax,
             "r": r}

    #load car size for plotting
    with open(modelparams_1) as file:
        params = yaml.load(file, Loader= yaml.FullLoader)
    lf = params['lf'] #[m]
    lr = params['lr'] #[m]
    lencar = (lf+lr)
    widthcar = lencar/2
    trk_plt = plotter(track_lu_table, smax, r, lencar)
    trk_plt.plot_track()
    plt.pause(0.1)
    #initialize both agents
    agent_1 = racer(track, Tsim, "agent_1", modelparams_1, solverparams_1)
    agent_1_info = {"phi_ob": 0,
                    "x_ob": 0,
                    "y_ob": 0,
                    "l_ob": lencar,
                    "w_ob": widthcar}
    agent_2 = racer(track, Tsim, "agent_2", modelparams_2, solverparams_2)
    agent_2_info = {"phi_ob": 0,
                    "x_ob": 0,
                    "y_ob": 0,
                    "l_ob": lencar,
                    "w_ob": widthcar}

    #agent 1 trajecotry initialization
    #starting position in track startidx = theta0[m] * 100 [pts/m]
    startidx_1 = 400
    xt0 = track_lu_table[startidx_1,trackvars.index('xtrack')]
    yt0 = track_lu_table[startidx_1,trackvars.index('ytrack')]
    phit0 = track_lu_table[startidx_1,trackvars.index('phitrack')]
    theta_hat0 = track_lu_table[startidx_1,trackvars.index('sval')]
    #initial condition
    xinit_1 = np.array([xt0, yt0, phit0, 0.2, 0.0, 0, 0, 0, theta_hat0])
    agent_1_info["phi_ob"] = xinit_1[xvars.index('phi')]
    agent_1_info["x_ob"] = xinit_1[xvars.index('posx')]
    agent_1_info["y_ob"] = xinit_1[xvars.index('posy')]

    #agent 2  trajectory initialization
    startidx_2 = 500
    xt0 = track_lu_table[startidx_2,trackvars.index('xtrack')]
    yt0 = track_lu_table[startidx_2,trackvars.index('ytrack')]
    phit0 = track_lu_table[startidx_2,trackvars.index('phitrack')]
    theta_hat0 = track_lu_table[startidx_2,trackvars.index('sval')]
    #initial condition
    xinit_2 = np.array([xt0, yt0, phit0, 0.2, 0.0, 0, 0, 0, theta_hat0])
    agent_2_info["phi_ob"] = xinit_2[xvars.index('phi')]
    agent_2_info["x_ob"] = xinit_2[xvars.index('posx')]
    agent_2_info["y_ob"] = xinit_2[xvars.index('posy')]

    #compute trajectory init
    z_current_1 = agent_1.initialize_trajectory(xinit_1, agent_2_info, startidx_1)
    z_current_2 = agent_2.initialize_trajectory(xinit_2, agent_1_info, startidx_2)

    trk_plt.plot_agents(z_current_1, z_current_2)
    trk_plt.plot_static_obstacle(agent_1_info["x_ob"],\
                                 agent_1_info["y_ob"],\
                                 agent_1_info["phi_ob"],\
                                 2*agent_1_info["l_ob"],\
                                 2*agent_1_info["w_ob"]
                                )
    trk_plt.plot_static_obstacle(agent_2_info["x_ob"],\
                                 agent_2_info["y_ob"],\
                                 agent_2_info["phi_ob"],\
                                 2*agent_2_info["l_ob"],\
                                 2*agent_2_info["w_ob"]
                                )
    plt.pause(0.1)
    input("press ENTER")
    trk_plt.clear_agents()
    trk_plt.clear_obstacles()
    ##########################SIMULATION#######################################
    for simidx in range(Nsim):
        agent_1_info["phi_ob"] = z_current_1[1, zvars.index('phi')]
        agent_1_info["x_ob"] = z_current_1[1, zvars.index('posx')]
        agent_1_info["y_ob"] = z_current_1[1, zvars.index('posy')]
        agent_2_info["phi_ob"] = z_current_2[1, zvars.index('phi')]
        agent_2_info["x_ob"] = z_current_2[1, zvars.index('posx')]
        agent_2_info["y_ob"] = z_current_2[1, zvars.index('posy')]


        z_current_1 = agent_1.update(agent_2_info)
        z_current_2 = agent_2.update(agent_1_info)

        trk_plt.plot_static_obstacle(agent_1_info["x_ob"],\
                                     agent_1_info["y_ob"],\
                                     agent_1_info["phi_ob"],\
                                     2*agent_1_info["l_ob"],\
                                     2*agent_1_info["w_ob"]
                                    )
        trk_plt.plot_static_obstacle(agent_2_info["x_ob"],\
                                     agent_2_info["y_ob"],\
                                     agent_2_info["phi_ob"],\
                                     2*agent_2_info["l_ob"],\
                                     2*agent_2_info["w_ob"]
                                    )
        trk_plt.plot_agents(z_current_1, z_current_2)
        plt.pause(0.01)
        #input("press ENTER")
        trk_plt.clear_agents()
        trk_plt.clear_obstacles()

    ###############################/SIMULATION##################################
    #trk_plt.animate_result(z_data,N,Tf, 'test.gif')
    zinit_vals_1, z_data_full_1 = agent_1.return_sim_data()
    zinit_vals_2, z_data_full_2 = agent_2.return_sim_data()

    agent_1_data = {"zinit": zinit_vals_1,
                    "zdata": z_data_full_1}
    agent_2_data = {"zinit": zinit_vals_2,
                    "zdata": z_data_full_2}

    simdata = {"agent1": agent_1_data,
               "agent2": agent_2_data}

    trk_plt.plot_traj(zinit_vals_1[:,3:])
    plt.show()
    #np.savetxt("full_sol_x_log.csv", )
        #print(simidx)
    return 0

if __name__ == "__main__":

    main()
