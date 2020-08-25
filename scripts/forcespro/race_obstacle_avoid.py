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
from python_sim_utils import   plotter, plot_pajecka
import matplotlib.pyplot as plt
import InterpolateTrack
import yaml
import pickle
import sys
from datetime import datetime
from pathlib import Path

def main(generatesolvers):
    np.set_printoptions(precision=6)
    np.set_printoptions(threshold=sys.maxsize)

    #state variable ordering
    trackvars = ['sval', 'tval', 'xtrack', 'ytrack', 'phitrack', 'cos(phi)', 'sin(phi)', 'g_upper', 'g_lower']
    zvars = ['ddot', 'deltadot', 'thetadot', 'posx', 'posy', 'phi', 'vx', 'vy', 'omega', 'd', 'delta', 'theta']
    xvars = ['posx', 'posy', 'phi', 'vx', 'vy', 'omega', 'd', 'delta', 'theta']
    debug_colors = ['b','g','r','c','m','y']

    # parameters for both agents
    modelparams_1 = "parameters/modelparams_1.yaml"
    solverparams_1 = "parameters/solverparams_1.yaml"
    modelparams_2 = "parameters/modelparams_2.yaml"
    solverparams_2 = "parameters/solverparams_2.yaml"

    #sim parameters
    with open(solverparams_1) as file:
        params = yaml.load(file, Loader= yaml.FullLoader)
    Tsim = 15
    Tf = params['Tf']
    N = params['N']
    Nsim = np.int(np.floor(N/Tf*Tsim))

    trackname = "slider"
    track_lu_table, smax = InterpolateTrack.generatelookuptable("tracks/"+trackname)
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
    agent_1 = racer(track, Tsim, "agent_1", modelparams_1, solverparams_1, generatesolvers)
    agent_1_info = {"const_dactive": 0,
                    "phi_ob": 0,
                    "x_ob": 0,
                    "y_ob": 0,
                    "l_ob": lencar*1.2,
                    "w_ob": widthcar*1.2}

    agent_2 = racer(track, Tsim, "agent_2", modelparams_2, solverparams_2, generatesolvers)
    agent_2_info = {"const_dactive": 0,
                    "phi_ob": 0,
                    "x_ob": 0,
                    "y_ob": 0,
                    "l_ob": lencar*1.2,
                    "w_ob": widthcar*1.2}

    #agent 1 trajecotry initialization
    #starting position in track startidx = theta0[m] * 100 [pts/m]
    startidx_1 = 400
    xt0 = track_lu_table[startidx_1,trackvars.index('xtrack')]
    yt0 = track_lu_table[startidx_1,trackvars.index('ytrack')]
    phit0 = track_lu_table[startidx_1,trackvars.index('phitrack')]
    theta_hat0 = track_lu_table[startidx_1,trackvars.index('sval')]
    #initial condition
    xinit_1 = np.array([xt0, yt0, phit0, 0.2, 0.0, 0, 0, 0, theta_hat0])
    agent_1_info["phi_ob"] = np.tile(xinit_1[xvars.index('phi')],(N,))
    agent_1_info["x_ob"] = np.tile(xinit_1[xvars.index('posx')],(N,))
    agent_1_info["y_ob"] = np.tile(xinit_1[xvars.index('posy')],(N,))

    #agent 2  trajectory initialization
    startidx_2 = 600
    xt0 = track_lu_table[startidx_2,trackvars.index('xtrack')]
    yt0 = track_lu_table[startidx_2,trackvars.index('ytrack')]
    phit0 = track_lu_table[startidx_2,trackvars.index('phitrack')]
    theta_hat0 = track_lu_table[startidx_2,trackvars.index('sval')]
    #initial condition
    xinit_2 = np.array([xt0, yt0, phit0, 0.2, 0.0, 0, 0, 0, theta_hat0])
    agent_2_info["phi_ob"] = np.tile(xinit_2[xvars.index('phi')],(N,))
    agent_2_info["x_ob"] = np.tile(xinit_2[xvars.index('posx')],(N,))
    agent_2_info["y_ob"] = np.tile(xinit_2[xvars.index('posy')],(N,))

    #compute trajectory init
    z_current_1 = agent_1.initialize_trajectory(xinit_1, agent_2_info, startidx_1)
    z_current_2 = agent_2.initialize_trajectory(xinit_2, agent_1_info, startidx_2)

    trk_plt.plot_agents(z_current_1, z_current_2)
    trk_plt.plot_static_obstacle(agent_1_info["x_ob"][0],\
                                 agent_1_info["y_ob"][0],\
                                 agent_1_info["phi_ob"][0],\
                                 2*agent_1_info["l_ob"],\
                                 2*agent_1_info["w_ob"],\
                                 debug_colors[np.mod(0,7)]
                                )
    trk_plt.plot_static_obstacle(agent_2_info["x_ob"][0],\
                                 agent_2_info["y_ob"][0],\
                                 agent_2_info["phi_ob"][0],\
                                 2*agent_2_info["l_ob"],\
                                 2*agent_2_info["w_ob"],\
                                 debug_colors[np.mod(0,7)]
                                )
    plt.pause(0.1)
    input("press ENTER")
    trk_plt.clear_agents()
    trk_plt.clear_obstacles()
    ##########################SIMULATION#######################################
    for simidx in range(Nsim):
        agent_1_info["phi_ob"] = z_current_1[:, zvars.index('phi')]
        agent_1_info["x_ob"] = z_current_1[:, zvars.index('posx')]
        agent_1_info["y_ob"] = z_current_1[:, zvars.index('posy')]
        agent_2_info["phi_ob"] = z_current_2[:, zvars.index('phi')]
        agent_2_info["x_ob"] = z_current_2[:, zvars.index('posx')]
        agent_2_info["y_ob"] = z_current_2[:, zvars.index('posy')]
        #activate constraints
        agent_1_info["const_dactive"] = 0
        agent_2_info["const_dactive"] = 0
        '''
        dist = np.sqrt(np.sum(np.square(z_current_1[0, zvars.index('posx'):zvars.index('posy')+1]\
                -z_current_2[0, zvars.index('posx'):zvars.index('posy')+1])))
        print("[INFO] Distance",dist)
        if dist>0.5:
            #deactivate constraints
            agent_1_info["const_dactive"] = 1
            agent_2_info["const_dactive"] = 1
        else:
        '''
        z_current_1 = agent_1.update(agent_2_info)
        z_current_2 = agent_2.update(agent_1_info)


        '''
        for stageidx in range(N):
            trk_plt.plot_static_obstacle(agent_1_info["x_ob"][stageidx],\
                                         agent_1_info["y_ob"][stageidx],\
                                         agent_1_info["phi_ob"][stageidx],\
                                         agent_1_info["l_ob"],\
                                         agent_1_info["w_ob"],\
                                         debug_colors[np.mod(stageidx,6)]
                                        )
            trk_plt.plot_static_obstacle(agent_2_info["x_ob"][stageidx],\
                                         agent_2_info["y_ob"][stageidx],\
                                         agent_2_info["phi_ob"][stageidx],\
                                         agent_2_info["l_ob"],\
                                         agent_2_info["w_ob"],\
                                         debug_colors[np.mod(stageidx,6)]
                                        )

        trk_plt.plot_agents(z_current_1, z_current_2)

        plt.pause(0.01)
        input("press ENTER")
        trk_plt.clear_agents()
        trk_plt.clear_obstacles()
        '''
    ###############################/SIMULATION##################################
    #trk_plt.animate_result(z_data,N,Tf, 'test.gif')
    zinit_vals_1, z_data_full_1 = agent_1.return_sim_data()
    zinit_vals_2, z_data_full_2 = agent_2.return_sim_data()

    agent_1_data = {"zinit": zinit_vals_1,
                    "zdata": z_data_full_1}
    agent_2_data = {"zinit": zinit_vals_2,
                    "zdata": z_data_full_2}

    simdata = {"agent1": agent_1_data,
               "agent2": agent_2_data,
               "smax" :  smax,
               "lencar": lencar,
               "track" : trackname,
               "r" : r}

    now = datetime.now()
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%m_%d_%H_%M_%S")
    name_data = "2_agent_race_"+dt_string+".pkl"
    root = Path(".")
    full_path = root / "simdata" / name_data
    with open(full_path,'wb') as f:
        pickle.dump(simdata, f)

    trk_plt.plot_traj(zinit_vals_1[:,3:])
    plt.show()
    #np.savetxt("full_sol_x_log.csv", )
        #print(simidx)
    return 0

if __name__ == "__main__":
    options = sys.argv[1]
    print("[INFO] Use the argument gen to generate solvers newly")
    if options == "gen":
        print("[INFO] REGENERATING SOLVERS")
        main(1)
    else :
        print("[INFO] LOADING OLD SOLVERS")
        main(0)
