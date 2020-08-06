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
from generate_solver_dynamic import get_forces_solver_dynamic
import forcespro.nlp
from python_sim_utils import   plotter, plot_pajecka, compute_objective
import matplotlib.pyplot as plt
import Bezier
import yaml
import sys


def main_dyn():
    np.set_printoptions(precision=4)
    np.set_printoptions(threshold=sys.maxsize)
    # model parameters
    paramfile = "modelparams.yaml"
    #load global constant model parameters
    with open(paramfile) as file:
        params = yaml.load(file, Loader= yaml.FullLoader)
    lf = params['lf'] #[m]
    lr = params['lr'] #[m]
    #pajecka and motor coefficients
    Bf = params['Bf']
    Br = params['Br']
    Cf = params['Cf']
    Cr = params['Cr']
    Cm1 = params['Cm1']
    Cm2 = params['Cm2']
    Croll = params['Croll']
    Cd = params['Cd']
    Df = params['Df']
    Dr = params['Dr']
    lencar = 2*(lf+lr)

    #sim parameters
    Tsim = 10
    Tf = 1.5
    N = 30
    Qc = 0.1
    Ql = 1000
    Q_theta = 100
    R_d = 0.01
    R_delta = 0.01
    Nsim = np.int(np.floor(N/Tf*Tsim))
    r = 0.15 #trackwidth

    solver = get_forces_solver_dynamic(N, Tf, paramfile)

    track_lu_table, smax = Bezier.generatelookuptable("tracks/sample_track")
    trk_plt = plotter(track_lu_table, smax, r, lencar)
    #plot_pajecka(paramfile)
    trk_plt.plot_track()

    #starting position in track startidx = theta0[m] * 100 [pts/m]
    startidx = 3000

    trackvars = ['sval', 'tval', 'xtrack', 'ytrack', 'phitrack', 'cos(phi)', 'sin(phi)', 'g_upper', 'g_lower']
    xvars = ['posx', 'posy', 'phi', 'vx', 'vy', 'omega', 'd', 'delta', 'theta']
    uvars = ['ddot', 'deltadot', 'thetadot']
    pvars = ['xt', 'yt', 'phit', 'sin_phit', 'cos_phit', 'theta_hat', 'Qc', 'Ql', 'Q_theta', 'R_d', 'R_delta', 'r']
    car_soln = []
    xt0 = track_lu_table[startidx,trackvars.index('xtrack')]
    yt0 = track_lu_table[startidx,trackvars.index('ytrack')]
    phit0 = track_lu_table[startidx,trackvars.index('phitrack')]
    theta_hat0 = track_lu_table[startidx,trackvars.index('sval')]

    #initial condition
    zvars = ['ddot', 'deltadot', 'thetadot', 'posx', 'posy', 'phi', 'vx', 'vy', 'omega', 'd', 'delta', 'theta']
    xinit = np.array([xt0, yt0, phit0, 1.8, 0.0, 0, 0, 0, theta_hat0])
    zinit = np.concatenate([np.array([0,0,0]), xinit])

    ############################################################################
    #initialization for theta values
    iter = 100
    z_current = np.tile(zinit,(N,1))
    #arbitrarily set theta  values and
    theta_old = theta_hat0*np.ones((N,)) + 0.001*np.arange(N)
    z_current[:,11] = theta_old
    index_lin_points = 100 * theta_old
    index_lin_points = index_lin_points.astype(np.int32)
    track_lin_points = track_lu_table[index_lin_points,:]

    #initialize x values on track
    z_current[:,3] = track_lin_points[:,trackvars.index('xtrack')]
    z_current[:,4] = track_lin_points[:,trackvars.index('ytrack')]
    z_current[:,5] = track_lin_points[:,trackvars.index('phitrack')]

    for idx in range(iter):

        all_parameters = []
        #get track linearization
        index_lin_points = 100 * theta_old
        index_lin_points = index_lin_points.astype(np.int32)
        track_lin_points = track_lu_table[index_lin_points,:]

        for stageidx in range(N):
            p_val = np.array([track_lin_points[stageidx,trackvars.index('xtrack')],
                                track_lin_points[stageidx,trackvars.index('ytrack')],
                                track_lin_points[stageidx,trackvars.index('phitrack')],
                                track_lin_points[stageidx,trackvars.index('sin(phi)')],
                                track_lin_points[stageidx,trackvars.index('cos(phi)')],
                                track_lin_points[stageidx,trackvars.index('sval')],  #aka theta_hat
                                Qc,
                                Ql,
                                Q_theta,
                                R_d,
                                R_delta,
                                r-lencar/2
                                ])
            all_parameters.append(p_val)

        all_parameters = np.array(all_parameters)

        #problem dictionary, arrays have to be flattened
        problem = {"x0": z_current.reshape(-1,),
                   "xinit": xinit,
                   "all_parameters": all_parameters.reshape(-1,)}
        #solve problem
        output, exitflag, info = solver.solve(problem)
        #print(info)
        #print(output)
        #input("hit [enter] to continue.")

        #extract theta values
        idx_sol = 0
        for key in output:
            #print(key)
            zsol = output[key]
            usol = zsol[0:3]
            z_current[idx_sol,:] = zsol
            idx_sol = idx_sol+1

        theta_current = z_current[:, zvars.index('theta')]

        #compute difference
        theta_diff = np.sum(np.abs(theta_current-theta_old))
        print("theta init difference: ", theta_diff)
        print("theta values", theta_current)
        theta_old = theta_current

    ############################################################################
    #setup using estimated initial trajectory
    theta_vals = theta_current

    step_sol_z_arr = z_current
    print("resulting z:")
    print(step_sol_z_arr)

    print("plotting initialization")
    trk_plt.plot_horizon(theta_vals, step_sol_z_arr[:, 3:6])
    trk_plt.plot_input_state_traj(step_sol_z_arr, zvars)
    plt.pause(0.1)
    input("hit [enter] to continue.")
    plt.pause(0.1)
    trk_plt.clear_horizion()
    trk_plt.clear_input_state_traj()


    #list storing visited states
    zinit_vals = np.zeros((Nsim, 12))
    z_data = np.zeros((Nsim,N,12))
    laps = 0
    ##########################SIMULATION#######################################
    for simidx in range(Nsim):
        step_sol_u = []
        all_parameters = []
        theta_old = theta_vals
        #get track linearization
        index_lin_points = 100 * theta_old
        index_lin_points = index_lin_points.astype(np.int32)
        print("track linearized around entries:", index_lin_points )
        track_lin_points = track_lu_table[index_lin_points,:]

        #######################################################################
        #set params and warmstart
        for stageidx in range(N-1):
            p_val = np.array([track_lin_points[stageidx,trackvars.index('xtrack')],
                                track_lin_points[stageidx,trackvars.index('ytrack')],
                                track_lin_points[stageidx,trackvars.index('phitrack')],
                                track_lin_points[stageidx,trackvars.index('sin(phi)')],
                                track_lin_points[stageidx,trackvars.index('cos(phi)')],
                                track_lin_points[stageidx,trackvars.index('sval')],  #aka theta_hat
                                Qc,
                                Ql,
                                Q_theta,
                                R_d,
                                R_delta,
                                r-lencar/2
                                ])
            #create parameter matrix
            all_parameters.append(p_val)
            #stack state initializations
            z_current[stageidx,:] = step_sol_z_arr[stageidx+1,:]

        #last stage copy old solution for init
        stageidx = N-1
        p_val = np.array([track_lin_points[stageidx,trackvars.index('xtrack')],
                            track_lin_points[stageidx,trackvars.index('ytrack')],
                            track_lin_points[stageidx,trackvars.index('phitrack')],
                            track_lin_points[stageidx,trackvars.index('sin(phi)')],
                            track_lin_points[stageidx,trackvars.index('cos(phi)')],
                            track_lin_points[stageidx,trackvars.index('sval')],  #aka theta_hat
                            Qc,
                            Ql,
                            Q_theta,
                            R_d,
                            R_delta,
                            r-lencar/2
                            ])
        all_parameters.append(p_val)
        all_parameters = np.array(all_parameters)
        z_current[N-1,:] = step_sol_z_arr[N-1,:]

        #######################################################################
        #problem dictionary, arrays have to be flattened
        problem = {"x0": z_current.reshape(-1,),
                   "xinit": xinit,
                   "all_parameters": all_parameters.reshape(-1,)}
        #solve problem
        output, exitflag, info = solver.solve(problem)

        #extract solution
        idx_sol = 0
        for key in output:
            zsol = output[key]
            usol = zsol[0:3]
            step_sol_u.append(usol)
            step_sol_z_arr[idx_sol,:] = zsol
            idx_sol = idx_sol+1

        #log solution
        z_data[simidx,:,:] = step_sol_z_arr

        print("theta: ", zinit[zvars.index('theta')])
        print("vx: ", zinit[zvars.index('vx')])
        print("vy: ", zinit[zvars.index('vy')])
        print("omega: ", zinit[zvars.index('omega')])
        print("phi: ", zinit[zvars.index('phi')]*180/3.1415)
        print("d: ", zinit[zvars.index('d')])
        print("delta: ", zinit[zvars.index('delta')])

        vx = zinit[zvars.index('vx')]
        vy = zinit[zvars.index('vy')]
        delta = zinit[zvars.index('delta')]
        omega = zinit[zvars.index('omega')]
        d = zinit[zvars.index('d')]
        alphaf = -np.arctan2((omega*lf + vy), vx) + delta
        Ffy = Df*np.sin(Cf*np.arctan(Bf*alphaf))
        alphar = np.arctan2((omega*lr - vy),vx)
        Fry = Dr*np.sin(Cr*np.arctan(Br*alphar))
        Frx = (Cm1-Cm2*vx) * d - Cr -Cd*vx*vx

        print("Ffy: ", Ffy)
        print("Fry: ", Fry)
        print("Frx: ", Frx)

        theta_vals = step_sol_z_arr[:,zvars.index('theta')]
        zinit = step_sol_z_arr[0,:]
        xinit = zinit[3:]

        zinit_vals[simidx,:] = zinit
        step_sol_u_arr = np.array(step_sol_u)


        #plotting result
        trk_plt.plot_horizon(theta_vals, step_sol_z_arr[:, 3:6])
        trk_plt.plot_input_state_traj(step_sol_z_arr, zvars)
        #plt.show()

        plt.pause(0.1)
        input("hit [enter] to continue.")
        plt.pause(0.1)
        trk_plt.clear_horizion()
        trk_plt.clear_input_state_traj()

        #preparation for next timestep
        theta_vals = np.hstack((step_sol_z_arr[1:, zvars.index('theta')], step_sol_z_arr[-1, zvars.index('theta')]+0.1))

        if theta_vals[0] > smax :
            print("#################################RESET###############################")
            laps = laps + 1
            print("lap:", laps)
            theta_vals = theta_vals - smax
            step_sol_z_arr[:,zvars.index('theta')] = theta_vals

    ###############################/SIMULATION##################################
    #trk_plt.animate_result(z_data,N,Tf, 'test.gif')
    trk_plt.plot_traj(zinit_vals[:,3:])
    plt.show()
    #np.savetxt("full_sol_x_log.csv", )
        #print(simidx)
    return 0

if __name__ == "__main__":

    main_dyn()
