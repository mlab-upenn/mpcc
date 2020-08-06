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
from acados_settings import *
from python_sim_utils import   plotter, plot_pajecka, compute_objective
import matplotlib.pyplot as plt
import Bezier
import yaml
import sys


def main_kin():
    np.set_printoptions(precision=3)

    # model parameters
    paramfile = "modelparams.yaml"
    #load global constant model parameters
    with open(paramfile) as file:
        params = yaml.load(file, Loader= yaml.FullLoader)
    lf = params['lf'] #[m]
    lr = params['lr'] #[m]
    lencar = lf + lr

    #sim parameters
    Tsim = 4
    Tf = 1
    N = 20
    Qc = 0.1
    Ql = 100
    Q_theta = 1
    R_d = 0.1
    R_delta = 0.1
    r = 0.2 #trackwidth
    Nsim = np.int(np.floor(N/Tf*Tsim))


    track_lu_table, smax = Bezier.generatelookuptable("tracks/simpleoval")
    trk_plt = plotter(track_lu_table, smax, r, lencar)
    trk_plt.plot_track()

    constraints, model, acados_solver, ocp = acados_settings_kin(Tf, N,  paramfile)

    #starting position in track startidx = theta0[m] * 100 [pts/m]
    startidx = 1050

    trackvars = ['sval', 'tval', 'xtrack', 'ytrack', 'phitrack', 'cos(phi)', 'sin(phi)', 'g_upper', 'g_lower']
    xvars = ['posx', 'posy', 'phi', 'vx', 'theta', 'd', 'delta']
    uvars = ['ddot', 'deltadot', 'thetadot']
    car_soln = []
    xt0 = track_lu_table[startidx,trackvars.index('xtrack')]
    yt0 = track_lu_table[startidx,trackvars.index('ytrack')]
    phit0 = track_lu_table[startidx,trackvars.index('phitrack')]
    theta_hat0 = track_lu_table[startidx,trackvars.index('sval')]
    x0 = np.array([xt0, yt0, phit0, 1, theta_hat0, 0, 0])

    ############################################################################
    #initialization for theta values
    iter = 30
    theta_old = theta_hat0*np.ones((N,))
    x_current = np.tile(x0,(N,1))

    for idx in range(iter):
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
                                r-0.5*lencar
                                ])
            #print("stage idx: ",stageidx,"pval: ",p_val[:-6])
            x_val = x_current[stageidx]
            acados_solver.set(stageidx,"p", p_val)
            acados_solver.set(stageidx,"x", x_val)
        #constrain x0
        acados_solver.set(0, "lbx", x0)
        acados_solver.set(0, "ubx", x0)

        #solve problem
        status = acados_solver.solve()
        #extract theta values
        for idx_sol in range(N):
            xsol = acados_solver.get(idx_sol,"x")
            x_current[idx_sol,:] = xsol

        theta_current = x_current[:,4]

        #compute difference
        theta_diff = np.sum(np.abs(theta_current-theta_old))
        print("theta init difference: ", theta_diff)
        print("theta values", theta_current)
        theta_old = theta_current
    ############################################################################
    #setup using estimated initial trajectory
    theta_vals = theta_current

    step_sol_x_arr = x_current

    print("plotting initialization")
    trk_plt.plot_horizon(theta_vals, step_sol_x_arr[:, :3])
    plt.pause(0.1)
    input("hit [enter] to continue.")
    plt.pause(0.1)
    trk_plt.clear_horizion()


    #list storing visited states
    x0vals = []
    laps = 0
    ##########################SIMULATION#######################################
    for simidx in range(Nsim):

        step_sol_x = []
        step_sol_u = []

        theta_old = theta_vals
        #get track linearization
        index_lin_points = 100 * theta_old - 100*laps*smax
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
                                track_lin_points[stageidx,trackvars.index('sval')] + laps*smax,  #aka theta_hat
                                Qc,
                                Ql,
                                Q_theta,
                                R_d,
                                R_delta,
                                r-0.5*lencar
                                ])
            #print("stage idx: ",stageidx)
            #print("pval: ",p_val[:-5])

            acados_solver.set(stageidx,"p", p_val)
            acados_solver.set(stageidx, "x", step_sol_x_arr[stageidx+1])

        #last stage copy old solution for init
        stageidx = N-1
        p_val = np.array([track_lin_points[stageidx,trackvars.index('xtrack')],
                            track_lin_points[stageidx,trackvars.index('ytrack')],
                            track_lin_points[stageidx,trackvars.index('phitrack')],
                            track_lin_points[stageidx,trackvars.index('sin(phi)')],
                            track_lin_points[stageidx,trackvars.index('cos(phi)')],
                            track_lin_points[stageidx,trackvars.index('sval')] + laps*smax,  #aka theta_hat
                            Qc,
                            Ql,
                            Q_theta,
                            R_d,
                            R_delta,
                            r-0.5*lencar
                            ])
        #print("stage idx: ",stageidx,"pval: ",p_val[:-6])
        acados_solver.set(stageidx,"p", p_val)
        acados_solver.set(stageidx, "x", step_sol_x_arr[stageidx])
        #######################################################################
        #constrain x0
        acados_solver.set(0, "lbx", x0)
        acados_solver.set(0, "ubx", x0)

        status = acados_solver.solve()
        #acados_solver.print_statistics()
        #print(status)

        #extract solution
        x0 = acados_solver.get(1,"x")
        x0vals.append(x0)

        for idx_sol in range(N):
            xsol = acados_solver.get(idx_sol,"x")
            #print("stage:",idx_sol," xsol:", xsol)
            step_sol_x.append(xsol)
            step_sol_u.append(acados_solver.get(idx_sol,"u"))

        step_sol_x_arr = np.array(step_sol_x)
        step_sol_u_arr = np.array(step_sol_u)

        theta_vals = step_sol_x_arr[:, 4]
        print("theta vals", theta_vals)

        objective = compute_objective(Tf/float(N),
                    Qc,
                    Ql,
                    Q_theta,
                    R_d,
                    R_delta,
                    theta_vals,
                    theta_old,
                    step_sol_x_arr[:, :2],
                    step_sol_u_arr,
                    track_lin_points[:,trackvars.index('xtrack'):trackvars.index('ytrack')+1],
                    track_lin_points[:,trackvars.index('phitrack')]
                    )


        #plotting result

        print("objective value", objective)
        trk_plt.plot_horizon(theta_vals, step_sol_x_arr[:, :3])
        trk_plt.plot_input_state_traj(step_sol_x_arr, step_sol_u_arr, xvars, uvars)
        #plt.show()

        plt.pause(0.1)
        input("hit [enter] to continue.")
        plt.pause(0.1)
        trk_plt.clear_horizion()
        trk_plt.clear_input_state_traj()

        #preparation for next timestep
        theta_vals = np.hstack((step_sol_x_arr[1:, 4], step_sol_x_arr[-1, 4]+0.1))

        if theta_vals[0] > (laps+1)*smax :
            print("#################################RESET###############################")
            laps = laps + 1


    ###############################/SIMULATION##################################
    trk_plt.plot_traj(np.array(x0vals))
    plt.show()
    #np.savetxt("full_sol_x_log.csv", )
        #print(simidx)
    return 0

def main_dyn():
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
    lencar = lf+lr

    #sim parameters
    Tsim = 20
    Tf = 1
    N = 20
    Qc = 0.1
    Ql = 1000
    Q_theta = 10
    R_d = 0.01
    R_delta = 0.01
    Nsim = np.int(np.floor(N/Tf*Tsim))
    r = 0.2 #trackwidth

    track_lu_table, smax = Bezier.generatelookuptable("tracks/simpleoval")
    trk_plt = plotter(track_lu_table, smax, r, lencar)
    trk_plt.plot_track()

    constraints, model, acados_solver, ocp = acados_settings_dyn(Tf, N,  paramfile)
    #plot_pajecka(paramfile)

    #starting position in track startidx = theta0[m] * 100 [pts/m]
    startidx = 20

    trackvars = ['sval', 'tval', 'xtrack', 'ytrack', 'phitrack', 'cos(phi)', 'sin(phi)', 'g_upper', 'g_lower']
    xvars = ['posx', 'posy', 'phi', 'vx', 'vy', 'omega', 'd', 'delta', 'theta']
    uvars = ['ddot', 'deltadot', 'thetadot']
    pvars = ['xt', 'yt', 'phit', 'sin_phit', 'cos_phit', 'theta_hat', 'Qc', 'Ql', 'Q_theta', 'R_d', 'R_delta', 'r']

    xt0 = track_lu_table[startidx,trackvars.index('xtrack')]
    yt0 = track_lu_table[startidx,trackvars.index('ytrack')]
    phit0 = track_lu_table[startidx,trackvars.index('phitrack')]
    theta_hat0 = track_lu_table[startidx,trackvars.index('sval')]
    xinit = np.array([xt0, yt0, phit0, 1, 0.01, 0, theta_hat0, 0, 0])

    ############################################################################
    #initialization for theta values
    iter = 3
    x_current = np.tile(xinit,(N,1))
    #arbitrarily set theta  values and
    theta_old = theta_hat0*np.ones((N,)) + 0.001*np.arange(N)
    x_current[:,xvars.index('theta')] = theta_old
    index_lin_points = 100 * theta_old
    index_lin_points = index_lin_points.astype(np.int32)
    track_lin_points = track_lu_table[index_lin_points,:]

    #initialize x values on track
    x_current[:,0] = track_lin_points[:,trackvars.index('xtrack')]
    x_current[:,1] = track_lin_points[:,trackvars.index('ytrack')]
    x_current[:,2] = track_lin_points[:,trackvars.index('phitrack')]

    for idx in range(iter):
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
                                r
                                ])
            #print("stage idx: ",stageidx,"pval: ",p_val[:-6])
            x_val = x_current[stageidx]
            acados_solver.set(stageidx,"p", p_val)
            acados_solver.set(stageidx,"x", x_val)

        #constrain x0
        acados_solver.set(0, "lbx", xinit)
        acados_solver.set(0, "ubx", xinit)

        #solve problem
        status = acados_solver.solve()
        #extract theta values
        for idx_sol in range(N):
            xsol = acados_solver.get(idx_sol,"x")
            x_current[idx_sol,:] = xsol

        theta_current = x_current[:,xvars.index('theta')]

        #compute difference
        theta_diff = np.sum(np.abs(theta_current-theta_old))
        print("theta init difference: ", theta_diff)
        print("theta values", theta_current)
        theta_old = theta_current
    ############################################################################
    #setup using estimated initial trajectory
    theta_vals = theta_current

    step_sol_x_arr = x_current

    print("plotting initialization")
    trk_plt.plot_horizon(theta_vals, step_sol_x_arr[:, :3])
    plt.pause(0.1)
    input("hit [enter] to continue.")
    plt.pause(0.1)
    trk_plt.clear_horizion()


    #list storing visited states
    x0vals = []
    laps = 0
    ##########################SIMULATION#######################################
    for simidx in range(Nsim):
        step_sol_x = []
        step_sol_u = []

        theta_old = theta_vals
        #get track linearization
        index_lin_points = 100 * theta_old - 100*laps*smax
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
                                track_lin_points[stageidx,trackvars.index('sval')]+ laps*smax,  #aka theta_hat
                                Qc,
                                Ql,
                                Q_theta,
                                R_d,
                                R_delta,
                                r
                                ])
            #print("stage idx: ",stageidx,"pval: ",p_val[:-6])
            acados_solver.set(stageidx,"p", p_val)
            acados_solver.set(stageidx, "x", step_sol_x_arr[stageidx+1])

        #last stage copy old solution for init
        stageidx = N-1
        p_val = np.array([track_lin_points[stageidx,trackvars.index('xtrack')],
                            track_lin_points[stageidx,trackvars.index('ytrack')],
                            track_lin_points[stageidx,trackvars.index('phitrack')],
                            track_lin_points[stageidx,trackvars.index('sin(phi)')],
                            track_lin_points[stageidx,trackvars.index('cos(phi)')],
                            track_lin_points[stageidx,trackvars.index('sval')]+ laps*smax,  #aka theta_hat
                            Qc,
                            Ql,
                            Q_theta,
                            R_d,
                            R_delta,
                            r
                            ])
        #print("stage idx: ",stageidx,"pval: ",p_val[:-6])
        acados_solver.set(stageidx,"p", p_val)
        acados_solver.set(stageidx, "x", step_sol_x_arr[stageidx])
        #######################################################################
        #constrain x0
        acados_solver.set(0, "lbx", x0)
        acados_solver.set(0, "ubx", x0)

        status = acados_solver.solve()
        #acados_solver.print_statistics()
        #print(status)

        #extract solution
        x0 = acados_solver.get(1,"x")
        x0vals.append(x0)

        for idx_sol in range(N):
            xsol = acados_solver.get(idx_sol,"x")
            #print("stage:",idx_sol," xsol:", xsol)
            step_sol_x.append(xsol)
            step_sol_u.append(acados_solver.get(idx_sol,"u"))

        step_sol_x_arr = np.array(step_sol_x)
        step_sol_u_arr = np.array(step_sol_u)

        theta_vals = step_sol_x_arr[:, 6]
        print("theta vals", theta_vals)

        objective = compute_objective(Tf/float(N),
                    Qc,
                    Ql,
                    Q_theta,
                    R_d,
                    R_delta,
                    theta_vals,
                    theta_old,
                    step_sol_x_arr[:, :2],
                    step_sol_u_arr,
                    track_lin_points[:,trackvars.index('xtrack'):vars.index('ytrack')+1],
                    track_lin_points[:,trackvars.index('phitrack')]
                    )


        #plotting result

        print("objective value", objective)
        trk_plt.plot_horizon(theta_vals, step_sol_x_arr[:, :2])
        trk_plt.plot_input_state_traj(step_sol_x_arr, step_sol_u_arr, xvars, uvars)
        #plt.show()

        plt.pause(0.1)
        input("hit [enter] to continue.")
        plt.pause(0.1)
        trk_plt.clear_horizion()
        trk_plt.clear_input_state_traj()

        #preparation for next timestep
        theta_vals = np.hstack((step_sol_x_arr[1:, 4], step_sol_x_arr[-1, 4]+0.1))

        if theta_vals[0] > smax:
            laps = laps + 1
            print("#################################RESET###############################")

    ###############################/SIMULATION##################################
    trk_plt.plot_traj(np.array(x0vals))
    plt.show()
    #np.savetxt("full_sol_x_log.csv", )
        #print(simidx)
    return 0

if __name__ == "__main__":
    arg = sys.argv[1]

    if arg == "kin":
        print("executing with kinematic_model")
        main_kin()
    elif arg == "dyn":
        print("executing with dynamic_model")
        main_dyn()
    else:
        print("Please run model with valid modeltype: [kin, dyn]")
        print("e.g. python3 python_sim.py kin")
